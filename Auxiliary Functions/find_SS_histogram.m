

function [residual,mHistogram_out,mAssetsPrime_out,mConsumption_out] = find_SS_histogram(x)

    %% Initialization (global variables)
    global ggamma aalpha ddelta mmu llambda AA aggEmployment tau_n ...
    	   mZTransition mEpsilonTransition vZGrid  ...
           nEpsilon nZ nAssetsFine vAssetsGridFine mAssetsGridFine mZGridFine

    %% Declare prices
        r = x(1); 
        p = x(2);
        X = (p*aggEmployment^(aalpha-1)/(AA*aalpha))^(1/(aalpha-mmu));
        w = (1-aalpha)*AA*(X/aggEmployment)^aalpha/(1-tau_n) ;
        
        
disp('Interest Rate     Optimal Capital')
disp([r (mmu*p*vZGrid.^mmu/(r+ddelta)).^(1/(1-mmu))])
        
    
    %% Agent's policy functions (by EGM)
        [AssetsPrime_EGM,Consumption_EGM] = find_PF_EGM(r,p,w) ;
    
    %% Stationary distribution (Histogram)
        % Rename sizes
            n_a = nAssetsFine ; 
            n_e = nEpsilon    ;
            n_z = nZ          ;
    
        % Vectorize saving function
        vS = AssetsPrime_EGM(:) ; 
    
        % Compute transition matrix associated with policy rules
        % Find index and value of nearest neighbors of Savings in A_grid
        vInd = knnsearch(vAssetsGridFine,vS);
        vVal = vAssetsGridFine(vInd);
        
        vIndBelow = max( vInd - (vVal > vS) , 1   ) ; vValBelow = vAssetsGridFine(vIndBelow);
        vIndAbove = min( vIndBelow + 1      , n_a ) ; vValAbove = vAssetsGridFine(vIndAbove);

        % Compute lottery matrices
        vPrBelow = min( max( (vValAbove - vS) ./ (vValAbove - vValBelow) , 0 ) , 1 );
        vPrAbove = min( max( (vS - vValBelow) ./ (vValAbove - vValBelow) , 0 ) , 1 );

        % Assets transition matrix
        mTransitionAbove = zeros(n_z*n_e*n_a,n_a);
        mTransitionBelow = zeros(n_z*n_e*n_a,n_a);
        for i_a = 1:n_a
            mTransitionBelow(vIndBelow == i_a,i_a) = vPrBelow(vIndBelow == i_a);
            mTransitionAbove(vIndAbove == i_a,i_a) = vPrAbove(vIndAbove == i_a);
        end
        mTransition = mTransitionBelow + mTransitionAbove ; 
        
        % mAssetsTransition = zeros(n_a*n_e*n_z,n_a*n_e*n_z);
        % for i = 1:n_z*n_e*n_a
        %     mAssetsTransition(i,:) = reshape(mTransition(i,:)'*ones(1,n_e*n_z),1,n_z*n_e*n_a);
        % end
        mAssetsTransition = repmat(mTransition, [1,n_e*n_z]) ;

        % Epsilon
        mEpsilonTransitionHistogram = repmat( kron(mEpsilonTransition,ones(n_a)) , [n_z,n_z] );
        
        % Z
        mZTransitionHistogram = kron( ggamma*eye(n_z) + (1-ggamma)*mZTransition , ones(n_a*n_e) ) ;

        % Full transition matrix
        mTransition = mAssetsTransition .* mEpsilonTransitionHistogram .* mZTransitionHistogram ;

        % Compute invariant histogram
        errHistogram = 100;	iterationHistogram = 0;
        G = ones(n_a*n_e*n_z,1) ./ (n_a*n_e*n_z);
        while errHistogram > 1e-12 && iterationHistogram < 1e4
            G_New = mTransition' * G;
            errHistogram = max(abs(G_New - G)); 
            iterationHistogram = iterationHistogram + 1;
            G = G_New;
        end

        % Expand histogram matrix
        mHistogram = reshape(full(G),n_a,n_e,n_z);    
        
        
    %% Update prices 
        A = sum(vAssetsGridFine.*sum(sum(mHistogram,3),2)) ;
        K = min( llambda*mAssetsGridFine , (mmu*p*mZGridFine.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
        X = sum(sum(sum( (mZGridFine.*K).^mmu .* mHistogram ))).^(1/mmu) ;
        % w_new = (1-aalpha)*AA*(X/aggEmployment)^aalpha/(1-tau_n) ;
        p_new = (aalpha)*AA*(X)^(aalpha-mmu)*aggEmployment^(1-aalpha) ;
        
        f = @(R) ( sum(sum(sum( min( llambda*mAssetsGridFine , (mmu*p_new*mZGridFine.^mmu/(R+ddelta)).^(1/(1-mmu)) ).*mHistogram )))/A-1 ) ;
        options = optimoptions('fsolve','Display','off');
        [r_new,err,exitflag] = fsolve(f, r ,options);
                        
    
    % Residual
        residual = [r_new/r-1 ; p_new/p-1] ;
    
    %% Optional output
        if nargout>1
            mHistogram_out = mHistogram ;
            if nargout>2
                mAssetsPrime_out = AssetsPrime_EGM ;
                mConsumption_out = Consumption_EGM ;
            end 
        end 


end