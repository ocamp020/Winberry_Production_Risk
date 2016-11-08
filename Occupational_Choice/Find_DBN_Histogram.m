% Production Risk Economy with Occupational Choice
% Find Stationary Distribution with the Histogram Method
% Juan David Herreno & Sergio Ocampo (2016)

function [residual,mDBN_W_out,mDBN_E_out,mAp_W_out,mAp_E_out,OC_W_out,OC_E_out,VW_out,VE_out] = Find_DBN_Histogram(x)


%% Initialization (global variables)
    global ggamma aalpha ddelta mmu llambda AA tau_n mKappa ...
    	   mZ_Transition mE_Transition_W mE_Transition_E vZ_Grid vE_Grid  ...
           n_E n_Z n_A vA_Grid mA_Grid mZ_Grid n_State

%% Declare prices
    r = x(1); 
    p = x(2);
    w = x(3);
        
        
% disp('Price Wage Interest Rate     Optimal Capital')
% disp([p w r (mmu*p*vZ_Grid.^mmu/(r+ddelta)).^(1/(1-mmu))])
        
    
%% Agent's value and policy functions (by discrete VFI)
    [Ap_W_VFI,Ap_E_VFI,OC_W_VFI,OC_E_VFI,VW_VFI,VE_VFI] = Discrete_VFI(r,p,w) ;
    % disp(VW_VFI(:)) 
%     figure; plot(vA_Grid,[VW_VFI(:,:,1),VW_VFI(:,:,2),VW_VFI(:,:,3)]); title('VW');
%     figure; plot(vA_Grid,[VE_VFI(:,:,1),VE_VFI(:,:,2),VE_VFI(:,:,3)]); title('VE');
%     figure; plot(vA_Grid,[Ap_W_VFI(:,:,1),Ap_W_VFI(:,:,2),Ap_W_VFI(:,:,3)]); title('Ap_W');
%     figure; plot(vA_Grid,[Ap_E_VFI(:,:,1),Ap_E_VFI(:,:,2),Ap_E_VFI(:,:,3)]); title('Ap_E');
%     figure; plot(vA_Grid,[OC_W_VFI(:,:,1),OC_W_VFI(:,:,2),OC_W_VFI(:,:,3)]); title('Ap_W');
%     figure; plot(vA_Grid,[OC_E_VFI(:,:,1),OC_E_VFI(:,:,2),OC_E_VFI(:,:,3)]); title('Ap_E');
    
    
%% Stationary distribution (Histogram)
    
    % Vectorize saving function for both types of agent (W,E)
    vS_W = Ap_W_VFI(:) ; % In Index
    vS_E = Ap_E_VFI(:) ; % In Index
    
    % Assets transition matrix (Before occupational choice)
    mA_W_Transition = zeros(n_Z*n_E*n_A,n_A);
    mA_E_Transition = zeros(n_Z*n_E*n_A,n_A);
    for i_a = 1:n_A
        mA_W_Transition(vS_W == i_a,i_a) = 1 ;
        mA_E_Transition(vS_E == i_a,i_a) = 1 ;
    end
    mA_Transition = [repmat(mA_W_Transition, [1,n_E*n_Z]) , zeros(n_State) ;
                     zeros(n_State) , repmat(mA_E_Transition, [1,n_E*n_Z])];
    mA_Transition = sparse(mA_Transition) ;
    
    % Epsilon
    mE_Transition_Histogram = [ repmat( kron(mE_Transition_W,ones(n_A)) , [n_Z,n_Z] ) , zeros(n_State) ; 
                                zeros(n_State) , repmat( kron(mE_Transition_E,ones(n_A)) , [n_Z,n_Z] ) ];
    mE_Transition_Histogram = sparse(mE_Transition_Histogram) ;

    % Z
    mZ_Transition_Histogram = kron( eye(2) , kron( ggamma*eye(n_Z) + (1-ggamma)*mZ_Transition , ones(n_A*n_E) ) ) ;
    mZ_Transition_Histogram = sparse(mZ_Transition_Histogram) ;
    
    % Full transition matrix (Before occupational choice)
    mTransition = mA_Transition .* mE_Transition_Histogram .* mZ_Transition_Histogram ;
    
    % Ocupational Choice
    OC_W = OC_W_VFI(:) ; 
    OC_E = OC_E_VFI(:) ; 
    
        % Adjust savings choice by switching cost
        vS_W_cost = vA_Grid(vS_W)-mKappa(:) ;

        % Find index and value of nearest neighbors of adjusted savings in A_grid
        vInd = knnsearch(vA_Grid,vS_W_cost);
        vVal = vA_Grid(vInd);

        vIndBelow = max( vInd - (vVal > vS_W_cost) , 1   ) ; vValBelow = vA_Grid(vIndBelow);
        vIndAbove = min( vIndBelow + 1             , n_A ) ; vValAbove = vA_Grid(vIndAbove);

        % Compute lottery matrices
        vPrBelow = min( max( (vValAbove - vS_W_cost) ./ (vValAbove - vValBelow) , 0 ) , 1 );
        vPrAbove = min( max( (vS_W_cost - vValBelow) ./ (vValAbove - vValBelow) , 0 ) , 1 );
        
        % Adjust mTransition for type W agents
        aux = [1:n_State];
        aux = aux(OC_W==1);
        for i = aux
            mTransition(1:n_State,n_State + vIndBelow(i)) = mTransition(1:n_State,n_State + vIndBelow(i)) + ...
                                                vPrBelow(i)*mTransition(1:n_State,i) ;
            mTransition(1:n_State,n_State + vIndAbove(i)) = mTransition(1:n_State,n_State + vIndAbove(i)) + ...
                                                vPrAbove(i)*mTransition(1:n_State,i) ;
            mTransition(1:n_State,i) = 0 ;
            
        end
                                                
        % Adjust mTransition for type E agents
        mTransition(n_State+1:end,[OC_E;zeros(n_State,1)]==1) = mTransition(n_State+1:end,[zeros(n_State,1);OC_E]==1) ;
        mTransition(n_State+1:end,[zeros(n_State,1);OC_E]==1) = 0 ;

    % Compute invariant histogram
    errHistogram = 100;	iterationHistogram = 0;
    G = ones(n_State*2,1) ./ (n_State*2);
    while errHistogram > 1e-12 && iterationHistogram < 1e4
        G_New = mTransition' * G;
        errHistogram = max(abs(G_New - G)); 
        iterationHistogram = iterationHistogram + 1;
        G = G_New;
    end

    % Expand histogram matrix
    G_W = G(1:n_State)     ; 
    G_E = G(n_State+1:end) ;
    mDBN_W = reshape(full(G_W),n_A,n_E,n_Z);    
    mDBN_E = reshape(full(G_E),n_A,n_E,n_Z);   
        
        
    %% Update prices 
        A = sum( vA_Grid.*sum(sum(mDBN_W+mDBN_E,3),2) ) ;
        K = min( llambda*mA_Grid , (mmu*p*mZ_Grid.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
        X = sum(sum(sum( (mZ_Grid.*K).^mmu .* mDBN_E ))).^(1/mmu) ;
        N = sum( vE_Grid.*squeeze(sum(sum(mDBN_W,3),1))' ) ;
        w_new = (1-aalpha)*AA*(X/N)^aalpha/(1-tau_n) ;
        p_new = (aalpha)*AA*(X)^(aalpha-mmu)*N^(1-aalpha) ;
        
        f = @(R) ( sum(sum(sum( min( llambda*mA_Grid , (mmu*p_new*mZ_Grid.^mmu/(R+ddelta)).^(1/(1-mmu)) ).*mDBN_E )))/A-1 ) ;
        options = optimoptions('fsolve','Display','off');
        [r_new,err,exitflag] = fsolve(f, r ,options);
                        
    
    % Residual
        residual = [r_new/r-1 ; p_new/p-1 ; w_new/w-1] ;
        residual = sqrt(sum(residual.^2)) ;
    
    %% Optional output
        if nargout>1
            mDBN_W_out = mDBN_W ;
            mDBN_E_out = mDBN_E ;
            if nargout>3
                mAp_W_out = Ap_W_VFI ;
                mAp_E_out = Ap_E_VFI ;
                OC_W_out  = OC_W_VFI ;
                OC_E_out  = OC_E_VFI ;
                VW_out    = VW_VFI   ;
                VE_out    = VE_VFI   ;
            end 
        end 




end 