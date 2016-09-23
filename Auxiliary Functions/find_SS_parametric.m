function [residual,mCoefficients_out,mParameters_out,mMoments_out,mHat_out] = find_SS_parametric(x,mMoments,mHat,aGridMoments)

    % Initialization (global variables)
    global ssigma bbeta ggamma aalpha ddelta mmu llambda AA aggEmployment tau_n tau_k aaBar ...
    	   mZTransition mEpsilonTransition vZInvariant vEpsilonInvariant ...
           nEpsilon nZ nAssets vAssetsGrid mAssetsGrid mEpsilonGrid mZGrid vZGrid ...
           vAssetsPoly vAssetsPolySquared vAssetsPolyBC ...
           nAssetsQuadrature vAssetsGridQuadrature vAssetsPolyQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature mZGridQuadrature ...
           nMeasure vQuadratureWeights ...
           dampening tolerance maxIterations splineOpt
       
    % Set error tolerance & max iteration depending on use
        if nargout == 1 
            err1 = tolerance;
            err2 = 1e-4;
            tol2 = 200;
        elseif nargout > 1 
            err1 = 1e-8;
            err2 = 1e-6;
            tol2 = 500;
        else
            error('Check # of inputs & outputs for computeMCResidualPolynomials');
        end

    % Declare prices
        r = x(1); 
        p = x(2);
        X = (p*aggEmployment^(aalpha-1)/(AA*aalpha))^(1/(aalpha-mmu));
        w = (1-aalpha)*AA*(X/aggEmployment)^aalpha/(1-tau_n) ;
        
    % Firm's capital demand and profits
        K       = min( llambda*mAssetsGrid , (mmu*p*mZGrid.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
        Pr      = p*(mZGrid.*K).^mmu - (r+ddelta)*K ;
        Mrg_Pr  = (mmu*p*(mZGrid*llambda).^mmu.*mAssetsGrid.^(mmu-1)-(r+ddelta)*llambda) .* (K==llambda*mAssetsGrid) ;
        Y       = Pr*(1-tau_k) + (1+r*(1-tau_k))*mAssetsGrid ;
        
        K_Q       = min( llambda*mAssetsGridQuadrature , (mmu*p*mZGridQuadrature.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
        Pr_Q      = p*(mZGridQuadrature.*K_Q).^mmu - (r+ddelta)*K_Q ;
        Y_Q       = Pr_Q*(1-tau_k) + (1+r*(1-tau_k))*mAssetsGridQuadrature ;
        
    
    %% Agent's policy functions (by polynomial or spline approximation)
    if splineOpt == 0	% approximate conditional expectation function using polynomials
        
        % Initialize coefficients using rule of thumb savings rule
            C      = w*mEpsilonGrid + (Pr + r*mAssetsGrid)*(1-tau_k) ;
            ECp    = NaN(nAssets,nEpsilon,nZ) ;
                for i_z=1:nZ
                   ECp(:,:,i_z) = C(:,:,i_z).^(-ssigma)*mEpsilonTransition' ; 
                end
            EVp    = (1+(r+Mrg_Pr)*(1-tau_k))*bbeta*ggamma.*ECp ;
        
        mGridInit = log(EVp);
        mCoefficients = zeros(nAssets,nEpsilon,nZ);
        for iZ       = 1:nZ
        for iEpsilon = 1:nEpsilon	% interpolate
            vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * mGridInit(:,iEpsilon,iZ)'),2);
            mCoefficients(:,iEpsilon,iZ) = (vCoefficients ./ vAssetsPolySquared)';
        end
        end 

        % Iterate
        err = 100; iteration = 1;
        while err > tolerance && iteration <= maxIterations

            mCoefficientsNew = updateCoefficients_polynomials(mCoefficients,Y,r,p,w);
            err = max(abs(mCoefficientsNew(:) - mCoefficients(:)));
            iteration = iteration + 1;
            mCoefficients = dampening * mCoefficients + (1 - dampening) * mCoefficientsNew;

        end

        mCoefficients_out = mCoefficients;
        
    else	% approximate savings decision using linear splines
        
        % Initialize coefficients
        mAssetsPrime = mAssetsGrid;

        % Iterate
        err = 100; iteration = 1;
        while err > tolerance && iteration <= maxIterations
            mAssetsPrimeNew = updateCoefficients_splines(mAssetsPrime,Y,r,p,w);
            err = max(abs(mAssetsPrimeNew(:) - mAssetsPrime(:)));
            iteration = iteration + 1;
            mAssetsPrime = dampening * mAssetsPrime + (1 - dampening) * mAssetsPrimeNew;
        end

        mCoefficients_out = mAssetsPrime;
        
    end
    
    %% Compute policies over quadrature grid for integration
    if splineOpt == 0
       
        % Compute conditional expectation
        mConditionalExpectation = NaN(nAssetsQuadrature,nEpsilon,nZ) ;
        for i_z=1:nZ
            mConditionalExpectation(:,:,i_z) = exp(vAssetsPolyQuadrature * mCoefficients(:,:,i_z));
        end 

        % Compute savings policy
        mAssetsPrimeStar = w * mEpsilonGridQuadrature + Y_Q  - (mConditionalExpectation .^ (-1 / ssigma));
        mAssetsPrimeQuadrature = max(mAssetsPrimeStar,aaBar );

    else

        % Compute weights
        [vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,vAssetsGridQuadrature);

        % Linear interpolation
        mAssetsPrimeQuadrature = mAssetsPrime(vIndicesBelow,:,:) .* repmat(vWeightBelow,[1 nEpsilon nZ]) + ...
            mAssetsPrime(vIndicesAbove,:,:) .* repmat(vWeightAbove,[1 nEpsilon nZ]);

    end
    
    %% Compute policies at borrowing constraint for integration
    if splineOpt == 0

        % Compute conditional expectation
        mConditionalExpectation = NaN(1,nEpsilon,nZ) ;
        for i_z=1:nZ
            mConditionalExpectation(:,:,i_z) = exp(vAssetsPolyBC * mCoefficients(:,:,i_z));
        end 

        % Compute savings policy
        mAssetsPrimeStar = w * mEpsilonGrid(1,:,:) + Y(1,:,:) - (mConditionalExpectation .^ (-1 / ssigma));
        mAssetsPrimeBC = max(mAssetsPrimeStar,aaBar );

    else

        % Compute weights
        [vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,assetsMin);

        % Linear interpolation
        mAssetsPrimeBC = mAssetsPrime(vIndicesBelow,:,:) .* repmat(vWeightBelow,[1 nEpsilon nZ]) + ...
            mAssetsPrime(vIndicesAbove,:,:) .* repmat(vWeightAbove,[1 nEpsilon nZ]);

    end

    
    %% Stationary distribution (Parametric Family)
    
    % Initialize iteration
        err = 100; iteration = 1; 
        options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
            'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000); %,'iter'

    % Iteration
    while err > err2 && iteration <= tol2

        %%%
        % Update density away from borrowing constraint
        %%%

        % Compute parameters of the distribution by minimization
        mParameters = zeros(nEpsilon,nZ,nMeasure+1);
        for iZ       = 1 : nZ
        for iEpsilon = 1 : nEpsilon
            objectiveFunction = @(vParametersTilde) parametersResidual(vParametersTilde,squeeze(aGridMoments(:,iEpsilon,iZ,:)));
            [vParameters,normalization] = fminunc(objectiveFunction,zeros(nMeasure,1),options);
            % [vParameters,normalization] = fminsearch(objectiveFunction,zeros(nMeasure,1));
            mParameters(iEpsilon,iZ,:) = [1 / normalization; vParameters];
        end
        end

        % Compute new moments and centered moments grid
        mMomentsNew = zeros(nEpsilon,nZ,nMeasure);
        aGridMomentsNew = zeros(nAssetsQuadrature,nEpsilon,nZ,nMeasure);

        for iZTilde       = 1 : nZ
        for iEpsilonTilde = 1 : nEpsilon

            % Compute first moment (uncentered)
            mMomentsNew(iEpsilonTilde,iZTilde,1) = 0;      
            for iEpsilon = 1 : nEpsilon
                % First count agents that do not die so that their Z is equal to iZTilde
                % Only their Epsilon changes according to mEpsilonTransition
                mMomentsNew(iEpsilonTilde,iZTilde,1) = mMomentsNew(iEpsilonTilde,iZTilde,1) + ggamma* vZInvariant(iZTilde) * vEpsilonInvariant(iEpsilon) * mEpsilonTransition(iEpsilon,iEpsilonTilde) *( ...
                    (1 - mHat(iEpsilon,iZTilde)) * ...
                    (mParameters(iEpsilon,iZTilde,1)*vQuadratureWeights'*( mAssetsPrimeQuadrature(:,iEpsilon,iZTilde).* exp(squeeze(aGridMoments(:,iEpsilon,iZTilde,:)) * squeeze(mParameters(iEpsilon,iZTilde,2:nMeasure+1)))) ) +...
                    mHat(iEpsilon,iZTilde) * mAssetsPrimeBC(1,iEpsilon,iZTilde) ...
                    ) ;
            for iZ       = 1 : nZ
                % Now count agents that die thier original Z is iZ
                % Epsilon changes according to mEpsilonTransition
                mMomentsNew(iEpsilonTilde,iZTilde,1) = mMomentsNew(iEpsilonTilde,iZTilde,1) + (1-ggamma)* vZInvariant(iZ)*mZTransition(iZ,iZTilde)* vEpsilonInvariant(iEpsilon) * mEpsilonTransition(iEpsilon,iEpsilonTilde)*( ...
                    (1 - mHat(iEpsilon,iZ)) * (mParameters(iEpsilon,iZ,1)*vQuadratureWeights'*( mAssetsPrimeQuadrature(:,iEpsilon,iZ).* exp(squeeze(aGridMoments(:,iEpsilon,iZ,:)) * squeeze(mParameters(iEpsilon,iZ,2:nMeasure+1)))) ) +...
                    mHat(iEpsilon,iZ)  * mAssetsPrimeBC(1,iEpsilon,iZ) ...
                    ) ;	

            end
            end 

            mMomentsNew(iEpsilonTilde,iZTilde,1) = mMomentsNew(iEpsilonTilde,iZTilde,1) / (vEpsilonInvariant(iEpsilonTilde)*vZInvariant(iZTilde));
            aGridMomentsNew(:,iEpsilonTilde,iZTilde,1) = vAssetsGridQuadrature - mMomentsNew(iEpsilonTilde,iZTilde,1);


            % Compute higher order moments (centered)
            for iMoment = 2 : nMeasure

                mMomentsNew(iEpsilonTilde,iZTilde,iMoment) = 0; 
                for iEpsilon = 1 : nEpsilon
                    % First count agents that do not die so that their Z is equal to iZTilde
                    % Only their Epsilon changes according to mEpsilonTransition
                    mMomentsNew(iEpsilonTilde,iZTilde,iMoment) = mMomentsNew(iEpsilonTilde,iZTilde,iMoment) + ggamma* vZInvariant(iZTilde)* vEpsilonInvariant(iEpsilon) * mEpsilonTransition(iEpsilon,iEpsilonTilde)*( ...
                        (1 - mHat(iEpsilon,iZTilde,1)) * (mParameters(iEpsilon,iZTilde,1)*vQuadratureWeights'*( ((mAssetsPrimeQuadrature(:,iEpsilon,iZTilde) - mMomentsNew(iEpsilonTilde,iZTilde,1)) .^ iMoment) ...
                        .* exp(squeeze(aGridMoments(:,iEpsilon,iZTilde,:)) * squeeze(mParameters(iEpsilon,iZTilde,2:nMeasure+1)))) ) +...
                        mHat(iEpsilon,iZTilde,1) * ((mAssetsPrimeBC(1,iEpsilon,iZTilde)-mMomentsNew(iEpsilonTilde,iZTilde,1)).^iMoment) ...
                        ) ;
                for iZ       = 1 : nZ
                    % Now count agents that dieo thier origina Z is iZ
                    % Epsilon changes according to mEpsilonTransition
                    mMomentsNew(iEpsilonTilde,iZTilde,iMoment) = mMomentsNew(iEpsilonTilde,iZTilde,iMoment) + (1-ggamma)* vZInvariant(iZ)*mZTransition(iZ,iZTilde)* vEpsilonInvariant(iEpsilon) * mEpsilonTransition(iEpsilon,iEpsilonTilde)*( ...
                        (1 - mHat(iEpsilon,iZ,1)) * (mParameters(iEpsilon,iZ,1)*vQuadratureWeights'*( ((mAssetsPrimeQuadrature(:,iEpsilon,iZ) - mMomentsNew(iEpsilonTilde,iZTilde,1)) .^ iMoment) ...
                        .* exp(squeeze(aGridMoments(:,iEpsilon,iZ,:)) * squeeze(mParameters(iEpsilon,iZ,2:nMeasure+1)))) ) +...
                        mHat(iEpsilon,iZ,1) * ((mAssetsPrimeBC(1,iEpsilon,iZ)-mMomentsNew(iEpsilonTilde,iZTilde,1)).^iMoment ) ...
                        ) ;	
                end
                end 

                mMomentsNew(iEpsilonTilde,iZTilde,iMoment) = mMomentsNew(iEpsilonTilde,iZTilde,iMoment) / (vEpsilonInvariant(iEpsilonTilde)*vZInvariant(iZTilde));
                aGridMomentsNew(:,iEpsilonTilde,iZTilde,iMoment) = (vAssetsGridQuadrature - mMomentsNew(iEpsilonTilde,iZTilde,1)).^iMoment-mMomentsNew(iEpsilonTilde,iZTilde,iMoment);

            end

        end
        end
        
        %%%
        % Update mass at borrowing constraint
        %%%

        mHatNew = zeros(nEpsilon,nZ);

        for iZTilde       = 1 : nZ
        for iEpsilonTilde = 1 : nEpsilon
            
            for iEpsilon = 1 : nEpsilon

                mHatNew(iEpsilonTilde,iZTilde) = mHatNew(iEpsilonTilde,iZTilde) + ggamma* vZInvariant(iZTilde)* vEpsilonInvariant(iEpsilon) * mEpsilonTransition(iEpsilon,iEpsilonTilde)*(...
                    (1 - mHat(iEpsilon,iZTilde)) * mParameters(iEpsilon,iZTilde,1) * vQuadratureWeights' * ...
                    ((mAssetsPrimeQuadrature(:,iEpsilon,iZTilde) <= aaBar + 1e-8) .* exp(squeeze(aGridMoments(:,iEpsilon,iZTilde,:)) * squeeze(mParameters(iEpsilon,iZTilde,2:nMeasure+1))))  + ...
                    mHat(iEpsilon,iZTilde) * (mAssetsPrimeBC(1,iEpsilon,iZTilde) <= aaBar + 1e-8) ...
                    );
                
                for iZ = 1 : nZ
                    mHatNew(iEpsilonTilde,iZTilde) = mHatNew(iEpsilonTilde,iZTilde) + (1-ggamma)* vZInvariant(iZ)*mZTransition(iZ,iZTilde)* vEpsilonInvariant(iEpsilon) * mEpsilonTransition(iEpsilon,iEpsilonTilde)*(...
                    (1 - mHat(iEpsilon,iZ)) * mParameters(iEpsilon,iZ,1) * vQuadratureWeights' * ...
                    ((mAssetsPrimeQuadrature(:,iEpsilon,iZ) <= aaBar + 1e-8) .* exp(squeeze(aGridMoments(:,iEpsilon,iZ,:)) * squeeze(mParameters(iEpsilon,iZ,2:nMeasure+1))))  + ...
                    mHat(iEpsilon,iZ) * (mAssetsPrimeBC(1,iEpsilon,iZ) <= aaBar + 1e-8) ...
                    );
                    
                end
            end

            mHatNew(iEpsilonTilde,iZTilde) = mHatNew(iEpsilonTilde,iZTilde) / (vEpsilonInvariant(iEpsilonTilde)*vZInvariant(iZTilde));

        end
        end
        

        %%%
        % Update iteration
        %%%

        err = max([max(abs(mMomentsNew(:) - mMoments(:))),max(abs(mHatNew(:) - mHat(:)))]);
        iteration = iteration + 1;
        mMoments = mMomentsNew;
        aGridMoments = aGridMomentsNew;
        mHat = mHatNew;           

    end

    
    
    %% Update prices and residual
        mEpsilonInvariant = repmat(vEpsilonInvariant,1,nZ)  ; 
        mZInvariant       = repmat(vZInvariant',nEpsilon,1) ; 
        A = sum(sum( (mEpsilonInvariant.*mZInvariant.* ( (1 - mHat)*squeeze(mMoments(:,:,1)) + aaBar * mHat) ) ));
        K = min( llambda*mAssetsGridQuadrature , (mmu*p*mZGridQuadrature.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
        X = 0 ;
        for iZ       = 1:nZ
        for iEpsilon = 1:nEpsilon
            X = X + vZInvariant(iZ)*vEpsilonInvariant(iEpsilon)*(...
                (1-mHat(iEpsilon,iZ))*mParameters(iEpsilon,iZ,1)*vQuadratureWeights'*( (vZGrid(iZ)*K(:,iEpsilon,iZ)).^mmu .* exp(squeeze(aGridMoments(:,iEpsilon,iZ,:)) * squeeze(mParameters(iEpsilon,iZ,2:nMeasure+1))))  + ...
                mHat(iEpsilon,iZ)*(vZGrid(iZ)*K(1,iEpsilon,iZ)).^mmu ...
                );
        end
        end
        X = X^(1/mmu) ;
        % w_new = (1-aalpha)*AA*(X/aggEmployment)^aalpha/(1-tau_n) ;
        p_new = (aalpha)*AA*(X)^(aalpha-mmu)*aggEmployment^(1-aalpha) ;
        
        f = @(R) Capital_Market_Clearing(R,A,p_new,mParameters,aGridMoments,mHat) ;
        options = optimoptions('fsolve','Display','off');
        [r_new,err,exitflag] = fsolve(f, r ,options);
                        
    
    % Residual
        residual = [r_new/r-1 ; p_new/p-1] ;
        
        
    %% Optional output
    if nargout > 2

        mParameters_out = mParameters;
        mMoments_out = mMoments;
        mHat_out = mHat;

    end





end 

%% Function for capital market clearing (Finding R)
function res = Capital_Market_Clearing(R,A,p_new,mParameters,aGridMoments,mHat)

    global llambda mmu ddelta nEpsilon nZ mZGridQuadrature mAssetsGridQuadrature vZInvariant vEpsilonInvariant vQuadratureWeights nMeasure
    
    K_mat = min( llambda*mAssetsGridQuadrature , (mmu*p_new*mZGridQuadrature.^mmu/(R+ddelta)).^(1/(1-mmu)) )  ;

    K = 0 ;
    for iZ       = 1:nZ
    for iEpsilon = 1:nEpsilon
        K = K + vZInvariant(iZ)*vEpsilonInvariant(iEpsilon)*(...
                (1-mHat(iEpsilon,iZ))*mParameters(iEpsilon,iZ,1)*vQuadratureWeights'*( K_mat(:,iEpsilon,iZ) .* exp(squeeze(aGridMoments(:,iEpsilon,iZ,:)) * squeeze(mParameters(iEpsilon,iZ,2:nMeasure+1)))) + ...
                mHat(iEpsilon,iZ)*K_mat(1,iEpsilon,iZ) ...
                );
    end
    end
    
    res = K/A-1 ;
end


