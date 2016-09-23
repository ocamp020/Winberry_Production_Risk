

function mAssetsPrimeNew = updateCoefficients_splines(mAssetsPrime,Y,r,p,w)

% Updates linear splines approximating savings decision in steady state
% 
% Inputs
%   (1) mAssetsPrime: nEpsilon x nAssets matrix, storing previous iteration's decision along grid
%
% Outputs
%   (1) mAssetsPrimeNew: nEpsilon x nAssets matrix, storing updated decisions
% 
% Thomas Winberry and Alp Tuncay, July 26th, 2016
% Modified by Sergio Ocampo, August 24, 2016

% Declare global variables
global bbeta ggamma ssigma ddelta llambda mmu tau_k aaBar mEpsilonTransition ...
	nEpsilon nAssets nZ nState ...
	vAssetsGrid mEpsilonGrid mAssetsGrid mZGrid ...
	mEpsilonPrimeGrid mZPrimeGrid vIndicesAbove vIndicesBelow vWeightAbove vWeightBelow

%---------------------------------------------------------------
% Compute savings policy next period (mAssetsPrimePrime)
%---------------------------------------------------------------

% Compute weights and indices for linear interpolation
[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,mAssetsPrime(:));

% Compute interpolation for each realization of epsilon
mAssetsPrimePrime = mAssetsPrime(vIndicesBelow,:,:) .* repmat(vWeightBelow,[1 nEpsilon nZ]) + ...
	mAssetsPrime(vIndicesAbove,:,:) .* repmat(vWeightAbove,[1 nEpsilon nZ]);

%---------------------------------------------------------------
% Compute updated savings policy 
%---------------------------------------------------------------

% Compute consumption in next period using mAssetsPrimePrime and budget constraint
mAssetsPrimeGrid = repmat(mAssetsPrime(:),[1 nEpsilon nZ]) ;
K       = min( llambda*mAssetsPrimeGrid , (mmu*p*mZPrimeGrid.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
Pr      = p*(mZPrimeGrid.*K).^mmu - (r+ddelta)*K ;
YPrime  = Pr*(1-tau_k) + (1+r*(1-tau_k))*mAssetsPrimeGrid ;
mConsumptionPrime = w*mEpsilonPrimeGrid + YPrime - mAssetsPrimePrime;
    ECp_aux = NaN(nState,nEpsilon,nZ);
    ECp     = NaN(nAssets,nEpsilon,nZ);
    i = 1; 
    for i_z=1:nZ 
        ECp_aux(:,:,i_z) = mConsumptionPrime(:,:,i_z).^(-ssigma) * mEpsilonTransition' ;
        for i_e=1:nEpsilon
            i_a_min = 1 + nAssets*(i-1) ;
            i_a_max = nAssets*i ;
            ECp(:,i_e,i_z) = ECp_aux(i_a_min:i_a_max,i_e,i_z) ;
            i = i+1 ;
        end
    end 
Mrg_Pr  = (mmu*p*(mZGrid*llambda).^mmu.*mAssetsPrime.^(mmu-1)-(r+ddelta)*llambda) .* (K==llambda*mAssetsPrime) ;
	
% Compute savings policy from Euler Equation(i.e. assuming not borrowing constrained)
mConsumption = ((1+(r+Mrg_Pr)*(1-tau_k))*bbeta*ggamma.*ECp).^(-1/ssigma) ;
mAssetsPrimeNewStar = w * mEpsilonGrid + Y - mConsumption;

% Enforce borrowing constraint
mAssetsPrimeNew = max(mAssetsPrimeNewStar,aaBar * ones(nEpsilon,nAssets));




end 