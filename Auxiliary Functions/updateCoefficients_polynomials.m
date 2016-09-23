function mCoefficientsNew = updateCoefficients_polynomials(mCoefficients,Y,r,p,w)

% Updates polynomial coefficients approximating the conditional expectation function in steady state
% 
% Inputs
%   (1) mCoefficients: nAssets x nEpsilon x nZ matrix, storing previous iteration's coefficients
%
% Outputs
%   (1) mCoefficientsNew: nAssets x nEpsilon x nZ matrix, storing updated coefficients
% 
% Thomas Winberry, January 19, 2016
% Modified by Sergio Ocampo, August 24, 2016

% Declare global variables
global  bbeta ggamma ssigma ddelta llambda aaBar mmu tau_k mEpsilonTransition ...
        nEpsilon nAssets nZ nState assetsMin assetsMax ...
        mEpsilonGrid mAssetsPrimeGrid mEpsilonPrimeGrid mZPrimeGrid ...
        vAssetsPoly vAssetsPolySquared 

%---------------------------------------------------------------
% Compute current period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
    mConditionalExpectation = NaN(nAssets,nEpsilon,nZ) ;
for i_z=1:nZ
    mConditionalExpectation(:,:,i_z) = exp(vAssetsPoly * mCoefficients(:,:,i_z));
end 

% Compute target saving
mAssetsPrimeStar = w*mEpsilonGrid + Y - (mConditionalExpectation .^ (-1 / ssigma));

% Compute actual saving
mAssetsPrime = max(mAssetsPrimeStar,aaBar );
mAssetsPrimeGrid = repmat(mAssetsPrime(:),[1 nEpsilon nZ]);

% Compute next period's polynomials
mAssetsPrimeZeros = scaleDown(mAssetsPrime,assetsMin,assetsMax);
mPolyAssetsPrime  = computeChebyshev(nAssets, mAssetsPrimeZeros(:));

%---------------------------------------------------------------
% Compute next period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectationPrime = NaN(nState,nEpsilon,nZ) ;
for i_z=1:nZ
mConditionalExpectationPrime(:,:,i_z) = exp(mPolyAssetsPrime * mCoefficients(:,:,i_z) );
end 

% Compute target saving
K       = min( llambda*mAssetsPrimeGrid , (mmu*p*mZPrimeGrid.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
Pr      = p*(mZPrimeGrid.*K).^mmu - (r+ddelta)*K ;
YPrime  = Pr*(1-tau_k) + (1+r*(1-tau_k))*mAssetsPrimeGrid ;
Mrg_Pr  = (mmu*p*(mZPrimeGrid*llambda).^mmu.*mAssetsPrimeGrid.^(mmu-1)-(r+ddelta)*llambda) .* (K==llambda*mAssetsPrimeGrid) ;
mAssetsPrimePrimeStar = w * mEpsilonPrimeGrid + YPrime - (mConditionalExpectationPrime .^ (-1 / ssigma));

% Compute actual savings
mAssetsPrimePrimeGrid = max(mAssetsPrimePrimeStar,aaBar);

%---------------------------------------------------------------
% Update conditional expectation function
%---------------------------------------------------------------

% Compute new conditional expectation function
mConsumptionPrime = w * mEpsilonPrimeGrid + YPrime - mAssetsPrimePrimeGrid;
    ECp    = NaN(nAssets*nEpsilon*nZ,nEpsilon,nZ) ;
        for i_z=1:nZ
           ECp(:,:,i_z) = mConsumptionPrime(:,:,i_z).^(-ssigma)*mEpsilonTransition' ; 
        end
EVp    = (1+(r+Mrg_Pr)*(1-tau_k))*bbeta*ggamma.*ECp ;

mConditionalExpectation = NaN(nAssets,nEpsilon,nZ);
i = 1 ;
for i_z = 1:nZ
for i_e = 1:nEpsilon
    i_a_min = 1 + nAssets*(i-1) ;
    i_a_max = nAssets*i ;
    mConditionalExpectation(:,i_e,i_z) = EVp(i_a_min:i_a_max,i_e,i_z) ;
    i = i+1;
end 
end 

% aConditionalExpectationTilde = reshape(EVp,nEpsilon,nEpsilon,nAssets);
% 
% % Extract the relevant entries
% mConditionalExpectation = zeros(nEpsilon,nAssets);
% for iEpsilon = 1:nEpsilon
% 	mConditionalExpectation(iEpsilon,:) = aConditionalExpectationTilde(iEpsilon,iEpsilon,:);
% end

% Update the coefficients
mCoefficientsNew = zeros(nAssets,nEpsilon,nZ);
for iZ       = 1:nZ
for iEpsilon = 1:nEpsilon	% interpolate
    vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * log(mConditionalExpectation(:,iEpsilon,iZ)')),2);
    mCoefficientsNew(:,iEpsilon,iZ) = (vCoefficients ./ vAssetsPolySquared)';
end
end



end 

