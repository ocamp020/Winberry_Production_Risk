% Creates grids to use in various approximations
% 
% Thomas Winberry, July 26th, 2016

%---------------------------------------------------------------
% Grids for approximating individual decisions
%---------------------------------------------------------------

if splineOpt == 0	% if using polynomials and approximating conditional expectation

	global vAssetsGridZeros vAssetsGrid mEpsilonGrid mZGrid mAssetsGrid mEpsilonPrimeGrid mZPrimeGrid vZGrid

	% Zeros of chebyshev polynomial
	vAssetsGridZeros = -cos(((2 * (1:nAssets)-1)' * pi) / (2 * nAssets));

	% Scale up to state space
	vAssetsGrid = scaleUp(vAssetsGridZeros,assetsMin,assetsMax);
	
else	% if using splines and approximating asset accumulation

	global vAssetsGrid mEpsilonGrid mZGrid mAssetsGrid mEpsilonPrimeGrid mZPrimeGrid vZGrid

	% Grid over assets
	vAssetsGrid = exp(linspace(log(assetsMin + .01),log(assetsMax + .01),nAssets)');
	vAssetsGrid = vAssetsGrid - .01;
	
end
	
% Make matrix versions of the grids
mEpsilonGrid = repmat(vEpsilonGrid',[nAssets 1 nZ]);
mAssetsGrid  = repmat(vAssetsGrid,[1 nEpsilon nZ]);
mZGrid       = repmat(reshape(vZGrid,[1,1,nZ]) , [nAssets,nEpsilon,1]) ;
mEpsilonPrimeGrid = repmat(vEpsilonGrid',[nState 1 nZ]) ;
mZPrimeGrid       = repmat(reshape(vZGrid,[1,1,nZ]) , [nState,nEpsilon,1]) ;

%---------------------------------------------------------------
% Fine grid, for histogram and plotting functions
%---------------------------------------------------------------

global vAssetsGridFine vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mZGridFine AssetsGridCurvature

% Assets grid
vAssetsGridFine = assetsMin + linspace(0,1,nAssetsFine)'.^AssetsGridCurvature.*(assetsMax-assetsMin) ;
% vAssetsGridFine = linspace(assetsMin,assetsMax,nAssetsFine)';

% Scale down to [-1,1]
vAssetsGridFineZeros = scaleDown(vAssetsGridFine,assetsMin,assetsMax);

% Make matrix versions of grids
mEpsilonGridFine = repmat(vEpsilonGrid',[nAssetsFine 1 nZ]);
mAssetsGridFine  = repmat(vAssetsGridFine,[1 nEpsilon nZ]);
mZGridFine       = repmat(reshape(vZGrid,[1,1,nZ]) , [nAssetsFine,nEpsilon,1]) ;

%---------------------------------------------------------------
% Quadrature grid, to integrate density (away from borrowing constraint)
%---------------------------------------------------------------

global vQuadratureWeights vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature ...
	mAssetsGridQuadrature mZGridQuadrature

% Compute grid in the interval [-1, 1]
[vAssetsGridQuadratureZeros,vQuadratureWeights] = computeGaussLegendreQuadrature(nAssetsQuadrature);

% Scale up grid
vAssetsGridQuadrature = scaleUp(vAssetsGridQuadratureZeros,assetsMin+1e-1,assetsMax);

% Make matrix versions of the grids
mEpsilonGridQuadrature = repmat(vEpsilonGrid',[nAssetsQuadrature 1 nZ]);
mAssetsGridQuadrature  = repmat(vAssetsGridQuadrature,[1 nEpsilon nZ]);
mZGridQuadrature       = repmat(reshape(vZGrid,[1,1,nZ]) , [nAssetsQuadrature,nEpsilon,1]) ;

