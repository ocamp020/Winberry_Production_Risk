% Sets parameter values 
%
% Thomas Winberry, July 26th, 2016

%----------------------------------------------------------------
% Set economic parameters 
%----------------------------------------------------------------

global bbeta ggamma ssigma aaBar aalpha ddelta mmu llambda AA vEpsilonGrid mEpsilonTransition vEpsilonInvariant vZInvariant aggEmployment ...
	   vZGrid mZTransition rrhoTFP ssigmaTFP tau_n tau_k
	
% Preferences
bbeta  = 0.96      ;								% discount factor (annual calibration)
ggamma = (1-1/40)  ;                                % Survival probability 
ssigma = 1         ;								% coefficient of relative risk aversion
aaBar  = 0.0001    ;								% borrowing constraint

% Technology
aalpha  = 0.36;										% capital share
ddelta  = 0.05;										% depreciation rate (annual calibration)
mmu     = 0.80;										% Substitution parameter for intermediate goods
llambda = 1.50;                                     % Financial constraint 
AA      = 0.5 ;                                     % TFP

% Idioynscratic Shocks
% vEpsilonGrid = [3;7];
% aggEmployment = .93; uDuration = 1;
% mEpsilonTransition = [uDuration / (1 + uDuration), 1 - (uDuration / (1 + uDuration));...
% 					 ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration))),...
% 					 1 - ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))];
% vEpsilonInvariant = [1 - aggEmployment;aggEmployment];

% [vEpsilonGrid,mEpsilonTransition] = MarkovAR_95(3,0.9,0.01) ;
% vEpsilonGrid = exp(vEpsilonGrid) ;

vEpsilonGrid = [0.1 ; 1 ; 10] ;
mEpsilonTransition = [0.3  0.7  0.0;
                      0.1  0.8  0.1;
                      0.05 0.65 0.3];
[vEpsilonInvariant,~]  = eig(mEpsilonTransition') ;
vEpsilonInvariant = vEpsilonInvariant(:,1)/sum(vEpsilonInvariant(:,1)) ;

aggEmployment = vEpsilonGrid'*vEpsilonInvariant ;

% Idiosyncratic Productivity Types
% vZGrid = [0.0 2.5]; 
% mZTransition = [0.98 0.02;
%                 0.95 0.05] ;

% [vZGrid,mZTransition] = MarkovAR_95(3,0.9,0.01) ;
% vZGrid = exp(vZGrid) ; 

vZGrid = [0 0.01 2.5] ;
mZTransition = [0.98 0.00  0.02;
                0.00 0.98  0.02;
                0.50 0.45  0.05];
                  
[vZInvariant,~]  = eig(mZTransition') ;
vZInvariant = vZInvariant(:,1)/sum(vZInvariant(:,1)) ;


% Aggregate Shocks
rrhoTFP = .859;										
ssigmaTFP = .014;

% Taxes 
tau_n = 0.2 ; % Payroll taxes
tau_k = 0.2 ; % Capital income taxes

%----------------------------------------------------------------
% Set approximation parameters
%----------------------------------------------------------------

global nEpsilon nZ nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature AssetsGridCurvature ...
	nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening splineOpt displayOpt
	
% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% if splineOpt = 1, use splines to approximate savings policy; if splineOpt = 0, use polynomials
				% to approximate conditional expectation function

% Whether to print out results from steady state computation
displayOpt = 'iter-detailed';       % 'iter-detailed' or 'off'

% Order of approximation
nEpsilon = numel(vEpsilonGrid);
nZ       = numel(vZGrid);
nAssets  = 25; % number of gridpoints in spline approximation or polynomials in polynomial approximation
nState   = nEpsilon * nZ * nAssets;

% Bounds on grid space
kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / (bbeta*ggamma)) - (1 - ddelta))) ^ (1 / (1 - aalpha));
assetsMin = aaBar;	assetsMax = 100;

% Finer grid for analyzing policy functions
nAssetsFine = 100;
AssetsGridCurvature = 3.0 ;
nStateFine = nEpsilon * nZ * nAssetsFine;

% Approximation of distribution
nMeasure = 3; % Original 3
nAssetsQuadrature = 8;
nStateQuadrature = nEpsilon * nZ * nAssetsQuadrature;
nMeasureCoefficients = nEpsilon * nZ * nMeasure;

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;