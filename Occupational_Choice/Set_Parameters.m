% Production Risk Economy with Occupational Choice
% Set Parameters 
% Juan David Herreno & Sergio Ocampo (2016)

%% Preferences and Technology

global bbeta ggamma ssigma aaBar aalpha ddelta mmu llambda AA tau_n tau_k
	
% Preferences
bbeta   = 0.96      ; % discount factor (annual calibration)
ggamma  = (1-1/40)  ; % Survival probability 
ssigma  = 1         ; % coefficient of relative risk aversion
aaBar   = 0.0001    ; % borrowing constraint

% Technology
aalpha  = 0.36      ; % capital share
ddelta  = 0.05      ; % depreciation rate (annual calibration)
mmu     = 0.80      ; % Substitution parameter for intermediate goods
llambda = 1.50      ; % Financial constraint 
AA      = 0.5       ; % TFP

% Taxes 
tau_n = 0.2 ; % Payroll taxes
tau_k = 0.2 ; % Capital income taxes


%% State Space
global n_E n_Z n_A n_State ...
       vA_Grid mA_Grid A_Min A_Max Grid_Curvature ...
       vE_Grid mE_Grid mE_Transition vE_Invariant aggEmployment ...
	   vZ_Grid mZ_Grid mZ_Transition vZ_Invariant vKappa mKappa
	
% Order of approximation
n_E     = 3   ; % number of gridpoints for labor productivity
n_Z     = 5   ; % number of gridpoints for entrepreneurial productivity
n_A     = 75  ; % number of gridpoints for assets
n_State = n_E * n_Z * n_A;

% Bounds on grid space
A_Min = aaBar;	
A_Max = 100  ;
Grid_Curvature = 3.0 ;
vA_Grid = A_Min + linspace(0,1,n_A)'.^Grid_Curvature.*(A_Max-A_Min) ;

% Labor Productivity Shocks
if n_E == 2 
    vE_Grid = [3;7];
    aggEmployment = .93; uDuration = 1;
    mE_Transition = [uDuration / (1 + uDuration), 1 - (uDuration / (1 + uDuration));...
                         ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration))),...
                         1 - ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))];
    vE_Invariant = [1 - aggEmployment;aggEmployment];
elseif n_E == 3
    vE_Grid = [0.1 ; 1 ; 10] ;
    mE_Transition = [0.3  0.7  0.0;
                          0.1  0.8  0.1;
                          0.05 0.65 0.3];
else
    [vE_Grid,mE_Transition] = MarkovAR(n_E,3,0.9,0.03) ;
    vE_Grid = exp(vE_Grid) ;
end 
    [vE_Invariant,~]  = eig(mE_Transition') ;
    vE_Invariant      = vE_Invariant(:,1)/sum(vE_Invariant(:,1)) ;

    aggEmployment = vE_Grid'*vE_Invariant ;

% Entrepreneurial Productivity Types
if n_Z == 2
    vZ_Grid = [0.0 2.5]; 
    mZ_Transition = [0.98 0.02;
                    0.95 0.05] ;
else
    [vZ_Grid,mZ_Transition] = MarkovAR(n_Z,4.5,0.5,0.11) ;
    vZ_Grid = exp(vZ_Grid)' ;

    %vZGrid = [0 1.2 2.3] ;
    %mZTransition = [0.97 0.02  0.01;
    %                0.90 0.08  0.02;
    %                0.93 0.06  0.01];
end                
    [vZ_Invariant,~]  = eig(mZ_Transition') ;
    vZ_Invariant = vZ_Invariant(:,1)/sum(vZ_Invariant(:,1)) ;
    disp('vZGrid vZInvariant'); disp([vZ_Grid' vZ_Invariant])

% Entrepreneurial cost
    vKappa = 0.00*ones(n_Z,1) ;
    mKappa = repmat(reshape(vKappa,[1,1,n_Z]),[n_A,n_E,1]) ;
    
    
% Make matrix versions of grids
mE_Grid = repmat(vE_Grid',[n_A 1 n_Z]);
mA_Grid = repmat(vA_Grid,[1 n_E n_Z]);
mZ_Grid = repmat(reshape(vZ_Grid,[1,1,n_Z]) , [n_A,n_E,1]) ;    

%% Set approximation parameters

global maxIterations tolerance dampening displayOpt

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;


% Whether to print out results from steady state computation
displayOpt = 'iter-detailed';       % 'iter-detailed' or 'off'


