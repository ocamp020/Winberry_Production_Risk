% Production Risk Economy with Occupational Choice
% Set Parameters 
% Juan David Herreno & Sergio Ocampo (2016)

%% Preferences and Technology

global bbeta ggamma ssigma aaBar aalpha ddelta mmu llambda AA tau_n tau_k
	
% Preferences
bbeta   = 0.98      ; % discount factor (annual calibration)
ggamma  = (1-1/40)  ; % Survival probability 
ssigma  = 1         ; % coefficient of relative risk aversion
aaBar   = 0.0001    ; % borrowing constraint

% Technology
aalpha  = 0.36      ; % capital share
ddelta  = 0.025     ; % depreciation rate (annual calibration)
mmu     = 0.80      ; % Substitution parameter for intermediate goods
llambda = 1.50      ; % Financial constraint 
AA      = 0.5       ; % TFP

% Labor frictions
b       = 0.1       ; % Unemployement benefits: Replacement rate
jfr_W   = 0.30      ; % Job finding rate from unemployment
jdr_W   = 0.04      ; % Job destruction rate
jfr_E   = 0.30      ; % Job finding rate from unemployment
jdr_E   = 0.04      ; % Job destruction rate

% Taxes 
tau_n = 0.2 ; % Payroll taxes
tau_k = 0.2 ; % Capital income taxes


%% State Space
global n_E n_Z n_A n_State ...
       vA_Grid mA_Grid A_Min A_Max Grid_Curvature ...
       vE_Grid mE_Grid mE_Transition_W mE_Transition_E ...
       mE_Transition_W_VFI mE_Transition_E_VFI ...
	   vZ_Grid mZ_Grid mZ_Transition vZ_Invariant vKappa mKappa
	
% Order of approximation
n_E     = 3   ; % number of gridpoints for labor productivity
n_Z     = 5   ; % number of gridpoints for entrepreneurial productivity
n_A     = 150  ; % number of gridpoints for assets
n_State = n_E * n_Z * n_A;

% Bounds on grid space
A_Min = aaBar;	
A_Max = 50  ;
Grid_Curvature = 3.0 ;
vA_Grid = A_Min + linspace(0,1,n_A)'.^Grid_Curvature.*(A_Max-A_Min) ;
% figure; plot(vA_Grid,vA_Grid,'-o'); title('Asset Grid - Curvature');

% Labor Productivity Shocks
if n_E == 2 
    vE_Grid = [3;7];
    aggEmployment = .93; uDuration = 1;
    mE_Transition = [uDuration / (1 + uDuration), 1 - (uDuration / (1 + uDuration));...
                         ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration))),...
                         1 - ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))];
    vE_Invariant = [1 - aggEmployment;aggEmployment];
elseif n_E == 3
    vE_Grid = [0.05 ; 0.9 ; 10] ;
    mE_Transition = [0.1  0.8  0.1;
                     0.1  0.8  0.1;
                     0.05 0.65 0.3];
    mE_Transition_W = [0.25  0.70 0.05;
                       0.08  0.79 0.13;
                       0.08  0.63 0.29];
    mE_Transition_E = [0.40  0.55 0.05;
                       0.08  0.81 0.11;
                       0.07  0.67 0.26];
else
    [vE_Grid,mE_Transition] = MarkovAR(n_E-1,3,0.885,0.1) ;
    vE_Grid = exp(vE_Grid) ;
    vE_Grid = [b ; vE_Grid] ;
    vE_Grid = [b ; 1 ; 10 ; 20 ] ;
    mE_Transition_W = [1-jfr_W , jfr_W , zeros(1,n_E-2);
                     jdr_W*ones(n_E-1,1) (1-jdr_W)*mE_Transition];
    mE_Transition_E = [1-jfr_E , jfr_E , zeros(1,n_E-2);
                     jdr_E*ones(n_E-1,1) (1-jdr_E)*mE_Transition];
end 
    [vE_Invariant_W,~]  = eig(mE_Transition_W') ;
    vE_Invariant_W      = vE_Invariant_W(:,1)/sum(vE_Invariant_W(:,1)) ;
    [vE_Invariant_E,~]  = eig(mE_Transition_E') ;
    vE_Invariant_E      = vE_Invariant_E(:,1)/sum(vE_Invariant_E(:,1)) ;
    
    % Matrix for VFI
    mE_Transition_W_VFI = NaN(n_A,n_E,n_Z,n_A,n_E) ; 
    mE_Transition_E_VFI = NaN(n_A,n_E,n_Z,n_A,n_E) ; 
    for i_e = 1:n_E 
    for i_ep = 1:n_E 
        mE_Transition_W_VFI(:,i_e,:,:,i_ep) = mE_Transition_W(i_e,i_ep) ;
        mE_Transition_E_VFI(:,i_e,:,:,i_ep) = mE_Transition_E(i_e,i_ep) ;
    end
    end
    
    disp('vEGrid vEInvariant_W vE_Invariant_E'); disp([vE_Grid vE_Invariant_W vE_Invariant_E])

% Entrepreneurial Productivity Types
if n_Z == 2
    vZ_Grid = [0.0 2.5]; 
    mZ_Transition = [0.98 0.02;
                    0.95 0.05] ;
else
    [vZ_Grid,mZ_Transition] = MarkovAR(n_Z,4.5,0.5,0.10) ;
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


