% Production Risk Economy with Occupational Choice
% Main File
% Juan David Herreno & Sergio Ocampo (2016)

clear all; close all; clc;

delete('Log_Occupational_Choice.txt')
diary('Log_Occupational_Choice.txt')

%% Set up problem 

    Set_Parameters;

%% Stationary Equilibrium (Histogram Method and discrete VFI)

    t0 = tic;
    x_0 = [0.416914281250000];
%     options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
%     [x,err,exitflag] = fsolve(@(x) Find_DBN_Histogram(x),x_0,options);
%     options = optimset('Display','iter','TolFun',1e-05,'TolX',1e-06);
%     [x,err,exitflag] = fminsearch(@(x) Find_DBN_Histogram(x),x_0,options);
%     fprintf('Done! Time to compute: %2.2f seconds, error=%2.2d \n\n',toc(t0),err)
    
    % Graph of excess supply function
%     w_grid = linspace(0.41,0.42,100) ;
%     for i=1:numel(w_grid)
%         es(i) = Find_DBN_Histogram(w_grid(i)) ;
%         disp([i,w_grid(i),es(i)])
%     end 
%     figure; hold on; plot(w_grid,es,'-o'); plot(w_grid,zeros(size(w_grid))); hold off; xlim([w_grid(1) w_grid(end)]);

    % fzero
    options = optimset('Display','iter','TolFun',5e-04);
    [x,err] = fzero(@(x)Find_DBN_Histogram(x),0.50,options);
  

    % Bisection 
%     w_low  = 0.41 ; 
%     w_high = 0.42 ; 
%     excess_supply_low  = Find_DBN_Histogram(w_low)  ; 
%     excess_supply_high = Find_DBN_Histogram(w_high) ;
%     disp('Starting Bisection for Equilibrium Wage')
%     disp({'iter','w_low','w','w_high','excess_supply','w_dist'})
%     if excess_supply_low*excess_supply_high<0 
%         excess_supply = 1 ;
%         iter = 0 ;
%         while abs(excess_supply)>5e-4 && iter<maxIterations
%             w = (w_low+w_high)/2 ;
%             excess_supply = Find_DBN_Histogram(w)  ; 
%             if excess_supply>0 
%                 w_high = w ; 
%             else 
%                 w_low  = w ;
%             end 
%             iter = iter + 1 ;
%             disp(num2cell([iter w_low w w_high excess_supply w_high-w_low]))
%         end 
%         x = w ;
%     else
%         disp('Equilibrium not bracketed')
%         return
%     end 
    

    r_ben = r    ;
    w_ben = x(1) ; 
    
    Mat = [r_ben w_ben] ; 
    Mat = [{'Aggregate Prices',' ';'Interest Rate','Wage'};num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    Mat = (AA*vZ_Grid.*(aalpha/(r_ben+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w_ben))^mmu).^(1/(1-aalpha-mmu)) ; 
    Mat = [{'Optimal Capital'} cell(1,n_Z-1);num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

%% Distribution, Value and Policy Functions 

    [price_residual,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E,Transition] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual);
    
    [A_ben,K_ben,N_ben,Y_ben,Earnings_W_ben,Earnings_E_ben] = ...
        Graphs_Tables('ben',r_ben,w_ben,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E,Transition) ;
 
return
%% Experiment 1: Change in household's borrowing constraint

    % Set new borrowing limit
    aaBar = -w_ben ;
    disp(' ')
    disp('Changing borrowing constraint to:')
    disp(-aaBar/w_ben);
    disp(' ')
    
    % Adjust Grid
    A_Min = aaBar;	
    A_Max = 50  ;
    Grid_Curvature = 3.0 ;
    vA_Grid = A_Min + linspace(0,1,n_A)'.^Grid_Curvature.*(A_Max-A_Min) ;
    mA_Grid = repmat(vA_Grid,[1 n_E n_Z]);



    % Solve the model
    t0 = tic;
    x_0 = x ;
%     options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
%     [x,err,exitflag] = fsolve(@(x) Find_DBN_Histogram(x),x_0,options);
%     options = optimset('Display','iter','TolFun',1e-05,'TolX',1e-06);
%     [x,err,exitflag] = fminsearch(@(x) Find_DBN_Histogram(x),x_0,options);
%     fprintf('Done! Time to compute: %2.2f seconds, error=%2.2d \n\n',toc(t0),err)

%     % Graph of excess supply function
%     w_grid = linspace(0.40,0.42,100) ;
%     for i=1:numel(w_grid)
%         es(i) = Find_DBN_Histogram(w_grid(i)) ;
%         disp([i,w_grid(i),es(i)])
%     end 
%     figure; hold on; plot(w_grid,es,'-o'); plot(w_grid,zeros(size(w_grid))); hold off; xlim([w_grid(1) w_grid(end)]);
% return
% 
% 
%     % Bisection 
%     w_low  = 0.4 ; 
%     w_high = 0.42 ; 
%     excess_supply_low  = Find_DBN_Histogram(w_low)  ; 
%     excess_supply_high = Find_DBN_Histogram(w_high) ;
%     disp('Starting Bisection for Equilibrium Wage')
%     disp({'iter','w_low','w','w_high','excess_supply','w_dist'})
%     if excess_supply_low*excess_supply_high<0 
%         excess_supply = 1 ;
%         iter = 0 ;
%         while abs(excess_supply)>5e-4 && iter<maxIterations
%             w = (w_low+w_high)/2 ;
%             excess_supply = Find_DBN_Histogram(w)  ; 
%             if excess_supply>0 
%                 w_high = w ; 
%             else 
%                 w_low  = w ;
%             end 
%             iter = iter + 1 ;
%             disp(num2cell([iter w_low w w_high excess_supply w_high-w_low]))
%         end 
%     else
%         disp('Equilibrium not bracketed')
%         return
%     end 

    options = optimset('Display','iter','TolFun',5e-04);
    [x,err] = fzero(@(x)Find_DBN_Histogram(x),0.42,options);
    
    r_exp = r    ;
    w_exp = x(1) ; 

    Mat = [r_exp w_exp] ; 
    Mat = [{'Aggregate Prices',' ';'Interest Rate','Wage'};num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    Mat = (AA*vZ_Grid.*(aalpha/(r_exp+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w_exp))^mmu).^(1/(1-aalpha-mmu)) ; 
    Mat = [{'Optimal Capital'} cell(1,n_Z-1);num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    % Get Distribution, Value and Policy Functions
    [price_residual_exp,mDBN_W_exp,mDBN_E_exp,mAp_W_exp,mAp_E_exp,OC_W_exp,OC_E_exp,V_W_exp,V_E_exp,Transition_exp] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual_exp);
    
    
    [A_exp,K_exp,N_exp,Y_exp,Earnings_W_exp,Earnings_E_exp] = ...
        Graphs_Tables('exp',r_exp,w_exp,mDBN_W_exp,mDBN_E_exp,mAp_W_exp,mAp_E_exp,OC_W_exp,OC_E_exp,V_W_exp,V_E_exp,Transition_exp) ;
    
    % Aggregate variables difference
    Mat = 100*[A_exp/A_ben-1 K_exp/K_ben-1  N_exp/N_ben-1  Y_exp/Y_ben-1];
    Mat = [{'% Change'} cell(1,3);{'A','K','N','Y'};num2cell(Mat)] ;
    disp(' '); disp(Mat); disp(' ');
    
    Mat = 100*[(1+r_exp)/(1+r_ben)-1 w_exp/w_ben-1] ;
    Mat = [{'% Change'} cell(1,1);{'Int. Rate','Wage'};num2cell(Mat)] ;
    disp(' '); disp(Mat); disp(' ');

    % Compute Welfare Gain
    Welfare_Gain(mDBN_W,mDBN_E,V_W,V_E,mDBN_W_exp,mDBN_E_exp,V_W_exp,V_E_exp,Earnings_W_ben,Earnings_E_ben)
    

    
    
%% Experiment 2: Change in labor taxes

    % Set back A Grid
    vA_Grid = vA_Grid_ben ;
    mA_Grid = repmat(vA_Grid,[1 n_E n_Z]);
    
    % Adjust Unemployment payment Grid
    disp('Set unemployment earnings to 20% of wage (Originally 5%)')
    vE_Grid(1) = 0.2 ;
    mE_Grid = repmat(vE_Grid',[n_A 1 n_Z]);

    % Solve the model
    t0 = tic;
    x_0 = x ;
%     options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
%     [x,err,exitflag] = fsolve(@(x) Find_DBN_Histogram(x),x_0,options);
%     options = optimset('Display','iter','TolFun',1e-05,'TolX',1e-06);
%     [x,err,exitflag] = fminsearch(@(x) Find_DBN_Histogram(x),x_0,options);
%     fprintf('Done! Time to compute: %2.2f seconds, error=%2.2d \n\n',toc(t0),err)
    
    options = optimset('Display','iter','TolFun',5e-04);
    [x,err] = fzero(@(x)Find_DBN_Histogram(x),0.42,options);

    r_exp2 = r    ;
    w_exp2 = x(1) ; 

    Mat = [r_exp2 w_exp2] ; 
    Mat = [{'Aggregate Prices',' ';'Interest Rate','Wage'};num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    Mat = (AA*vZ_Grid.*(aalpha/(r_exp2+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w_exp2))^mmu).^(1/(1-aalpha-mmu)) ; 
    Mat = [{'Optimal Capital'} cell(1,n_Z-1);num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    % Get Distribution, Value and Policy Functions
    [price_residual_exp2,mDBN_W_exp2,mDBN_E_exp2,mAp_W_exp2,mAp_E_exp2,OC_W_exp2,OC_E_exp2,V_W_exp2,V_E_exp2,Transition_exp2] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual_exp2);
    
    
    [A_exp2,K_exp2,N_exp2,Y_exp2,Earnings_W_exp2,Earnings_E_exp2] = ...
        Graphs_Tables('exp2',r_exp2,w_exp2,mDBN_W_exp2,mDBN_E_exp2,mAp_W_exp2,mAp_E_exp2,OC_W_exp2,OC_E_exp2,V_W_exp2,V_E_exp2,Transition_exp2) ;
    
    % Aggregate variables difference
    Mat = 100*[A_exp2/A_ben-1 K_exp2/K_ben-1  N_exp2/N_ben-1  Y_exp2/Y_ben-1];
    Mat = [{'% Change'} cell(1,3);{'A','K','N','Y'};num2cell(Mat)] ;
    disp(' '); disp(Mat); disp(' ');
    
    Mat = 100*[(1+r_exp2)/(1+r_ben)-1 w_exp2/w_ben-1] ;
    Mat = [{'% Change'} cell(1,1);{'Int. Rate','Wage'};num2cell(Mat)] ;
    disp(' '); disp(Mat); disp(' ');

    % Compute Welfare Gain
    Welfare_Gain(mDBN_W,mDBN_E,V_W,V_E,mDBN_W_exp2,mDBN_E_exp2,V_W_exp2,V_E_exp2,Earnings_W_ben,Earnings_E_ben)


%% 

diary off