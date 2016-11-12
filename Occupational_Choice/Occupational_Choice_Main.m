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
    x_0 = [0.01 0.23 0.34] ;
    x_0 = [0.010008959812246   0.217382233344189   0.379041511277377];
%     options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
%     [x,err,exitflag] = fsolve(@(x) Find_DBN_Histogram(x),x_0,options);
    options = optimset('Display','iter','TolFun',1e-03);
    [x,err,exitflag] = fminsearch(@(x) Find_DBN_Histogram(x),x_0,options);
    fprintf('Done! Time to compute: %2.2f seconds, error=%2.2d \n\n',toc(t0),err)

    r_ben = x(1) ;
    p_ben = x(2) ;
    w_ben = x(3) ; 
    
    Mat = [r_ben p_ben w_ben] ; 
    Mat = [{'Aggregate Prices',' ',' ';'Interest Rate','Int. Goods Price','Wage'};num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    Mat = (mmu*p_ben*vZ_Grid.^mmu/(r_ben+ddelta)).^(1/(1-mmu)) ; 
    Mat = [{'Optimal Capital'} cell(1,n_Z-1);num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

%% Distribution, Value and Policy Functions 

    [price_residual,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual);
    
    [A_ben,X_ben,N_ben,Y_ben,Earnings_W_ben,Earnings_E_ben] = ...
        Graphs_Tables('ben',r_ben,p_ben,w_ben,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E) ;
    

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
    options = optimset('Display','iter','TolFun',1e-03);
    [x,err,exitflag] = fminsearch(@(x) Find_DBN_Histogram(x),x_0,options);
    fprintf('Done! Time to compute: %2.2f seconds, error=%2.2d \n\n',toc(t0),err)

    r_exp = x(1) ;
    p_exp = x(2) ;
    w_exp = x(3) ; 

    Mat = [r_exp p_exp w_exp] ; 
    Mat = [{'Aggregate Prices',' ',' ';'Interest Rate','Int. Goods Price','Wage'};num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    Mat = (mmu*p_exp*vZ_Grid.^mmu/(r_exp+ddelta)).^(1/(1-mmu)) ; 
    Mat = [{'Optimal Capital'} cell(1,n_Z-1);num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

    % Get Distribution, Value and Policy Functions
    [price_residual_exp,mDBN_W_exp,mDBN_E_exp,mAp_W_exp,mAp_E_exp,OC_W_exp,OC_E_exp,V_W_exp,V_E_exp] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual_exp);
    
    
    [A_exp,X_exp,N_exp,Y_exp,Earnings_W_exp,Earnings_E_exp] = ...
        Graphs_Tables('exp',r_exp,p_exp,w_exp,mDBN_W_exp,mDBN_E_exp,mAp_W_exp,mAp_E_exp,OC_W_exp,OC_E_exp,V_W_exp,V_E_exp) ;
    
    % Aggregate variables difference
    Mat = 100*[A_exp/A_ben-1 X_exp/X_ben-1  N_exp/N_ben-1  Y_exp/Y_ben-1];
    Mat = [{'% Change'} cell(1,3);{'A','X','N','Y'};num2cell(Mat)] ;
    disp(' '); disp(Mat); disp(' ');
    
    Mat = 100*[(1+r_exp)/(1+r_ben)-1 p_exp/p_ben-1 w_exp/w_ben-1] ;
    Mat = [{'% Change'} cell(1,2);{'Int. Rate','Int. Good Price','Wage'};num2cell(Mat)] ;
    disp(' '); disp(Mat); disp(' ');

    % Compute Welfare Gain
    Welfare_Gain(mDBN_W,mDBN_E,V_W,V_E,mDBN_W_exp,mDBN_E_exp,V_W_exp,V_E_exp,Earnings_W_ben,Earnings_E_ben)
    

    
    
%% Experiment 2: Change in labor taxes



%% 

diary off