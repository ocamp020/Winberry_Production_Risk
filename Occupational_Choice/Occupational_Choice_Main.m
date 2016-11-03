% Production Risk Economy with Occupational Choice
% Main File
% Juan David Herreno & Sergio Ocampo (2016)

clear all; close all; clc;
%% Set up problem 

    Set_Parameters;

%% Stationary Equilibrium (Histogram Method and discrete VFI)

    t0 = tic;
    x_0 = [0.05 0.12 1] ;
%     options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
%     [x,err,exitflag] = fsolve(@(x) Find_DBN_Histogram(x),x_0,options);
    options = optimset('Display','iter','TolFun',1e-03);
    [x,err,exitflag] = fminsearch(@(x) Find_DBN_Histogram(x),x_0,options);
    fprintf('Done! Time to compute: %2.2f seconds, error=%2.2d \n\n',toc(t0),err)

    r_ss = x(1) ;
    p_ss = x(2) ;
    w_ss = x(3) ; 

    disp('Interest Rate')
    disp(r_ss)
    disp('Price')
    disp(p_ss)
    disp('Wage')
    disp(w_ss)
    disp('Optimal Capital')
    disp((mmu*p_ss*vZ_Grid.^mmu/(r_ss+ddelta)).^(1/(1-mmu)))


%% Distribution, Value and Policy Functions 

    [~,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E] = Find_DBN_Histogram(x) ;



%% Graphs 



%% Experiment 1: Change in household's borrowing constraint



%% Experiment 2: Change in labor taxes



