% Production Risk Economy with Occupational Choice
% Main File
% Juan David Herreno & Sergio Ocampo (2016)

clear all; close all; clc;
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

    [price_residual,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual);
    
    

    DBN_A   = sum(sum(mDBN_W+mDBN_E,3),2) ;
    C_DBN_A = cumsum(DBN_A) ; 
    pct_ind = knnsearch(C_DBN_A,[0.1 0.25 0.5 0.75 0.9 0.99]');
    pct_A   = vA_Grid(pct_ind) ;
    for i=1:numel(pct_ind)
        top_shares(i) = sum(vA_Grid(pct_ind(i):end).*DBN_A(pct_ind(i):end))/(vA_Grid'*DBN_A) ;
    end
    
    disp('Top Wealth Shares')
    disp(100*[0.1 0.25 0.5 0.75 0.9 0.99])
    disp(100*top_shares)

    
    % Table with Composition
        % Workers
        W_share    = sum(sum(sum(mDBN_W))) ;
        W_share_U  = sum(sum(sum(mDBN_W(:,1,:))));
        W_share_EZ = squeeze(sum(mDBN_W,1))/W_share ;
        % Entrepreneurs
        E_share    = sum(sum(sum(mDBN_E))) ;
        E_share_EZ = squeeze(sum(mDBN_E,1))/E_share ;
        EZ_share_E = squeeze(sum(mDBN_E,1))./squeeze(sum(mDBN_W+mDBN_E,1));
        
        disp(['Workers','Self-Employed','Unemployed'])
        disp(100*[W_share-W_share_U E_share W_share_U])
        disp(['Self-Employed by EZ'])
        disp(100*EZ_share_E)
        disp('Composition of workers')
        disp(100*W_share_EZ)
        disp('Composition of self-employed')
        disp(100*E_share_EZ)
        
        
	% Earnings 
        % Firm's capital demand and profits
        K            = min( llambda*mA_Grid , (mmu*p_ss*mZ_Grid.^mmu/(r_ss+ddelta)).^(1/(1-mmu)) )  ;
        Pr           = p_ss*(mZ_Grid.*K).^mmu - (r_ss+ddelta)*K ;
        Earnings_E   = r_ss*(1-tau_k)*mA_Grid + (1-tau_k)*Pr ;
        % Labor Earnings
        Earnings_W   = r_ss*(1-tau_k)*mA_Grid(:,2:n_E,:) + w_ss*mE_Grid(:,2:n_E,:) ; 
        % Earnings
        Earnings_vec = [Earnings_W(:) ; Earnings_E(:) ] ;
        [Earnings_vec,Earnings_ind] = sort(Earnings_vec) ;
        % Occupation indicator
        Occupation_W = zeros(n_A,n_E-1,n_Z) ;
        Occupation_E = ones(n_A,n_E,n_Z)    ;
        Occupation   = [Occupation_W(:) ; Occupation_E(:)] ;
        Occupation   = Occupation(Earnings_ind) ;
        % DBN_all
        DBN_all_W    = mDBN_W(:,2:n_E,:) ; 
        DBN_all_E    = mDBN_E ; 
        DBN_all_vec  = [DBN_all_W(:) ; DBN_all_E(:)]  ;
        DBN_all_vec  = DBN_all_vec(Earnings_ind)/sum(DBN_all_vec) ;
        C_DBN_all_vec= cumsum(DBN_all_vec) ;
        % Percentiles 
        pct_list = [0.05:0.05:1]' ;
        pct_ind = knnsearch(C_DBN_all_vec,pct_list);
        E_share_pct = NaN(numel(pct_list),1) ;
        low = 0 ;
        for i = 1:numel(pct_list)
            high = Earnings_vec(pct_ind(i)) ;
            ind = Earnings_vec>low & Earnings_vec<=high ;
            mass_E   = sum( DBN_all_vec.*Occupation.*ind ) ; 
            mass_tot = sum( DBN_all_vec.*ind ) ;
            E_share_pct(i) = mass_E/mass_tot ; 
            low = high ;
        end 
        figure; plot(pct_list,E_share_pct); title('Share of self-employed by earnings pct'); 
        set(gcf,'color','w')
        print('-depsc','SS_SE_Rate_by_pct.eps')

    % Saving Functions
    figure;
    i_plot  = 1;
    for i_e = 1:n_E
    for i_z = 1:n_Z
        subplot(n_E,n_Z,i_plot);
        hold on; 
        plot(vA_Grid,vA_Grid(mAp_W(:,i_e,i_z)))
        plot(vA_Grid,vA_Grid,':')
        hold off;
        %title('Saving W')
        i_plot = i_plot +1 ;
    end
    end
    set(gcf,'color','w')
    
    figure;
    i_plot  = 1;
    for i_e = 1:n_E
    for i_z = 1:n_Z
        subplot(n_E,n_Z,i_plot);
        hold on; 
        plot(vA_Grid,vA_Grid(mAp_E(:,i_e,i_z)))
        plot(vA_Grid,vA_Grid,':')
        hold off;
        %title('Saving E')
        i_plot = i_plot +1 ;
    end
    end
    set(gcf,'color','w')
        
    % Aggregates (output, capital, etc)
    A_ss = vA_Grid'*DBN_A ;
    K_ss = min( llambda*mA_Grid , (mmu*p_ss*mZ_Grid.^mmu/(r_ss+ddelta)).^(1/(1-mmu)) )  ;
    X_ss = sum(sum(sum( (mZ_Grid.*K_ss).^mmu .* mDBN_E ))).^(1/mmu) ;
    N_ss = sum( vE_Grid.*squeeze(sum(sum(mDBN_W,3),1))' ) ;
    Y_ss = AA*X_ss^(aalpha)*N_ss^(1-aalpha) ;   
    
    % Transitions
    
    
    % Distribution Graph
    figure; 
    plot(vA_Grid,[sum(sum(mDBN_W,3),2)/sum(sum(sum(mDBN_W))),sum(sum(mDBN_E,3),2)/sum(sum(sum(mDBN_E))),sum(sum(mDBN_W+mDBN_E,3),2)],'-')
    xlabel('Wealth')
    title('Wealth Distribution'); legend('Workers','Entrepreneurs','All','location','southoutside','orientation','horizontal')
    set(gcf,'color','w')
    print('-depsc','SS_Distribution_Wealth.eps')

    



%% Experiment 1: Change in household's borrowing constraint

vA_Grid_old = vA_Grid ;

aaBar = -0.1*w_ss ;
% Bounds on grid space
A_Min = aaBar;	
A_Max = 50  ;
Grid_Curvature = 3.0 ;
vA_Grid = A_Min + linspace(0,1,n_A)'.^Grid_Curvature.*(A_Max-A_Min) ;
mA_Grid = repmat(vA_Grid,[1 n_E n_Z]);




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

    disp('Interest Rate')
    disp(r_exp)
    disp('Price')
    disp(p_exp)
    disp('Wage')
    disp(w_exp)
    disp('Optimal Capital')
    disp((mmu*p_exp*vZ_Grid.^mmu/(r_exp+ddelta)).^(1/(1-mmu)))


[price_residual_exp,mDBN_W_exp,mDBN_E_exp,mAp_W_exp,mAp_E_exp,OC_W_exp,OC_E_exp,V_W_exp,V_E_exp] = Find_DBN_Histogram(x) ;
    disp('price residual'); disp(price_residual_exp);
    
    

    DBN_A_exp   = sum(sum(mDBN_W_exp+mDBN_E_exp,3),2) ;
    C_DBN_A = cumsum(DBN_A_exp) ; 
    pct_ind = knnsearch(C_DBN_A,[0.1 0.25 0.5 0.75 0.9 0.99]');
    pct_A   = vA_Grid(pct_ind) ;
    for i=1:numel(pct_ind)
        top_shares_exp(i) = sum(vA_Grid(pct_ind(i):end).*DBN_A_exp(pct_ind(i):end))/(vA_Grid'*DBN_A_exp) ;
    end
    
    disp('Top Wealth Shares')
    disp(100*[0.1 0.25 0.5 0.75 0.9 0.99])
    disp(100*top_shares_exp)

    
    % Table with Composition
        % Workers
        W_share_exp    = sum(sum(sum(mDBN_W_exp))) ;
        W_share_U_exp  = sum(sum(sum(mDBN_W_exp(:,1,:))));
        W_share_EZ_exp = squeeze(sum(mDBN_W_exp,1))/W_share_exp ;
        % Entrepreneurs
        E_share_exp    = sum(sum(sum(mDBN_E_exp))) ;
        E_share_EZ_exp = squeeze(sum(mDBN_E_exp,1))/E_share_exp ;
        EZ_share_E_exp = squeeze(sum(mDBN_E_exp,1))./squeeze(sum(mDBN_W_exp+mDBN_E_exp,1));
        
        disp(['Workers','Self-Employed','Unemployed'])
        disp(100*[W_share_exp-W_share_U_exp E_share_exp W_share_U_exp])
        disp(['Self-Employed by EZ'])
        disp(100*EZ_share_E_exp)
        disp('Composition of workers')
        disp(100*W_share_EZ_exp)
        disp('Composition of self-employed')
        disp(100*E_share_EZ_exp)
        
        
	% Earnings 
        % Firm's capital demand and profits
        K            = min( llambda*mA_Grid , (mmu*p_exp*mZ_Grid.^mmu/(r_exp+ddelta)).^(1/(1-mmu)) )  ;
        Pr           = p_exp*(mZ_Grid.*K).^mmu - (r_exp+ddelta)*K ;
        Earnings_E_exp   = r_exp*(1-tau_k)*mA_Grid + (1-tau_k)*Pr ;
        % Labor Earnings
        Earnings_W_exp   = r_exp*(1-tau_k)*mA_Grid(:,2:n_E,:) + w_exp*mE_Grid(:,2:n_E,:) ; 
        % Earnings
        Earnings_vec = [Earnings_W_exp(:) ; Earnings_E_exp(:) ] ;
        [Earnings_vec,Earnings_ind] = sort(Earnings_vec) ;
        % Occupation indicator
        Occupation_W = zeros(n_A,n_E-1,n_Z) ;
        Occupation_E = ones(n_A,n_E,n_Z)    ;
        Occupation   = [Occupation_W(:) ; Occupation_E(:)] ;
        Occupation   = Occupation(Earnings_ind) ;
        % DBN_all
        DBN_all_W    = mDBN_W_exp(:,2:n_E,:) ; 
        DBN_all_E    = mDBN_E_exp ; 
        DBN_all_vec  = [DBN_all_W(:) ; DBN_all_E(:)]  ;
        DBN_all_vec  = DBN_all_vec(Earnings_ind)/sum(DBN_all_vec) ;
        C_DBN_all_vec= cumsum(DBN_all_vec) ;
        % Percentiles 
        pct_list = [0.05:0.05:1]' ;
        pct_ind = knnsearch(C_DBN_all_vec,pct_list);
        E_share_pct_exp = NaN(numel(pct_list),1) ;
        low = 0 ;
        for i = 1:numel(pct_list)
            high = Earnings_vec(pct_ind(i)) ;
            ind = Earnings_vec>low & Earnings_vec<=high ;
            mass_E   = sum( DBN_all_vec.*Occupation.*ind ) ; 
            mass_tot = sum( DBN_all_vec.*ind ) ;
            E_share_pct_exp(i) = mass_E/mass_tot ; 
            low = high ;
        end 
        figure; plot(pct_list,E_share_pct_exp); title('Share of self-employed by earnings pct'); 
        set(gcf,'color','w')
        print('-depsc','SS_SE_Rate_by_pct_exp.eps')

    % Saving Functions
    figure;
    i_plot  = 1;
    for i_e = 1:n_E
    for i_z = 1:n_Z
        subplot(n_E,n_Z,i_plot);
        hold on; 
        plot(vA_Grid,vA_Grid(mAp_W_exp(:,i_e,i_z)))
        plot(vA_Grid,vA_Grid,':')
        hold off;
        %title('Saving W')
        i_plot = i_plot +1 ;
    end
    end
    set(gcf,'color','w')
    
    figure;
    i_plot  = 1;
    for i_e = 1:n_E
    for i_z = 1:n_Z
        subplot(n_E,n_Z,i_plot);
        hold on; 
        plot(vA_Grid,vA_Grid(mAp_E_exp(:,i_e,i_z)))
        plot(vA_Grid,vA_Grid,':')
        hold off;
        %title('Saving E')
        i_plot = i_plot +1 ;
    end
    end
    set(gcf,'color','w')
        
    % Aggregates (output, capital, etc)
    A_exp = vA_Grid'*DBN_A_exp ;
    K_exp = min( llambda*mA_Grid , (mmu*p_exp*mZ_Grid.^mmu/(r_exp+ddelta)).^(1/(1-mmu)) )  ;
    X_exp = sum(sum(sum( (mZ_Grid.*K_exp).^mmu .* mDBN_E_exp ))).^(1/mmu) ;
    N_exp = sum( vE_Grid.*squeeze(sum(sum(mDBN_W_exp,3),1))' ) ;
    Y_exp = AA*X_exp^(aalpha)*N_exp^(1-aalpha) ;  
    
    
        
    figure; 
    plot(vA_Grid,[sum(sum(mDBN_W_exp,3),2)/sum(sum(sum(mDBN_W_exp))),sum(sum(mDBN_E_exp,3),2)/sum(sum(sum(mDBN_E_exp))),sum(sum(mDBN_W_exp+mDBN_E_exp,3),2)],'-')
    xlabel('Wealth')
    title('Wealth Distribution'); legend('Workers','Entrepreneurs','All','location','southoutside','orientation','horizontal')
    set(gcf,'color','w')
    print('-depsc','SS_Distribution_Wealth_exp.eps')
    
    % Welfare Gain by state
        % Approximate value functions
        V_W_aux = NaN(n_A,n_E,n_Z);
        V_E_aux = NaN(n_A,n_E,n_Z);
        for i_z=1:n_Z
        for i_e=1:n_E
            V_W_aux(:,i_e,i_z) = interp1(vA_Grid,V_W_exp(:,i_e,i_z),vA_Grid_old) ;
            V_E_aux(:,i_e,i_z) = interp1(vA_Grid,V_E_exp(:,i_e,i_z),vA_Grid_old) ;
        end 
        end
    CE_W = 100*( exp((V_W_aux - V_W )/(1-bbeta*ggamma)) - 1 ) ;
    CE_E = 100*( exp((V_E_aux - V_E )/(1-bbeta*ggamma)) - 1 ) ;
    CE_W(mDBN_W<=1e-5) = 0 ;
    CE_E(mDBN_E<=1e-5) = 0 ;
    
    
    % Average Welfare Gains
        Av_CE   = sum(sum(sum( CE_W.*mDBN_W + CE_E.*mDBN_E ))) ;
        Av_CE_W = sum(sum(sum( CE_W(:,2:n_E,:).*mDBN_W(:,2:n_E,:) ))) / sum(sum(sum( mDBN_W(:,2:n_E,:) ))) ;
        Av_CE_U = sum(sum(sum( CE_W(:,1,:).*mDBN_W(:,1,:) ))) / sum(sum(sum( mDBN_W(:,1,:) ))) ;
        Av_CE_E = sum(sum(sum( CE_E.*mDBN_E ))) / sum(sum(sum( mDBN_E ))) ;
        
        mDBN_E_exp(V_E_exp==-Inf) = 0;
        V_E_exp(V_E_exp==-Inf) = 0 ;
        
        mDBN_E_exp = mDBN_E_exp/sum(mDBN_E_exp(:)) ;
        V_E_exp(mDBN_E_exp<1e-4) = 0 ;
        V_E_exp(V_E_exp==-Inf)= 0 ; 
        
        Av_V_W_exp = sum(sum(sum( V_W_exp.*mDBN_W_exp ))) / sum(sum(sum(mDBN_W_exp))) ;
        Av_V_W     = sum(sum(sum( V_W    .*mDBN_W     ))) / sum(sum(sum(mDBN_W))) ;
        
        Av_CE2   = 100*(exp((sum(sum(sum(V_W_exp.*mDBN_W_exp + V_E_exp.*mDBN_E_exp)))-sum(sum(sum(V_W.*mDBN_W + V_E.*mDBN_E))))/(1-bbeta*ggamma) )-1) ;
        Av_CE2_W = sum(sum(sum( CE_W(:,2:n_E,:).*mDBN_W(:,2:n_E,:) ))) / sum(sum(sum( mDBN_W(:,2:n_E,:) ))) ;
        Av_CE2_U = sum(sum(sum( CE_W(:,1,:).*mDBN_W(:,1,:) ))) / sum(sum(sum( mDBN_W(:,1,:) ))) ;
        Av_CE2_E = sum(sum(sum( CE_E.*mDBN_E ))) / sum(sum(sum( mDBN_E ))) ;
        
    % Welfare Gain by Earning Percentile
        CE_W_aux = CE_W(:,2:n_E,:) ;
        CE_vec = [ CE_W_aux(:) ; CE_E(:) ] ;
        % Earnings
        Earnings_vec = [Earnings_W(:) ; Earnings_E(:) ] ;
        [Earnings_vec,Earnings_ind] = sort(Earnings_vec) ;
        % DBN_all
        DBN_all_W    = mDBN_W(:,2:n_E,:) ; 
        DBN_all_E    = mDBN_E ; 
        DBN_all_vec  = [DBN_all_W(:) ; DBN_all_E(:)]  ;
        DBN_all_vec  = DBN_all_vec(Earnings_ind)/sum(DBN_all_vec) ;
        C_DBN_all_vec= cumsum(DBN_all_vec) ;
        % Percentiles 
        pct_list = [0.05:0.05:1]' ;
        pct_ind = knnsearch(C_DBN_all_vec,pct_list);
        CE_pct = NaN(numel(pct_list),1) ;
        low = 0 ;
        for i = 1:numel(pct_list)
            high = Earnings_vec(pct_ind(i)) ;
            ind = Earnings_vec>low & Earnings_vec<=high ;
            mass_E   = sum( DBN_all_vec.*CE_vec.*ind ) ; 
            mass_tot = sum( DBN_all_vec.*ind ) ;
            CE_pct(i) = mass_E/mass_tot ; 
            low = high ;
        end 
        figure; plot(pct_list,CE_pct); title('CE by earnings pct'); 
        set(gcf,'color','w')
        print('-depsc','SS_CE_by_pct.eps')

        % Tables
        Mat = [E_share          W_share-W_share_U        W_share_U     ;
               E_share_exp   W_share_exp-W_share_U_exp   W_share_U_exp ];
        Mat = 100*Mat ;
        Mat = [{'SE Share','W Share','U Share'};num2cell(Mat)] ;
        disp(Mat)
        
        
        Mat = [X_exp/X_ss-1  N_exp/N_ss-1  Y_exp/Y_ss-1 (1+r_exp)/(1+r_ss)-1 p_exp/p_ss-1 w_exp/w_ss-1];
        Mat = 100*Mat ;
        Mat = [{'X','N','Y','r','p','w'};num2cell(Mat)] ;
        disp(Mat)
        
        Mat = [Av_CE_E Av_CE_W Av_CE_U Av_CE];
        Mat = [{'CE - SE','CE - W','CE - U','CE - All'};num2cell(Mat)] ;
        disp(Mat)
    
    
%% Experiment 2: Change in labor taxes



