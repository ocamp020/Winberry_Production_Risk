% Production Risk Economy with Occupational Choice
% Graphs And Tables
% Juan David Herreno & Sergio Ocampo (2016)

function [] = Welfare_Gain(mDBN_W,mDBN_E,V_W,V_E,mDBN_W_exp,mDBN_E_exp,V_W_exp,V_E_exp,Earnings_W_ben,Earnings_E_ben)

global n_A n_E n_Z vA_Grid vA_Grid_ben bbeta ggamma


    mDBN_W_exp(V_W_exp==-Inf) = 0 ;
    V_W_exp(V_W_exp==-Inf) = 0 ;
    mDBN_E_exp(V_E_exp==-Inf) = 0 ; 
    V_E_exp(V_E_exp==-Inf) = 0 ;

%% Welfare Gain by state
    % Approximate value functions
        V_W_aux = NaN(n_A,n_E,n_Z);
        V_E_aux = NaN(n_A, 1 ,n_Z);
        for i_z=1:n_Z           
            V_E_aux(:,1,i_z) = interp1(vA_Grid,V_E_exp(:,1,i_z),vA_Grid_ben) ;
        for i_e=1:n_E
            V_W_aux(:,i_e,i_z) = interp1(vA_Grid,V_W_exp(:,i_e,i_z),vA_Grid_ben) ;
        end 
        end
        
    % Compute Welfare Gain by state
        CE_W = 100*( exp((V_W_aux - V_W )*(1-bbeta*ggamma)) - 1 ) ;
        CE_E = 100*( exp((V_E_aux - V_E )*(1-bbeta*ggamma)) - 1 ) ;

    
%% Average Welfare Gains
    Av_CE   = sum(sum(sum( CE_W.*mDBN_W ))) + sum(sum(sum( CE_E.*mDBN_E ))) ;
    Av_CE_W = sum(sum(sum( CE_W(:,2:n_E,:).*mDBN_W(:,2:n_E,:) ))) / sum(sum(sum( mDBN_W(:,2:n_E,:) ))) ;
    Av_CE_U = sum(sum(sum( CE_W(:,1,:).*mDBN_W(:,1,:) ))) / sum(sum(sum( mDBN_W(:,1,:) ))) ;
    Av_CE_E = sum(sum(sum( CE_E.*mDBN_E ))) / sum(sum(sum( mDBN_E ))) ;
    
    Mat = [Av_CE_E Av_CE_W Av_CE_U Av_CE];
    Mat = [{'CE - SE','CE - W','CE - U','CE - All'};num2cell(Mat)] ;
    disp(' '); disp('Welfare Gain by state of benchmark'); disp(Mat); disp(' ');
 
    Mat = squeeze(sum( CE_W.*mDBN_W ,1))./squeeze(sum( mDBN_W ,1)) ;
    disp(' '); disp('Welfare Gain for Workers by EZ'); disp(Mat); disp(' ');
    
    Mat = squeeze(sum( CE_E.*mDBN_E ,1))./squeeze(sum( mDBN_E ,1)) ;
    disp(' '); disp('Welfare Gain for Self-Employed by EZ'); disp(Mat); disp(' ');
    
%% Welfare Gain of Average agent
    Av_V_W_exp = sum(sum(sum( V_W_exp(:,2:n_E,:).*mDBN_W_exp(:,2:n_E,:) ))) / sum(sum(sum(mDBN_W_exp(:,2:n_E,:)))) ;
    Av_V_W_ben = sum(sum(sum( V_W(:,2:n_E,:)    .*mDBN_W(:,2:n_E,:)     ))) / sum(sum(sum(mDBN_W(:,2:n_E,:))))     ;

    Av_V_U_exp = sum(sum(sum( V_W_exp(:,1,:).*mDBN_W_exp(:,1,:) ))) / sum(sum(sum(mDBN_W_exp(:,1,:)))) ;
    Av_V_U_ben = sum(sum(sum( V_W(:,1,:)    .*mDBN_W(:,1,:)     ))) / sum(sum(sum(mDBN_W(:,1,:))))     ;

    Av_V_E_exp = sum(sum(sum( V_E_exp.*mDBN_E_exp ))) / sum(sum(sum(mDBN_E_exp))) ;
    Av_V_E_ben = sum(sum(sum( V_E    .*mDBN_E     ))) / sum(sum(sum(mDBN_E))) ;

    Av_V_exp   = sum(sum(sum( V_W_exp.*mDBN_W_exp + V_E_exp.*mDBN_E_exp )))  ;
    Av_V_ben   = sum(sum(sum( V_W    .*mDBN_W     + V_E    .*mDBN_E     )))  ;

    Av_CE2   = 100*(exp((Av_V_exp   - Av_V_ben  )*(1-bbeta*ggamma))-1) ;
    Av_CE2_W = 100*(exp((Av_V_W_exp - Av_V_W_ben)*(1-bbeta*ggamma))-1) ;
    Av_CE2_U = 100*(exp((Av_V_U_exp - Av_V_U_ben)*(1-bbeta*ggamma))-1) ;
    Av_CE2_E = 100*(exp((Av_V_E_exp - Av_V_E_ben)*(1-bbeta*ggamma))-1) ;
    
    Mat = [Av_CE2_E Av_CE2_W Av_CE2_U Av_CE2];
    Mat = [{'CE - SE','CE - W','CE - U','CE - All'};num2cell(Mat)] ;
    disp(' '); disp('Welfare Gain by Averages'); disp(Mat); disp(' ');
    

%% Welfare Gain by Earning Percentile
    CE_W_aux = CE_W(:,2:n_E,:) ;
    CE_vec   = [ CE_W_aux(:) ; CE_E(:) ] ;
    % Earnings
    Earnings_vec = [Earnings_W_ben(:) ; Earnings_E_ben(:) ] ;
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
    file_name_eps = ['CE_by_pct.eps'] ;
    file_name_fig = ['CE_by_pct.fig'] ;
    print('-depsc',file_name_eps)
    savefig(file_name_fig)

end