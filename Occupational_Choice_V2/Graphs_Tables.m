% Production Risk Economy with Occupational Choice
% Graphs And Tables
% Juan David Herreno & Sergio Ocampo (2016)

function [A_ss,X_ss,N_ss,Y_ss,Earnings_W,Earnings_E] = Graphs_Tables(model,r,w,mDBN_W,mDBN_E,mAp_W,mAp_E,OC_W,OC_E,V_W,V_E,Transition)

global mmu ddelta llambda tau_k tau_n AA aalpha...
       n_A n_E n_Z vA_Grid vE_Grid mA_Grid mZ_Grid mE_Grid

%% Capital and Profits and Earnings

    K  = min( max( 0 , llambda*mA_grid ) , ...
               (AA*mZ_Grid.*(aalpha/(r+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w))^mmu).^(1/(1-aalpha-mmu)) )  ;
    N  = (mmu*AA*mZ_Grid.*K.^aalpha/((1+tau_n)*w)).^(1/(1-mmu)) ;
    Pr = AA*mZ_Grid.*K.^aalpha.*N.^mmu - (r+ddelta)*K - (1+tau_n)*w*N;
    
    % Self-Employed Earnings
    Earnings_E   = r*(1-tau_k)*mA_Grid + (1-tau_k)*Pr ;
    Earnings_E_vec = Earnings_E(:); 
    [Earnings_E_vec,Earnings_E_ind] = sort(Earnings_E_vec) ;
    % Labor Earnings
    Earnings_W   = r*(1-tau_k)*mA_Grid(:,2:n_E,:) + w*mE_Grid(:,2:n_E,:) ; 
    Earnings_W_vec = Earnings_W(:); 
    [Earnings_W_vec,Earnings_W_ind] = sort(Earnings_W_vec) ;
    % Earnings
    Earnings_vec = [Earnings_W(:) ; Earnings_E(:) ] ;
    [Earnings_vec,Earnings_ind] = sort(Earnings_vec) ;
       
%% Policy and Value Functions

    % Value Functions
        figure; 
        for i_z=1:n_Z
            subplot(2,3,i_z); plot(vA_Grid,V_W(:,:,i_z),'linewidth',2); title('V_W');
            xlim([min(vA_Grid) max(vA_Grid)])
        end 
        set(gcf,'color','w')
        file_name_eps = ['V_W_',model,'.eps'] ;
        file_name_fig = ['V_W_',model,'.fig'] ;
        print('-depsc',file_name_eps)
        savefig(file_name_fig)

        figure; 
        for i_z=1:n_Z
            subplot(2,3,i_z); plot(vA_Grid,V_E(:,:,i_z),'linewidth',2); title('V_E');
            xlim([min(vA_Grid) max(vA_Grid)])
        end 
        set(gcf,'color','w')
        file_name_eps = ['V_E_',model,'.eps'] ;
        file_name_fig = ['V_E_',model,'.fig'] ;
        print('-depsc',file_name_eps)
        savefig(file_name_fig)
        
        
    % Saving Functions
        figure;
        i_plot  = 1;
        for i_e = 1:n_E
        for i_z = 1:n_Z
            subplot(n_E,n_Z,i_plot); hold on; 
            plot(vA_Grid,vA_Grid(mAp_W(:,i_e,i_z)),'linewidth',2)
            plot(vA_Grid,vA_Grid,':','linewidth',2)
            hold off; xlim([min(vA_Grid) max(vA_Grid)])
            %title('Saving W')
            i_plot = i_plot +1 ;
        end
        end
        set(gcf,'color','w')
        file_name_eps = ['Ap_W_',model,'.eps'] ;
        file_name_fig = ['Ap_W_',model,'.fig'] ;
        print('-depsc',file_name_eps)
        savefig(file_name_fig)

        figure;
        i_plot  = 1;
        for i_e = 1:n_E
        for i_z = 1:n_Z
            subplot(n_E,n_Z,i_plot); hold on; 
            plot(vA_Grid,vA_Grid(mAp_E(:,i_e,i_z)),'linewidth',2)
            plot(vA_Grid,vA_Grid,':')
            hold off; xlim([min(vA_Grid) max(vA_Grid)])
            %title('Saving E')
            i_plot = i_plot +1 ;
        end
        end
        set(gcf,'color','w')
        file_name_eps = ['Ap_E_',model,'.eps'] ;
        file_name_fig = ['Ap_E_',model,'.fig'] ;
        print('-depsc',file_name_eps)
        savefig(file_name_fig)
        
    % Occupational Choice
        figure;
        i_plot  = 1;
        for i_e = 1:n_E
        for i_z = 1:n_Z
            subplot(n_E,n_Z,i_plot);
            plot(vA_Grid,OC_W(:,i_e,i_z),'linewidth',2)
            xlim([min(vA_Grid) max(vA_Grid)]); ylim([-1,2]);
            %title('OC W')
            i_plot = i_plot +1 ;
        end
        end
        set(gcf,'color','w')
        file_name_eps = ['OC_W_',model,'.eps'] ;
        file_name_fig = ['OC_W_',model,'.fig'] ;
        print('-depsc',file_name_eps)
        savefig(file_name_fig)

        figure;
        i_plot  = 1;
        for i_e = 1:n_E
        for i_z = 1:n_Z
            subplot(n_E,n_Z,i_plot);
            plot(vA_Grid,OC_E(:,i_e,i_z),'linewidth',2)
            xlim([min(vA_Grid) max(vA_Grid)]); ylim([-1,2]);
            %title('OC E')
            i_plot = i_plot +1 ;
        end
        end
        set(gcf,'color','w')
        file_name_eps = ['OC_E_',model,'.eps'] ;
        file_name_fig = ['OC_E_',model,'.fig'] ;
        print('-depsc',file_name_eps)
        savefig(file_name_fig)


%% Distribution and Composition
    figure; 
    plot(vA_Grid,[sum(sum(mDBN_W,3),2)/sum(sum(sum(mDBN_W))),sum(sum(mDBN_E,3),2)/sum(sum(sum(mDBN_E))),sum(sum(mDBN_W+mDBN_E,3),2)],'-')
    xlabel('Wealth')
    title('Wealth Distribution'); legend('Workers','Self-Employed','All','location','southoutside','orientation','horizontal')
    set(gcf,'color','w')
    file_name_eps = ['Wealth_Dist_',model,'.eps'] ;
    file_name_fig = ['Wealth_Dist_',model,'.fig'] ;
    print('-depsc',file_name_eps)
    savefig(file_name_fig)
    
     % Table with Composition
        % Workers
        W_share    = sum(sum(sum(mDBN_W))) ;
        W_share_U  = sum(sum(sum(mDBN_W(:,1,:))));
        W_share_EZ = squeeze(sum(mDBN_W,1))/W_share ;
        % Entrepreneurs
        E_share    = sum(sum(sum(mDBN_E))) ;
        E_share_EZ = squeeze(sum(mDBN_E,1))/E_share ;
        EZ_share_E = squeeze(sum(mDBN_E,1))./squeeze(sum(mDBN_W+mDBN_E,1));
        
        Mat = 100*[W_share-W_share_U E_share W_share_U] ; 
        Mat = [{'Agent Composition',' ',' ';'Workers','Self-Employed','Unemployed'};num2cell(Mat)];
        disp(' '); disp(Mat); disp(' ');
        
        Mat = 100*EZ_share_E ; 
        Mat = [{'Self-Employed by EZ'} cell(1,n_Z-1);num2cell(Mat)];
        disp(' '); disp(Mat); disp(' ');
        
        Mat = 100*W_share_EZ ; 
        Mat = [{'Composition Workers'} cell(1,n_Z-1);num2cell(Mat)];
        disp(' '); disp(Mat); disp(' ');
        
        
        Mat = 100*E_share_EZ ; 
        Mat = [{'Composition Self-Employed'} cell(1,n_Z-1);num2cell(Mat)];
        disp(' '); disp(Mat); disp(' ');

%% Transitions
    Mat = Transition ; 
    Mat = [{' ','Up','Wp','Ep'};{'U';'W';'E'} num2cell(Mat)] ;
    disp(' '); disp('Transitions'); disp(Mat); disp(' ');


%% Top Wealth Shares

    DBN_A   = sum(sum(mDBN_W+mDBN_E,3),2) ;
    C_DBN_A = cumsum(DBN_A) ; 
    pct_ind = knnsearch(C_DBN_A,[0.1 0.25 0.5 0.75 0.9 0.99]');
    pct_A   = vA_Grid(pct_ind) ;
    for i=1:numel(pct_ind)
        top_shares(i) = sum(vA_Grid(pct_ind(i):end).*DBN_A(pct_ind(i):end))/(vA_Grid'*DBN_A) ;
    end
    
    Mat = 100*top_shares ; 
    Mat = [{'Top Wealth Shares'} cell(1,numel(pct_ind));
           {'pct'},num2cell(100*[0.1 0.25 0.5 0.75 0.9 0.99]);
           {'Share'},num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');

%% Earnings

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
    pct_list = [0.05:0.05:0.95]' ;
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
    figure; plot(pct_list,E_share_pct,'-o','linewidth',2); title('Share of self-employed by earnings pct'); 
    set(gcf,'color','w')
    file_name_eps = ['SE_Rate_by_pct_',model,'.eps'] ;
    file_name_fig = ['SE_Rate_by_pct_',model,'.fig'] ;
    print('-depsc',file_name_eps)
    savefig(file_name_fig)
    
    
    % Distribution of Earnings 
    DBN_all_E_vec = DBN_all_E(:); DBN_all_E_vec = DBN_all_E_vec(Earnings_E_ind); DBN_all_E_vec = DBN_all_E_vec/sum(DBN_all_E_vec);
    DBN_all_W_vec = DBN_all_W(:); DBN_all_W_vec = DBN_all_W_vec(Earnings_W_ind); DBN_all_W_vec = DBN_all_W_vec/sum(DBN_all_W_vec);
    figure;
    plot(Earnings_vec,DBN_all_vec,Earnings_E_vec,DBN_all_E_vec,Earnings_W_vec,DBN_all_W_vec,'linewidth',2)
    legend('All Agents','Workers','Self-Employed','location','southoutside','orientation','horizontal')
    xlim([min(Earnings_vec),max(Earnings_vec)]);
    set(gcf,'color','w')
    file_name_eps = ['Earnings_Dist_',model,'.eps'] ;
    file_name_fig = ['Earnings_Dist_',model,'.fig'] ;
    print('-depsc',file_name_eps)
    savefig(file_name_fig)
    


%% Aggregates

    A_ss = vA_Grid'*sum(sum(mDBN_W+mDBN_E,3),2) ;
    K_ss = sum(sum(sum( K.*mDBN_E )))  ;
    N_ss = sum(sum(sum( N.*mDBN_E )))  ;
    Y_ss = sum(sum(sum( AA*mZ_Grid.*K^(aalpha)*N^mmu.*mDBN_E ))) ;   
    
    Mat = [A_ss K_ss N_ss Y_ss] ; 
    Mat = [{'Aggregate Variables',' ',' ',' ';'Assets','Capital','Labor','Output'};num2cell(Mat)];
    disp(' '); disp(Mat); disp(' ');
    


end