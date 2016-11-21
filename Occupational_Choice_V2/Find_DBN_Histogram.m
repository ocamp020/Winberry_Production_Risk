% Production Risk Economy with Occupational Choice
% Find Stationary Distribution with the Histogram Method
% Juan David Herreno & Sergio Ocampo (2016)

function [residual,mDBN_W_out,mDBN_E_out,mAp_W_out,mAp_E_out,OC_W_out,OC_E_out,VW_out,VE_out,mTransition_out] = Find_DBN_Histogram(x,solver)


%% Initialization (global variables)
    global ggamma aalpha ddelta mmu llambda AA tau_n mKappa r jfr_E ...
    	   mZ_Transition mE_Transition_W mE_Transition_E vE_Grid  ...
           n_E n_Z n_A vA_Grid mA_Grid_E mA_Grid_W mZ_Grid mZ_Ind_W n_State

%% Declare prices
    w = x(1);
        
    
%% Agent's value and policy functions (by discrete VFI)
    [Ap_W_VFI,Ap_E_VFI,OC_W_VFI,OC_E1_VFI,OC_E2_VFI,VW_VFI,VE_VFI] = Discrete_VFI(r,w) ;
    % disp(VW_VFI(:)) 
%     figure; plot(vA_Grid,[VW_VFI(:,:,1),VW_VFI(:,:,2),VW_VFI(:,:,3)]); title('VW');
%     figure; plot(vA_Grid,[VE_VFI(:,:,1),VE_VFI(:,:,2),VE_VFI(:,:,3)]); title('VE');
%     figure; plot(vA_Grid,[Ap_W_VFI(:,:,1),Ap_W_VFI(:,:,2),Ap_W_VFI(:,:,3)]); title('Ap_W');
%     figure; plot(vA_Grid,[Ap_E_VFI(:,:,1),Ap_E_VFI(:,:,2),Ap_E_VFI(:,:,3)]); title('Ap_E');
%     figure; plot(vA_Grid,[OC_W_VFI(:,:,1),OC_W_VFI(:,:,2),OC_W_VFI(:,:,3)]); title('Ap_W');
%     figure; plot(vA_Grid,[OC_E_VFI(:,:,1),OC_E_VFI(:,:,2),OC_E_VFI(:,:,3)]); title('Ap_E');
    
    
%% Stationary distribution (Histogram)
    
    % Vectorize saving function for both types of agent (W,E)
    vS_W = Ap_W_VFI(:) ; % In Index
    vS_E = Ap_E_VFI(:) ; % In Index
    
    % Assets transition matrix (Before occupational choice)
    mA_W_Transition = zeros(n_Z*n_E*n_A,n_A);
    mA_E_Transition = zeros(n_Z* 1 *n_A,n_A);
    for i_a = 1:n_A
        mA_W_Transition(vS_W == i_a,i_a) = 1 ;
        mA_E_Transition(vS_E == i_a,i_a) = 1 ;
    end
    mA_Transition = [repmat(mA_W_Transition, [1,n_E*n_Z]) , zeros(n_State,n_A*n_Z) ;
                     zeros(n_A*n_Z,n_State) , repmat(mA_E_Transition, [1,n_Z])];
    mA_Transition = sparse(mA_Transition) ;
    
    % Epsilon
    mE_Transition_Histogram = [ repmat( kron(mE_Transition_W,ones(n_A)) , [n_Z,n_Z] ) , zeros(n_State,n_A*n_Z) ; 
                                zeros(n_A*n_Z,n_State) , ones(n_A*n_Z) ];
    mE_Transition_Histogram = sparse(mE_Transition_Histogram) ;

    % Z
    mZ_Transition_Histogram = [kron( ggamma*eye(n_Z) + (1-ggamma)*mZ_Transition , ones(n_A*n_E) ) , zeros(n_State,n_A*n_Z) ;
                                zeros(n_A*n_Z,n_State)  , kron( ggamma*eye(n_Z) + (1-ggamma)*mZ_Transition , ones(n_A) )];
    mZ_Transition_Histogram = sparse(mZ_Transition_Histogram) ;
    
    % Full transition matrix (Before occupational choice)
    mTransition = mA_Transition .* mE_Transition_Histogram .* mZ_Transition_Histogram ;
    if abs(min(sum(mTransition,2))-1)>1e-10 || abs(max(sum(mTransition,2))-1)>1e-10
        disp([min(sum(mTransition,2)) max(sum(mTransition,2))])
        error('mTransition not working - Before occupational choice')
    end 
    % Ocupational Choice
    OC_W = OC_W_VFI(:) ; 
    OC_E1 = OC_E1_VFI(:) ; 
    OC_E2 = OC_E2_VFI(:) ; 
    
        % Adjust savings choice by switching cost
        vS_W_cost = vA_Grid(vS_W)-mKappa(:) ;

        % Find index and value of nearest neighbors of adjusted savings in A_grid
        vInd = knnsearch(vA_Grid,vS_W_cost);
        vVal = vA_Grid(vInd);

        vIndBelow = max( vInd - (vVal > vS_W_cost) , 1   ) ; vValBelow = vA_Grid(vIndBelow);
        vIndAbove = min( vIndBelow + 1             , n_A ) ; vValAbove = vA_Grid(vIndAbove);

        % Compute lottery matrices
        vPrBelow = min( max( (vValAbove - vS_W_cost) ./ (vValAbove - vValBelow) , 0 ) , 1 );
        vPrAbove = min( max( (vS_W_cost - vValBelow) ./ (vValAbove - vValBelow) , 0 ) , 1 );
        vPrBelow(vIndBelow==vIndAbove) = 1 ; 
        vPrAbove(vIndBelow==vIndAbove) = 0 ; 
        
        % Adjust mTransition for type W agents
        aux = [1:n_State];
        aux = aux(OC_W==1);
        z_ind = mZ_Ind_W(:) ;
        for i = aux
            mTransition(1:n_State,n_State + vIndBelow(i)+n_A*(z_ind(i)-1)) = ...
                mTransition(1:n_State,n_State + vIndBelow(i)+n_A*(z_ind(i)-1)) + ...
                vPrBelow(i)*mTransition(1:n_State,i) ;
            mTransition(1:n_State,n_State + vIndAbove(i)+n_A*(z_ind(i)-1)) = ...
                mTransition(1:n_State,n_State + vIndAbove(i)+n_A*(z_ind(i)-1)) + ...
                vPrAbove(i)*mTransition(1:n_State,i) ;
            mTransition(1:n_State,i) = 0 ;     
        end
        if abs(min(sum(mTransition,2))-1)>1e-10 || abs(max(sum(mTransition,2))-1)>1e-10
            disp([min(sum(mTransition,2)) max(sum(mTransition,2))])
            error('mTransition not working - After worker occupational choice')
        end 
                                                
        % Adjust mTransition for type E agents
        i=1;
        %stay_prob = ones(1,n_A*n_Z) ;
        mTransition_aux = mTransition;
        for i_z = 1:n_Z
        for i_ep = 1:n_E
        for i_a = 1:n_A
            if i_ep ==1 
                % Move to unemployment 
                mTransition(n_State+1:end,i) = mTransition(n_State+1:end,i) + ...
                    (1-jfr_E)*OC_E1_VFI(i_a,1,i_z)*mTransition_aux(n_State+1:end,i_a+n_A*(i_z-1));
                mTransition(n_State+1:end,i_a+n_A*(i_z-1)) = mTransition(n_State+1:end,i_a+n_A*(i_z-1)) - ...
                    (1-jfr_E)*OC_E1_VFI(i_a,1,i_z)*mTransition_aux(n_State+1:end,i_a+n_A*(i_z-1)); 
                %disp([i_ep stay_prob(i_a+n_A*(i_z-1)) (1-jfr_E)*OC_E1_VFI(i_a,1,i_z)])
                %stay_prob(i_a+n_A*(i_z-1)) = stay_prob(i_a+n_A*(i_z-1)) - (1-jfr_E)*OC_E1_VFI(i_a,1,i_z) ;
            else 
                % Move to employment 
                mTransition(n_State+1:end,i) = mTransition(n_State+1:end,i) + ...
                    jfr_E*mE_Transition_E(i_ep)*OC_E2_VFI(i_a,i_ep,i_z)*mTransition_aux(n_State+1:end,i_a+n_A*(i_z-1));
                mTransition(n_State+1:end,i_a+n_A*(i_z-1)) = mTransition(n_State+1:end,i_a+n_A*(i_z-1)) - ...
                    jfr_E*mE_Transition_E(i_ep)*OC_E2_VFI(i_a,i_ep,i_z)*mTransition_aux(n_State+1:end,i_a+n_A*(i_z-1));
                %disp([i_ep stay_prob(i_a+n_A*(i_z-1)) jfr_E*mE_Transition_E(i_ep)*OC_E2_VFI(i_a,i_ep,i_z)])
                %stay_prob(i_a+n_A*(i_z-1)) = stay_prob(i_a+n_A*(i_z-1)) - jfr_E*mE_Transition_E(i_ep)*OC_E2_VFI(i_a,i_ep,i_z) ;
            end 
            i = i+1 ;
        end 
        end 
        end 
        % Multiply columns by stay probability
        % mTransition(n_State+1:end,n_State+1:end) = repmat(max(stay_prob,0),[n_A*n_Z,1]).*mTransition(n_State+1:end,n_State+1:end);
        
    % Compute invariant histogram
    errHistogram = 100;	iterationHistogram = 0;
    G = ones(n_State+n_A*n_Z,1) ./ (n_State+n_A*n_Z);
    while errHistogram > 1e-12 && iterationHistogram < 1e4
        G_New = mTransition' * G;
        errHistogram = max(abs(G_New - G)); 
        iterationHistogram = iterationHistogram + 1;
        G = G_New;
    end

    % Expand histogram matrix
    G_W = G(1:n_State)     ; 
    G_E = G(n_State+1:end) ;
    mDBN_W = reshape(full(G_W),n_A,n_E,n_Z);    
    mDBN_E = reshape(full(G_E),n_A, 1 ,n_Z);   
        
        
%% Update prices 
    % Labor supply and new taxes
    K_demand = min( max( 0 , llambda*mA_Grid_E ) , (AA*mZ_Grid.*(aalpha/(r+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w))^mmu).^(1/(1-aalpha-mmu)) ) ;
    N_supply = sum( vE_Grid(2:n_E).*squeeze(sum(sum(mDBN_W(:,2:n_E,:),3),1))' ) ;
    N_demand = sum(sum(sum( mDBN_E .* (mmu*AA*mZ_Grid.*K_demand.^aalpha/((1+tau_n)*w)).^(1/(1-mmu)) ))) ;
    
    % tau_n = vE_Grid(1)*sum(sum(mDBN_W(:,1,:))) / N_supply ; 
        
%         % Update Wage
%         f = @(W) ( sum(sum(sum( mDBN_E .* ...
%             (mmu*AA*mZ_Grid.*(  min( max( 0 , llambda*mA_Grid ) , (AA*mZ_Grid.*(aalpha/(r+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*W))^mmu).^(1/(1-aalpha-mmu)) )  ).^aalpha/((1+tau_n)*W)).^(1/(1-mmu)) ...
%                     ))) /N_supply-1 ) ;
%         options = optimoptions('fsolve','Display','off');
%         [w_new,err,exitflag] = fsolve(f, w ,options);
%         
%         N_demand = sum(sum(sum( mDBN_E .* ...
%             (mmu*AA*mZ_Grid.*(  min( max( 0 , llambda*mA_Grid ) , (AA*mZ_Grid.*(aalpha/(r+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w_new))^mmu).^(1/(1-aalpha-mmu)) )  ).^aalpha/((1+tau_n)*w_new)).^(1/(1-mmu)) ...
%                     ))) ;
        
% disp([err N_supply N_demand sum(sum(mDBN_W(:,2,:))) sum(sum(mDBN_W(:,3,:))) w w_new])

%% Residual

        if strcmp(solver,'fsolve')
            residual = (N_supply/N_demand-1) ;
        elseif strcmp(solver,'fminsearch')
            residual = abs(N_supply-N_demand) ;
        elseif strcmp(solver,'fzero')
            residual = (N_supply-N_demand) ;
        elseif strcmp(solver,'bisection')
            residual = (N_supply-N_demand) ;
        else
            error('Invalid solver. It must be fsole, fminsearch, fzero or bisection')
        end
        % disp([w w_new])
    
%% Optional output
    if nargout>1
        mDBN_W_out = mDBN_W ;
        mDBN_E_out = mDBN_E ;
        if nargout>3
            mAp_W_out = Ap_W_VFI ;
            mAp_E_out = Ap_E_VFI ;
            OC_W_out  = OC_W_VFI ;
            OC_E_out  = OC_E_VFI ;
            VW_out    = VW_VFI   ;
            VE_out    = VE_VFI   ;

            mTransition_out = NaN(3,3) ; 
            U_ind = [ repmat(kron([1 ; 0 ; 0],ones(n_A,1)),n_Z,1) ; zeros(n_State,1) ]==1 ;
            W_ind = [ repmat(kron([0 ; 1 ; 1],ones(n_A,1)),n_Z,1) ; zeros(n_State,1) ]==1 ;
            E_ind = [zeros(n_State,1) ; ones(n_State,1)]==1 ;
            % From Unemployment
            mTransition_out(1,1) = sum(sum(mTransition(U_ind,U_ind),2).*G(U_ind))/sum(G(U_ind)) ;
            mTransition_out(1,2) = sum(sum(mTransition(U_ind,W_ind),2).*G(U_ind))/sum(G(U_ind)) ;
            mTransition_out(1,3) = sum(sum(mTransition(U_ind,E_ind),2).*G(U_ind))/sum(G(U_ind)) ;
            % From Employment
            mTransition_out(2,1) = sum(sum(mTransition(W_ind,U_ind),2).*G(W_ind))/sum(G(W_ind)) ;
            mTransition_out(2,2) = sum(sum(mTransition(W_ind,W_ind),2).*G(W_ind))/sum(G(W_ind)) ;
            mTransition_out(2,3) = sum(sum(mTransition(W_ind,E_ind),2).*G(W_ind))/sum(G(W_ind)) ;
            % From Unemployment
            mTransition_out(3,1) = sum(sum(mTransition(E_ind,U_ind),2).*G(E_ind))/sum(G(E_ind)) ;
            mTransition_out(3,2) = sum(sum(mTransition(E_ind,W_ind),2).*G(E_ind))/sum(G(E_ind)) ;
            mTransition_out(3,3) = sum(sum(mTransition(E_ind,E_ind),2).*G(E_ind))/sum(G(E_ind)) ;

            mTransition_out = 100*mTransition_out ;
        end 
    end 




end 