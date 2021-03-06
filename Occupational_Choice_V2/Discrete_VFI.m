% Production Risk Economy with Occupational Choice
% Find Stationary Distribution with the Histogram Method
% Juan David Herreno & Sergio Ocampo (2016)

function [mAp_W_out,mAp_E_out,OC_W_out,OC_E_out,VW_out,VE_out] = Discrete_VFI(r,w)

%% Initialization (global variables)
    global bbeta ggamma ddelta aalpha mmu AA llambda tau_k tau_n vKappa ...
    	   mE_Transition_W_VFI mE_Transition_E_VFI  ...
           n_E n_Z n_A vA_Grid mA_Grid mZ_Grid mE_Grid...
           maxIterations tolerance 

%% Build Potential Consumption (size n_A x n_E x n_Z x N_Ap)

    % States
    mA_aux  = repmat(mA_Grid,[1, 1, 1,n_A]) ;
    mE_aux  = repmat(mE_Grid,[1, 1, 1,n_A]) ;
    mZ_aux  = repmat(mZ_Grid,[1, 1, 1,n_A]) ;
    
    % Future Assets 
    mAp_aux = permute(mA_aux,[4,2,3,1])     ;
      
    % Firm's capital demand and profits
    K       = min( max( 0 , llambda*mA_aux ) , ...
               (AA*mZ_aux.*(aalpha/(r+ddelta)).^(1-mmu).*(mmu/((1+tau_n)*w))^mmu).^(1/(1-aalpha-mmu)) )  ;
    N       = (mmu*AA*mZ_aux.*K.^aalpha/((1+tau_n)*w)).^(1/(1-mmu)) ;
    Pr      = AA*mZ_aux.*K.^aalpha.*N.^mmu - (r+ddelta)*K - (1+tau_n)*w*N;

    % Consumption
    mC_W =  (1+(1-tau_k)*r)*mA_aux + w*mE_aux     - mAp_aux ;
    mC_E =  (1+(1-tau_k)*r)*mA_aux + (1-tau_k)*Pr - mAp_aux ;
    

%%  Value Functions - Initial Conditions (size n_A x n_E x n_Z)
    mVW = Utility((1-tau_k)*r*mA_Grid + w*mE_Grid)/(1-bbeta*ggamma) ; 
    mVE = Utility((1-tau_k)*r*mA_Grid + (1-tau_k)*Pr(:,:,:,1))/(1-bbeta*ggamma) ; 
    
    % Allocate variables
    mVE_interp = NaN(size(mVW)) ;
    
%% VFI 
    iter   = 0 ;
    V_Dist = 1 ;
    while V_Dist>=tolerance && iter<=maxIterations
        
        % Get Interpolation for mVE with cost
        for i_e = 1:n_E 
        for i_z = 1:n_Z
            mVE_interp(:,i_e,i_z) = interp1(vA_Grid,mVE(:,i_e,i_z),vA_Grid-vKappa(i_z)) ;
            % Correct for agents that cannot afford the switching cost
            mVE_interp((vA_Grid-vKappa(i_z))<0,i_e,i_z) = -Inf ;
        end
        end
        
        % Get max for each future epsilon
        mVW_max = max(mVW,mVE_interp) ;
        mVE_max = max(mVW,mVE       ) ;
        
        % Replicate value to get new dimensions Ap and Ep 
        mVW_p   = permute( repmat( mVW_max , [1,1,1,n_A,n_E] ) , [4,5,3,1,2] ) ; % dimension: A,E,Z,Ap,Ep
        mVE_p   = permute( repmat( mVE_max , [1,1,1,n_A,n_E] ) , [4,5,3,1,2] ) ; % dimension: A,E,Z,Ap,Ep
        
        % Expected value over Ep
        E_V_W   = sum( mE_Transition_W_VFI.*mVW_p , 5 );
        E_V_E   = sum( mE_Transition_E_VFI.*mVE_p , 5 );
        
        % Maximize over Ap
        mVW_new = max(  Utility(mC_W) + bbeta*ggamma*E_V_W , [] , 4 );
        mVE_new = max(  Utility(mC_E) + bbeta*ggamma*E_V_E , [] , 4 );
        
        % Get distance
        V_Dist  = max( max( abs( mVW_new(:)./mVW(:) - 1  ) ) , max( abs( mVE_new(:)./mVE(:) - 1  ) ) ) ;
        
        % Update Value and Iteration
        mVW     = mVW_new ;
        mVE     = mVE_new ;
        iter    = iter+1  ;

    end

%% Policy Functions
    [~,mAp_W_out] = max(  Utility(mC_W) + bbeta*ggamma*E_V_W , [] , 4 ) ;
    [~,mAp_E_out] = max(  Utility(mC_E) + bbeta*ggamma*E_V_E , [] , 4 ) ;
    OC_W_out      = mVW<mVE_interp ; 
    OC_E_out      = mVW>mVE        ; 


%% Output 
    VE_out   = mVE ;
    VW_out   = mVW ;
    

end 


function U = Utility(C)
    global ssigma 
    if ssigma ~= 1
        U =   C.^(1-ssigma)/(1-ssigma) ;
    else
        U =   log(C) ;
    end
    U(C<0) = -Inf ;

end 