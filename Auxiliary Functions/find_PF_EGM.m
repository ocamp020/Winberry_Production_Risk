function [Ap,C] = find_PF_EGM(r,p,w)

    % Tolerance
    Max_Iter = 10000; 
    Tol_EGM  = 1e-6 ;

    % Initialization (global variables)
    global bbeta ggamma ssigma aaBar ddelta mmu llambda tau_k ...
    	   mEpsilonTransition   ...
           nEpsilon nZ nAssetsFine  ...
           vAssetsGridFine  mEpsilonGridFine mAssetsGridFine mZGridFine


    %% Firm's capital demand and profits
        K       = min( llambda*mAssetsGridFine , (mmu*p*mZGridFine.^mmu/(r+ddelta)).^(1/(1-mmu)) )  ;
        Pr      = p*(mZGridFine.*K).^mmu - (r+ddelta)*K ;
        Mrg_Pr  = (mmu*p*(mZGridFine*llambda).^mmu.*mAssetsGridFine.^(mmu-1)-(r+ddelta)*llambda) .* (K==llambda*mAssetsGridFine) ;
        Y       =  Pr*(1-tau_k) + (1+r*(1-tau_k))*mAssetsGridFine ;
        
    %% Initial Values 
        C      = w*mEpsilonGridFine + (Pr + r*mAssetsGridFine)*(1-tau_k) ;
        ECp    = NaN(nAssetsFine,nEpsilon,nZ) ;
            for i_z=1:nZ
               ECp(:,:,i_z) = C(:,:,i_z).^(-ssigma)*mEpsilonTransition' ; 
            end
        EVp    = (1+(r+Mrg_Pr)*(1-tau_k))*bbeta*ggamma.*ECp ;
        C_exo  = C ;
 
    %% EGM 
    Dist = 1 ;
    iter   = 1 ;
    while iter<=Max_Iter && Dist>Tol_EGM     
        % Solve for optimal savings with EGM
            C_endo = EVp.^(-1/ssigma)  ;
            Y_endo = C_endo + mAssetsGridFine - w*mEpsilonGridFine ;
            
        % Get optimal savings on exogenous grid
            % Interpolate to get exogenous consumption
            for i_z = 1:nZ
            for i_e = 1:nEpsilon 
                % Use interpolation for points such that Y_mat >= min(Y_endo)
                % ind = (A_grid>=min(A_endo(:,i_e))) & (A_grid<=max(A_endo(:,i_e))) ;
                ind = (Y(:,i_e,i_z)>=min(Y_endo(:,i_e,i_z))) ;
                C_exo(ind,i_e,i_z)  = spline(Y_endo(:,i_e,i_z),C_endo(:,i_e,i_z),vAssetsGridFine(ind)) ;
                % If Y_mat<min(Y_endo) then solve the FOC non-linearly
                if sum(~ind)>0
                    pp = spline(Y(:,i_e,i_z),EVp(:,i_e,i_z)) ;
                    % A_aux = linspace(1,max(A_grid)+2,10000);
                    % EVp_spline = ppval(pp,A_aux);
                    aux = 1:nAssetsFine ;
                    for i_a = aux(~ind)
                        f = @(ap) Euler_Residual(ap,pp,w,i_a,i_e,i_z,Y(i_a,i_e,i_z));
                        options = optimoptions('fsolve','Display','off');
                        [Ap_foc,err,exitflag] = fsolve(f, vAssetsGridFine(i_a) ,options);
                        C_exo(i_a,i_e,i_z) =  w*mEpsilonGridFine(i_a,i_e,i_z) + Y(i_a,i_e,i_z) - Ap_foc ;
                        if err>1e-5 
                            disp('Error in fsolve')
%                             disp([A_grid(i_a) Ap_foc exitflag min(A_endo(:,i_e)) max(A_endo(:,i_e))])
%                             disp(err)
%                             figure; plot(A_grid,EVp(:,i_e),A_aux,EVp_spline,':'); legend('real','spline')
%                             figure; plot(A_grid,C(:,i_e)); legend('real','spline')
%                             pause(5)
                            return 
                        end 
                    end 
                end 
            end
            end 
            
            % Check for bounds and correct consumption
            Ap    = max( w*mEpsilonGridFine + Y - C_exo , aaBar) ;
            C_exo =      w*mEpsilonGridFine + Y - Ap ;
            
            
            
            
        % Updates and distance
            Dist = max(abs(C(:)-C_exo(:)))    ;
            ECp    = NaN(nAssetsFine,nEpsilon,nZ) ;
                for i_z=1:nZ
                   ECp(:,:,i_z) = C_exo(:,:,i_z).^(-ssigma)*mEpsilonTransition' ; 
                end
            EVp    = (1+(r+Mrg_Pr)*(1-tau_k))*bbeta*ggamma.*ECp ;
            C      = C_exo                      ;
            iter   = iter + 1                   ;
    end 
    
    Ap = w*mEpsilonGridFine + Y - C ;


end 


function res = Euler_Residual(ap,pp,w,i_a,i_e,i_z,Y)
    global ssigma mEpsilonGridFine

    EVp_spline = ppval(pp,ap) ; 
    cons       = w*mEpsilonGridFine(i_a,i_e,i_z) + Y - ap ;
    res = ( (cons)^(-ssigma) - EVp_spline );
    
end