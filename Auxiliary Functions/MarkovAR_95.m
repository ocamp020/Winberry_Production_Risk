% =============================================================================
% MarkovAR_95: Approximates a contiuous AR(1) process with a discrete Markov process - Rowenhorst 95
%
% Usage: Call MarkovAR_95(n_z,rho,sigma,z_grid,P)
%
% Input: n_z    , integer , Grid dimension for Markov process
%        rho    , real(dp), Persistence of x where x'=rho*x+sigma*e and e~N(0,1)
%        sigma  , real(dp), Standard Deviation of shocks to x where x'=rho*x+sigma*e and e~N(0,1)
%
% Output: z_grid , real(dp), dimension(n_z),     Grid of values of the Markov process. Size n_z.
%         P      , real(dp), dimension(n_z,n_z), Transition Matrix of z. Prob sum across rows.
%
% Remarks: This routine generates a grid for z (z_grid) and a transition probability matrix (P)
%          The transition probability matrix is organized (z,z') so that it sums to 1 by rows
% 		   P(i,j) = p_ij is the transition prob of going from z=z_i to z'=z_j
%          Note that this creates matrixes in the opposite order than MarkovAR that follows Tauchen
%          The method follows Rouwenhorst (1995) as shown in Kopecky and Suen (2010)

function [z_grid,P_z] = MarkovAR_95(n_z,rho,sigma)

    z_grid = NaN(n_z,1)   ;
    P_z    = NaN(n_z,n_z) ;
    P_b    = NaN(n_z,n_z);
    P_half = NaN(n_z,n_z);

    % Parameters p, q and psi
        p = (1+rho)/2 ;
        q = (1+rho)/2 ;
        psi = sqrt(real(n_z)-1)*sigma/sqrt(1-rho^2) ;

    % Step of grid
        step = 2*psi/(n_z-1) ;

    % Compute z_grid
        z_grid(1) = -psi ;
        for i=2:n_z
            z_grid(i) = z_grid(i-1) + step ;
        end

    % Compute transition matrix for n_z=2
        P_2 = [p,1-p;1-q,q] ;

    % Compute transition matrix for arbitrary n_z 
        if n_z>2
            [~,P_aux] = MarkovAR_95(n_z-1,rho,sigma) ;
            % To create P_n you take P which is n_z -1 by n_z - 1 and create
            % 4 matrixes which have P in one corner and 0's in the row and column
            % touching the other corner.  For example,
            % [1 2; 3 4] => [1 2 0; 3 4 0; 0 0 0], [ 0 1 2; 0 3 4; 0 0 0] ...
            % plus
            P_a = zeros(n_z+1,n_z+1) ;
            P_a(2:n_z, 2:n_z) = P_aux ;
            P_b = ( p*P_a(2:n_z+1,2:n_z+1) + (1-p)*P_a(2:n_z+1,1:n_z) + ...
                    (1-q)*P_a(1:n_z,2:n_z+1) + q*P_a(1:n_z,1:n_z) ) ;
            P_half(1,:)       =  1.0  ;
            P_half(2:n_z-1,:) =  0.50 ;
            P_half(n_z,:)     =  1.0  ;
            P_z = P_b.*P_half ;
        else
            P_z = P_2 ;
        end 


end