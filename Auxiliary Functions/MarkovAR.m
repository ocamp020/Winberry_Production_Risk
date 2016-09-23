% =============================================================================
% MarkovAR: Approximates a contiuous AR(1) process with a discrete Markov process - Tauchen 86
%
% Usage: Call MarkovAR(n_z,lambda,rho,sigma,z_grid,P)
%
% Input: n_z    , integer , Grid dimension for Markov process
%        lambda , real(dp), Expansion of the grid: [-lambda*std(x) , lambda*std(x)] where x is AR(1)
%        rho    , real(dp), Persistence of x where x'=rho*x+sigma*e and e~N(0,1)
%        sigma  , real(dp), Standard Deviation of shocks to x where x'=rho*x+sigma*e and e~N(0,1)
%
% Output: z_grid , real(dp), dimension(n_z),     Grid of values of the Markov process. Size n_z.
%         P      , real(dp), dimension(n_z,n_z), Transition Matrix of z. Prob sum across cols.
%
% Remarks: This routine generates a grid for z (z_grid) and a transition probability matrix (P)
%          The transition probability matrix is organized (z',z) so that it sums to 1 by cols
% 		   P(j,i) = p_ij is the transition prob of going from z=z_i to z'=z_j
%          This method is taken from Tauchen (1986)

function [z_grid,P] = MarkovAR(n_z,lambda,rho,sigma)
	    
    z_grid = NaN(n_z,1)   ; 
    P      = NaN(n_z,n_z) ;
    
    % Step of grid
    step = 2*(lambda*sigma/sqrt(1-rho^2))/(n_z-1) ;

    % Compute z_grid
    z_grid(1) = -1*(lambda*sigma/sqrt(1-rho^2)) ;
    for i=2:n_z
        z_grid(i) = z_grid(i-1) + step ;
    end

    % Compute transition matrix    
    for i=1:n_z
        P(1,i) = normcdf( (z_grid(1)-rho*z_grid(i)+step/2)/sigma ) ;
        for j=2:n_z-1
            P(j,i) = normcdf( (z_grid(j)-rho*z_grid(i)+step/2)/sigma ) - normcdf( (z_grid(j)-rho*z_grid(i)-step/2)/sigma ) ;
        end
        P(n_z,i)=1-sum(P(1:n_z-1,i)) ;
    end

end