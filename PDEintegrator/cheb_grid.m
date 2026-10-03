function [t,theta] = cheb_grid(K,intvaltest)
% returns the Chebyshev grid with K+1 grid points
% if intvaltest is true t is interval-valued
% the optional second output theta contains the angles of the grid points,
% so that t = cos(theta)

if exist('intvaltest','var') && altisintval(intvaltest)
   ipi = intval('pi'); 
else
   ipi = pi; 
end

theta = (0:K)*ipi/K;
t = cos(theta);
t(end) = -1;
    
end

