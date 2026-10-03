function [ y ] = cheb_eval(a,t,theta)
% evaluates the Chebyshev polynomial 
% corresponding to the Chebyshev coefficients a
% the optional input theta contains the angles of the points t,
% so that t = cos(theta) and acos is not needed

K=size(a,2)-1;
k=(0:K)';
if exist('theta','var')
    cosktheta=cos(k*theta);
else
    cosktheta=cos(k*acos(t));
end

a(:,2:end)=2*a(:,2:end);
% absorb the factor 2 in the definition of the Chebyshev coefficients 
y=a*cosktheta;

end
