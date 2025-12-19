function Lambda = lebesgueconstant(K,intvaltest)
% Computes the Lebesque constant (with interval arithmetic if intvaltest is an interval) 
% When K is odd, we use an exact formula. When K is even, we use an
% upper-bound, except for small values (namely K = 2 and K = 4), for which
% we computed the exact value

if exist('intvaltest','var') && exist('intval.m','file') && isintval(intvaltest(1))
    ipi = intval('pi');
    sqrt2 = sqrt(intval(2));
    % computation needed for the case K = 4
    x = infsup(inf(intval('0.37227692859622')),sup(intval('0.37227692859624')));
%     x = verifynlss('8*x^3-6*(1+sqrt(0*x+2))*x^2-6*x+1+2*sqrt(0*x+2)',0.38);    
else
    ipi = pi;
    sqrt2 = sqrt(2);
    % computation needed for the case K = 4
    x = 0.37227692859623;
%     x = fzero('8*x^3-6*(1+sqrt(2))*x^2-6*x+1+2*sqrt(2)',0.38);
end

if K == 2
    Lambda = 0*x+5/4;
elseif K == 4
    Lambda = 2*x^4-2*(1+sqrt2)*x^3-3*x^2+(1+2*sqrt2)*x+1;
else
    Lambda = sum(cot((1:2:2*K-1)*ipi/(4*K)))/K; %sharp if K odd, possible overestimation by at most 1/K^2 if K even
    % Lambda = 1+2/ipi*log(K+1); % valid upperbound, but less sharp than the one above
end
    
end 