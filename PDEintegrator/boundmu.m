function mu = boundmu(d,l,N,problem)
% computes the bound mu needed to estimate 
% the inverse of the preconditioner for domain decomposition

% Warning: the definition changes slightly in the new paper.
% This code does not include the exponential.

if l == d
    mu = 0;
    return
elseif l > d
    error('Invalid inputs d and l')
end

tau=problem.domains/2;
P=problem.pde.order;
ToepV=problem.semigroup.ToepV;
NQ=problem.semigroup.N;

% only the even derivatives play a role
% since the odd ones lead to purely imaginary contributions to lambda
alpha=reshape(real(ToepV(NQ+1,1:d,1:2:P-1)),d,[]);
Nalpha=ceil(altsup(sqrt(sum(abs(alpha),2))));
Nmax=max(Nalpha,N);

% for n>Nmax the lambdas are definitely negative and decreasing
% (by simply estimating the leading term versus the sum of the others;
% this improves for the derivative)
Nmax=max(Nmax(l+1:d));

% compute the relevant eigenvalues
nvalues=N+1:Nmax+1;
rlambda=repmat(-nvalues'.^P,1,d);
if altisintval(ToepV(1))
    rlambda=intval(rlambda);
end
p=0:(P/2)-1;
n2psign=(nvalues').^(2*p).*(-1).^p;
rlambda=rlambda+n2psign*alpha';

% add all the terms 
if problem.infinitetime && d==length(tau)
    tauRlambda=0;
else
    tauRlambda=2*tau(d)*max(rlambda(:,d),0);
end
mtau = (l+1:d-1);
tauRlambda =  tauRlambda + 2*rlambda(:,mtau)*tau(mtau)';

% find maximum
mu=max(tauRlambda);

end

