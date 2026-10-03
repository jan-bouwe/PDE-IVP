function lambda = computelambdatail(N,problem,d)

% Computes lambda_n, for n = -N,...,N or 0:N or 1:N (depending on symmetry), 
% given by
% lambda_n = -n^{2R} + \sum_{j=0}^{2R-1} (in)^j v^{(j)}_0.
% Note that these are NOT the correct eigenvalues of the semigroup for |n|
% <= NQ.
%
% One can either get the eigenvalue for a given subdomain, selected by the
% third input variable d, or the eigenvalues on all subdomains, by not
% inputting d (or by using d='all').
%
% The output is of intval type if problem.semigroup.Lambda consists of intervals

NQ=problem.semigroup.N;
iv0=(size(problem.semigroup.ToepV,1)+1)/2;
D = problem.sol.D;
der = 1i*symNrange(N,problem.symmetry)';
if isa(problem.semigroup.Lambda(1),'intval')
    der=intval(der);
end

if nargin < 3 || strcmp(d,'all') % eigenvalues on all subdomains
    lambda = altzeros([length(der),D],problem.semigroup.Lambda(1));
    lambdadom = -abs(der).^problem.pde.order;
    for d = 1:D
        lambdad = lambdadom;
        for p = 0:length(problem.pde.polynomials)-1
            if ~isempty(problem.pde.polynomials{p+1})
                lambdad = lambdad + der.^(p)*problem.semigroup.ToepV(iv0,d,p+1);
            end
        end
        lambda(:,d) = lambdad;
    end
else % eigenvalues on a given subdomain
    lambda = -abs(der).^problem.pde.order;
    for p = 0:length(problem.pde.polynomials)-1
        if ~isempty(problem.pde.polynomials{p+1})
            lambda = lambda + der.^(p)*problem.semigroup.ToepV(iv0,d,p+1);
        end
    end
end