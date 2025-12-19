function Dphifull = getDphifull(x,problem)
% determines Dphifull, which contains all nonzero terms of the derivatives 
% of the nonlinearities in the integrand without the exponential 
% and without the semigroup residue
% (Fourier-Chebyshev coefficient format)

K=problem.sol.K;
D=problem.sol.D;

Dphifull=cell(length(problem.pde.polynomials),D);
  

for d=1:D
    
    % extend x so that we can use convolutions to compute nonlinearities
    xx=[x{d}(:,K(d)+1:-1:2),x{d}];

    % nonlinear terms 
    for p=0:length(problem.pde.polynomials)-1
        if ~isempty(problem.pde.polynomials{p+1})
            % terms with p-th derivative
            [~,Dphip]=convnonlinearity(xx,p,1,problem);
            Dphifull{p+1,d}=Dphip;
        end
    end

end % end of loop over domain

end









