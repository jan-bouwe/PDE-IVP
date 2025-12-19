function [fullintegrand,integrand] = getintegrand(x,problem)
% determines nonlinearity in the integrand without the exponential 
% (Fourier-Chebyshev coefficient format)
% 
% fullintegrand : all nonzero terms   
% integrand : truncated in the Fourier direction to match dimension of x

NQ=problem.semigroup.N;
NX=problem.sol.N;
indQ=NX+1+(-NQ:NQ);
K=problem.sol.K;
D=problem.sol.D;

% the integrand (before multiplying by the exponential)
% will be fully stored (all nonzero coefficients)
maxorder=max([1,nonlinmaxdegree(problem)]);
Nmax=maxorder*NX;
ToepV=problem.semigroup.ToepV;
fullintegrand=cell(D,1);

reducedintegrandrequested=(nargout>1);
if reducedintegrandrequested
    integrand=cell(D,1);
end

for d=1:D

    Kmax=maxorder*K(d);
    % extend x so that we can use convolutions to compute nonlinearities
    xx=[x{d}(:,K(d)+1:-1:2),x{d}];
                
    integranddomaind=altzeros([2*Nmax+1,Kmax+1],x(1));

    % nonlinear terms
    for p=0:length(problem.pde.polynomials)-1
        % terms with p-th derivative
        if ~isempty(problem.pde.polynomials{p+1})
            [~,phip]=convnonlinearity(xx,p,0,problem);
            phip=settensorsize(phip,[Nmax,Kmax]);
            phip=phip(:,Kmax+1:2*Kmax+1);
            phip=diag((1i.*(-Nmax:Nmax)).^p)*phip;
            integranddomaind=integranddomaind+phip;
        end
    end

    % subtract the linear semigroup term
    Tv=altzeros([2*NQ+1,2*NQ+1],xx(1));
    v0=altzeros([1,1],xx(1));
    DNQ=(1i*(-NQ:NQ));
    DNX=(1i*(-NX:NX));
    for p=0:length(problem.pde.polynomials)-1
        if ~isempty(problem.pde.polynomials{p+1})
            Tv=Tv+diag(DNQ.^p)*toeplitzshift(ToepV(:,d,p+1)); 
            v0=v0+DNX.^p*ToepV(NQ+1,d,p+1);
        end
    end
    linearpart=diag(v0)*xx;
    linearpart(indQ,:)=Tv*xx(indQ,:);
    linearpart=settensorsize(linearpart,[Nmax,Kmax]);
    % linear part has minus sign (it is subtracted since part of the semigroup)
    integranddomaind=integranddomaind - linearpart(:,Kmax+1:2*Kmax+1);

    % enforce symmetry
    if isfield(problem,'symmetry')
        if strcmp(problem.symmetry,'cosineseries')
            % cosine series: both x and integrand have real coefficients
            integranddomaind=real(integranddomaind);
        elseif strcmp(problem.symmetry,'sineseries')
            % cosine series: both x and integrand are purely imaginary
            integranddomaind=1i*imag(integranddomaind);
        end
    end

    % deal with the semigroup residue
    semigroupresidue = problem.semigroup.Residue(:,:,d)*xx(indQ,:);
    semigroupresidue = settensorsize(semigroupresidue,[Nmax,Kmax]);
    integranddomaind = integranddomaind - semigroupresidue(:,Kmax+1:2*Kmax+1);

    fullintegrand{d}=integranddomaind;
    if reducedintegrandrequested
        integrand{d}=integranddomaind(Nmax+1+(-NX:NX),:);
    end
    
end



end









