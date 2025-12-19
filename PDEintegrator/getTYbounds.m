function [T,Y,extraouput] = getTYbounds(x,problem)
% computes the map T and the Y bounds
% The map T is used in the numerics, the bound Y is used in the proof
%
% The extraoutput 
% (used by generatecheb.m to determine number if Chebyshev nodes)
% contains
% extraouput.bounds1 : the residue from evaluating at interpolation nodes
%                  including bound on quadrature errors
%                  (or reather the corresponding interpolation polynomial)
% extraouput.bounds2 : a bound on the interpolation error

% We estimate the C0 norm by computing the C0 norm of a high degree
% interpolation polynomial and then adding the interpolation error
% estimate. For the higher Fourier modes, we also take the minimum 
% with the cruder bound obtained by taking the sup of the integrand (which 
% does not require interpolation error estimates, and might therefore be 
% better behaved for rather large values of n, for which we have to 
% interpolate very steep exponentials).

%%%%%%%%%%%%%%%%%%%
%  Initialization %
%%%%%%%%%%%%%%%%%%%

Ybounds=(nargout>1);

NQ=problem.semigroup.N;
NX=problem.sol.N;
K=problem.sol.K;
D=problem.sol.D;
symmetry=problem.symmetry;

fullintegrand = getintegrand(x,problem);

x0 = x;
for d = 1:D
    x0{d} = x0{d}(:,1);
end
problem0 = problem;
problem0.sol.K = zeros(size(problem.sol.K));
fullintegrand0 = getintegrand(x0,problem0);

% number of Fourier modes is the same for all domains
Nmax = (size(fullintegrand{1},1)-1)/2;
indQsym = symNindex(NQ,Nmax,symmetry);
indX = Nmax+1+(-NX:NX);

K0 = problem.proof.interpolation.K0;
K1 = problem.proof.integrals.K1;

T=cell(D,1);
nu = problem.proof.nu; 
Yinterp = altzeros([D,1],x{1}(1));
Yinterperror = altzeros([D,1],x{1}(1));
Ytailbound = altzeros([D,1],x{1}(1));
% thetau (already multiplied by Qinv)
thetau = altzeros([2*Nmax+1,1],x{1}(1));
thetau(indX) = problem.initial.data;

thetau=symQinvproject(thetau,problem.semigroup.Qinv(:,:,1),symmetry);

indL=symNindex(NQ,Nmax,symmetry);

integrationdata=[];

for d = 1:D
    tau=problem.domains(d)/2;
    Q=problem.semigroup.Q(:,:,d);
    Qinv=problem.semigroup.Qinv(:,:,d); 
    t1 = cheb_grid(K0(d));

    lambda = computelambdatail(Nmax,problem,d);
    % Putting in the correct lambdas for |n|<=N 
    lambda(indL) = problem.semigroup.Lambda(:,d); 

    % full nonlinearity 
    phi = fullintegrand{d}; 

    phi0 = fullintegrand0{d}; 
    phir = phi;
    phir(:,1) = phir(:,1) - phi0;
    LinvQinvphi0 = (1./lambda).*symQinvproject(phi0,Qinv,symmetry);
    Linvphi0 = symQproject(LinvQinvphi0,Q,symmetry);
    extraterm = (tau*(t1+1).*expm1div(tau*lambda*(t1+1))) .* symQinvproject(phi0,Qinv,symmetry);
    
    phi = symQinvproject(phi,Qinv,symmetry);
    phir = symQinvproject(phir,Qinv,symmetry);

    % We start by T 
    expinterp_val = exp(tau*lambda*(t1+1)).*thetau;

    [integralinterp_val,integrationdata] = computeintegral(phir,tau*lambda,K0(d),K1(d),problem,integrationdata);
    Tval = expinterp_val + extraterm + tau*integralinterp_val;
    Tval=symQproject(Tval,Q,symmetry);
    Tinterp = cheb_coeffs(Tval);

    % the map T for the numerics
    T{d}=Tinterp(indX,1:K(d)+1);

    if problem.infinitetime && d==D
        % fixed point problem for equilibrium
        Tval=-symQproject((1./lambda).*phi(:,1),Q,symmetry);
        T{D}(:,1)=Tval(indX);
        T{D}(:,2:end)=0;
    end

    if Ybounds
        if ~isfield(problem,'nodisplay') || ~problem.nodisplay
            disp(['Computing the Y-bound on domain ',num2str(d)])
        end
        
        %%% We estimate, on each subdomain, T(\bu) - \bu by splitting it into
        %%% P_{\tilde K}(T(\bu) - \bu) + (I-P_{\tilde K})T(\bu). We have two different estimates
        %%% for the first term, and take the minimum among those. The
        %%% approximate solution \bu is called x in the code, and \tilde K is K0.

        Tminusuinterp = Tinterp; %P_{\tilde K}T(\bu)
        Tminusuinterp(indX,1:K(d)+1) = Tinterp(indX,1:K(d)+1) - x{d}; %P_{\tilde K}T(\bu) - \bu
        
        % The C0 norm of the interpolation of T-I on domain d, bounded
        % using the ell1 norm of.
        TminusuinterpC0 = normC0(Tminusuinterp);
        
        % Alternative one can bound the C0 norm of the interpolation polynomial
        % by the maximum of the interpolation values times the Lebesgue constant
        if not(problem.infinitetime && d==D)
            Tminusuval = Tval;
            Tminusuval(indX,:) = Tminusuval(indX,:) - cheb_eval(x{d}, t1);
            TminusuinterpC0bis = lebesgueconstant(K0(d)) * max(abs(Tminusuval),[],2);                
            TminusuinterpC0 = min(TminusuinterpC0, TminusuinterpC0bis);
        end
        
        % The interpolation error estimate for (I-P_{\tilde K})T(\bu).
        % T(\bu) is the sum of an exponential term and an intergral term,
        % and we compute an interpolation error bound for each term
        % separately.
        Nlambda=length(lambda);
        interperrorexp = altzeros([Nlambda,1],x{1}(1));
        interperrorint = altzeros([Nlambda,1],x{1}(1));
        prefactor = @(rho) 4*rho.^(-K0(d))./(rho-1);
        prefactor_int = @(rho) 2*rho.^(-K0(d)).*(rho+rho^(-1)+2)/(rho-1);

        for n = 1:Nlambda
            % analytic interpolation estimate for the exponential
            interperrorexp(n) = analyticerrorestimate(prefactor,tau*lambda(n),1,K0(d));
            % analytic interpolation estimate for the integral
            interperrorint(n) = tau*analyticerrorestimate(prefactor_int,...
                                        tau*lambda(n),phir(n,:),K0(d),[],@expm1div);
        end   
        % the exponential term
        expterm = abs(thetau+LinvQinvphi0);
        expterm(lambda==0) = 0;
        interperrorexp = interperrorexp.*expterm; 
        
        % add the two terms
        interperror = interperrorexp+interperrorint;
        interperror= abs(symQproject(interperror,abs(Q),symmetry));
        
        % the crude error bound obtained by taking separate C^0 norms
        % (once could do this first, and not compute the other estimates if
        % this is already below some threshold, this might speed up the
        % whole proof a bit)
        crude1 = exp(2*tau*max(real(lambda),0)) .* expterm ;
        crude2 = 2*tau.*expm1div(2*tau*real(lambda)) .* normC0(phir);
        xdr = altzeros([length(Linvphi0),size(x{d},2)],x{d}(1));
        xdr(indX,:) = x{d};
        xdr(:,1) = xdr(:,1) + Linvphi0;
        crude3 = normC0(xdr);

        crude = abs(symQproject(crude1+crude2,abs(Q),symmetry));
        crude(:) = crude(:) + crude3;
        
        % Checking whether the crude error bound might be better
        crudeisbetter = crude < altsup(TminusuinterpC0+interperror);
        TminusuinterpC0(crudeisbetter) = crude(crudeisbetter);
        interperror(crudeisbetter) = 0;
        
        % Final norm for the interpolation of T-I on domain d
        Yinterp(d) = normL1(TminusuinterpC0,nu);
        
        % Final interpolation error bound for T-I on domain d
        Yinterperror(d) = normL1(interperror,nu);

        if problem.infinitetime && d==D
            % overwrite Yinterperror(D) and Yinterp(D)
            % phi is Qinvgamma
            if all(real(lambda)<0)
                LQinvgamma=(1./lambda).*phi(:,1);
                constantterm=-symQproject(LQinvgamma,Q,symmetry);
                constantterm(indX)=constantterm(indX)-x{D}(:,1);

                timedependent=LQinvgamma+thetau;
                timedependent=symQproject(abs(timedependent),abs(Q),symmetry);

                Yinterp(D)=normL1(constantterm)+normL1(timedependent);
            else
                Yinterp(D)=Inf;
            end
            Yinterperror(D)=0;
        end
        mu = boundmu(d,0,NX,problem);
        Ytailbound(d) = exp(mu)*problem.initial.tailbound;
    end %% of Ybounds part

    % Updating thetau for the next subinterval
    if d < D
        % recurrence relation
        thetau = expinterp_val(:,1) + extraterm(:,1) + tau*integralinterp_val(:,1); 
        % basis transformation
        thetau(indQsym) = (problem.semigroup.Qinv(:,:,d+1) * Q) * thetau(indQsym);
        % It would be more efficient cost-wise to do Qinv*(Q*thetau)
        % instead of (Qinv*Q)*thetau, but the latter incurs less wrapping
    
    end
end
    
if Ybounds   
    Y = Yinterp + Yinterperror + Ytailbound;
    extraouput.bounds1=Yinterp;
    extraouput.bounds2=Yinterperror;
    if problem.infinitetime
        extraouput.bounds1(end)=[];
        extraouput.bounds2(end)=[];
        Ytest = max(altsup(Yinterperror(1:end-1)))/max(altsup(Y(1:end-1)));
        if ~isfield(problem,'nodisplay') || ~problem.nodisplay
            disp(['size of Y on finite domains is ',num2str(altsup(max(Y(1:D-1))))]);
            disp(['size of Y on final infinite domain is ',num2str(altsup(Y(D)))]);
        end
    else
        Ytest = max(altsup(Yinterperror))/max(altsup(Y));
        if ~isfield(problem,'nodisplay') || ~problem.nodisplay
            disp(['size of Y is ',num2str(altsup(max(Y)))]);
        end
    end
    if ~isfield(problem,'nodisplay') || ~problem.nodisplay            
        if Ytest > 10^-1
            fprintf("You may want to increase K0 to reduce the size of the residue.\n")
            fprintf("There is a non-negligible part (around %g percent) of the Y bound coming from interpolation errors.\n",100*Ytest)
        elseif Ytest < 10^-3
            fprintf("You could try decreasing K0 to reduce the cost of the proof.\n")
            fprintf("The ratio between interpolation error and C^0 error is %g.\n",Ytest)
        end
    end
end

if ~Ybounds
    Y=[];
    extraouput=[];
end

