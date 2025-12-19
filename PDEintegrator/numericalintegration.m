function [x,alldata]=numericalintegration(problem)
% Numerical integration of the Fourier-truncated PDE
% to find a guess for the solution
% The number of Chebyshev nodes is K+1
% The initial data are stored in problem.initial.data

K=problem.sol.K;
D=problem.sol.D;

domains=diff(problem.timegrid);
x=cell(D,1);

% no intervals; this is numerics
for p=1:length(problem.pde.polynomials)
    problem.pde.polynomials{p}=altmid(problem.pde.polynomials{p});
end

if problem.infinitetime
    D=D-1;
end
initial=problem.initial.data;

alltimes=zeros(1+sum(K),1);
counter=1;
% generate all nodes
for d=1:D
    % in each domain we want the solution at the Chebyshev nodes
     times=(fliplr(cheb_grid(K(d)))+1)/2*domains(d);
     alltimes(counter+(1:K(d)))=problem.timegrid(d)+times(2:end)';
     counter=counter+K(d);
end

% numerical integration
[~,allorbits]=ode15s(@(t,a) vectorfieldpde(a,problem),alltimes,initial);
% [~,allorbits]=ode45(@(t,a) vectorfield(a,problem),times,initial);
% opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
% [~,allorbits]=ode15s(@(t,a) vectorfield(a,problem),times,initial,opts);

if D==1 && K(1)==1
     % extracting begin and end point
     allorbits=allorbits([1 end],:);
end
alldata.times=alltimes;
alldata.orbit=allorbits;

% store the Chebyshev coefficients
counter=1;
for d=1:D
 %  orbit=allorbits(1+(d-1)*K:1+d*K,:);
    orbit=allorbits(counter+(0:K(d)),:);
    x{d}=cheb_coeffs(fliplr(orbit.'));
    counter=counter+K(d);
end

if problem.infinitetime
    x{D+1}(:,1)=allorbits(end,:).';
end

end

function f=vectorfieldpde(a,problem)
% vectorfield of the ODE which results from truncating the PDE
    N=(length(a)-1)/2;
    % leading term
    f=(-((-N:N)').^problem.pde.order).*a;
    for p=0:length(problem.pde.polynomials)-1
        % terms with p-th derivative
        cp=convnonlinearity(a,p,0,problem);
        cp=diag((1i.*(-N:N)).^p)*cp;
        f=f+cp;
    end
end
