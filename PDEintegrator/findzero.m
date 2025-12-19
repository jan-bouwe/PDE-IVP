function [xout,residue] = findzero(xin,problem,updatesemigroup)
% performs iterations of the map T 
% to find a zero of the problem F=0 
% starting at xin

x=xin;

problem.type='numerics';

if nargin<=2
    updatesemigroup=false;
end
    

% initialize iterations
nsteps=24; 
tolerance=1e-13;
normdx=ones(nsteps+1,1);
n=1;

% delay is used to check it does not aimlessly wander
% with small numbers larger than tolerance but smaller than relaxedtol
delay=2;
relaxedtol=1e-11;
% for default behaviour choose delay=0; relaxedtol=tolerance;
% delay=0;
% relaxedtol=tolerance;

if ~isfield(problem,'nodisplay') || ~problem.nodisplay
    disp('Start of iterations of (contracting) map')
end

while n<=nsteps && normdx(n)>=tolerance && ~(n>delay && all(normdx(n-delay:n)<relaxedtol))
    xnew=getTYbounds(x,problem);
    % enforce symmetry
    for d=1:length(xnew)
        xnew{d}=(xnew{d}+conj(xnew{d}(end:-1:1,:)))/2;
        if isfield(problem,'symmetry') 
            if strcmp(problem.symmetry,'cosineseries')
                xnew{d}=real(xnew{d});
            elseif strcmp(problem.symmetry,'sineseries')
                xnew{d}=1i*imag(xnew{d});
            end
        end
    end
    n=n+1;
    normdx(n)=distancex(xnew,x,problem.proof.nu);
    if ~isfield(problem,'nodisplay') || ~problem.nodisplay
        disp(['stepsize is ',num2str(normdx(n))]);
    end
    x=xnew;
    if updatesemigroup
        problem.semigroup = fixsemigroup(x,problem);
    end
end

xout=x;
residue=distancex(getTYbounds(x,problem),x,problem.proof.nu);

if ~isfield(problem,'nodisplay') || ~problem.nodisplay
    disp('End of iterations')
end

end

function nx=distancex(x1,x2,nu)
% norm of x
D=length(x1);
nx=altzeros([D,1],x1{1}(1));
for d=1:D
    nx(d)=normL1(normC0(x1{d}-x2{d}),nu);
end
nx=max(nx);
end



