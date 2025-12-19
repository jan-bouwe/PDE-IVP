function problem = preprocessdata(problem,initialdata,K,NQ,griddata,symmetry,scale)
% preprocessing of the input data to define
% problem.initial.data
% problem.initial.tailbound
% problem.timegrid
% problem.timegridfloats
% problem.symmetry
% problem.scale
% problem.scalefloats
% problem.sol.N
% problem.sol.K
% problem.sol.D
% problem.semigroup.N
%
% keeps track of rounding errors when converting to intervals
%
% if griddata.variable is true it will try to generate a good grid
% if griddata.optimizecheb is true it will try to guess good number of Chebyshec modes

if all(imag(initialdata)==0)
    initialdatareal=true;
else
    initialdatareal=false;
end
if all(real(initialdata)==0)
    initialdataimag=true;
else
    initialdataimag=false;
end

% turn into intervals if appropriate
if exist('intval.m','file') && ~isintval(initialdata(1))
    % factor sqrt(2) for complex-valued data
    if initialdatareal || initialdataimag
        initialdataradii=eps(initialdata);
    else
        initialdataradii=sup(intval('sqrt2')*eps(initialdata));
    end
    zeroinitialdata=(initialdata==0);
    initialdataradii(zeroinitialdata)=0;
    initialdata=midrad(initialdata,initialdataradii);
end

% symmetrize to make initial data real; note the factor 1/2
initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

% store initial data
problem.initial.data=initialdata;
problem.initial.tailbound=0; % in the nu-norm
problem.sol.N=(size(initialdata,1)-1)/2;
problem.sol.K=K;
problem.semigroup.N=NQ;

% determine symmetry
if initialdatareal && contains(symmetry,'even')
    problem.symmetry='cosineseries';
elseif initialdataimag && contains(symmetry,'odd')
    problem.symmetry='sineseries';
else
    problem.symmetry='none';
end

% set scales
if ~exist('scale','var') || isempty(scale)
    scale.space=1;
    scale.time=1;
end
problem.scale=scale;
problem.scalefloats.space=altmid(problem.scale.space);
problem.scalefloats.time=altmid(problem.scale.time);

if isfield(griddata,'infinite') && griddata.infinite==true
    problem.infinitetime=true;
else
    problem.infinitetime=false;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% time grid %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(griddata,'variable') && griddata.variable
    % grid based on Z-error estimate
    problem.timegridfloats=generategrid(problem,griddata);
else
    % traditional constructed grid
    D=griddata.D-problem.infinitetime;
    T=altmid(griddata.T);
    if ~isfield(griddata,'skew')
        griddata.skew=0;
    end
    if ~isfield(griddata,'dip')
        griddata.dip=0;
    end
    s=linspace(0,1,D+1);
    % define the grid
    problem.timegridfloats=T*s.*(1+griddata.skew*(s-1)+griddata.dip*(s-1).*(s-0.5));
    if ~all((problem.timegridfloats(2:D+1)-problem.timegridfloats(1:D))>0) 
        disp('Negative step sizes: need to adjust skew-dip');
        error('Negative step sizes: need to adjust skew-dip');
    end
end
if altisintval(griddata.T)
    problem.timegrid=intval(problem.timegridfloats);
    % replace final gridpoint by the integration time as an interval
    problem.timegrid(end)=griddata.T;
else
    problem.timegrid=problem.timegridfloats;
end
if exist('intval.m','file')
    intvaltimegrid=intval(problem.timegrid);
    problem.domains=(intvaltimegrid(2:end)-intvaltimegrid(1:end-1));
else
    problem.domains=diff(problem.timegrid);
end

%%%%%% end time grid %%%%%%%

D=length(problem.timegrid)-1;
problem.sol.K(1:D)=K;

if problem.infinitetime
    D=D+1;
    problem.timegrid(D+1)=Inf;
    problem.domains(D)=Inf;
    problem.sol.K(D)=0;
end

problem.sol.D=D;

problem.proof.interpolation.K0=repmat(problem.proof.interpolation.K0,D,1);  
problem.proof.integrals.K1=repmat(problem.proof.integrals.K1,D,1);  

if isfield(griddata,'optimizecheb') && griddata.optimizecheb
   % optimize the number of Chebyshev modes used
   optimizedchebyshev=generatecheb(problem,griddata);
   problem.sol.K=optimizedchebyshev.K;
   problem.proof.interpolation.K0=optimizedchebyshev.K0;
   problem.proof.integrals.K1=optimizedchebyshev.K1;
end

% rstar may be different on each domain but if not, then we need to repeat it
if length(problem.proof.rstar)==1
    problem.proof.rstar=repmat(problem.proof.rstar,D,1);
end
if problem.infinitetime && isfield(problem.proof,'rstarinf')
    problem.proof.rstar(D)=problem.proof.rstarinf;
end

end