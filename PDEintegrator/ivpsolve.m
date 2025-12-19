function [y,problem] = ivpsolve(problem)
% solve the IVP numerically using an integrator 
% and then improve using Newton iterations
%
% The default is to use the more elaborate semigroup described in the paper.
% If problem.semigroup.type='naive' then 
% the problem setup uses the naive semigroup
% based on the highest order derivative only
% and ignoring any nonlinear terms.
%
% outputs symmetric data if the problem is symmetric

% numerical integation of IVP
disp('Numerical integration');
problemfloat=problem;
problemfloat.initial.data=altmid(problemfloat.initial.data);
problemfloat.timegrid=altmid(problemfloat.timegrid);
problemfloat.domains=altmid(problemfloat.domains);
x=numericalintegration(problemfloat);

if isfield(problem,'symmetry') 
    if strcmp(problem.symmetry,'cosineseries')
        for d=1:length(x)
            x{d}=real(x{d});
        end
    elseif strcmp(problem.symmetry,'sineseries')
        for d=1:length(x)
            x{d}=1i*imag(x{d});
        end
    end
end
% real-valued solution
for d=1:length(x)
    x{d}=(x{d}+conj(x{d}(end:-1:1,:)))/2;
end

% Now iterate the map with the semigroup
if ~isfield(problem.semigroup,'type') || ~isequal(problem.semigroup.type,'naive') 
    % Newton iteration while improving semigroup
    disp('Determine the semigroup');
    problemfloat.semigroup = fixsemigroup(x,problemfloat);
    y=findzero(x,problemfloat,true);
    % final iterations with fixed semigroup 
    disp('Update the semigroup');
    problemfloat.semigroup = fixsemigroup(y,problemfloat);
    [y,residue]=findzero(y,problemfloat,false);
else
    % replace with naive semigroup data
    disp('Use naive semigroup');
    problemfloat.semigroup=naivesemigroup(x,problemfloat);
    [y,residue]=findzero(x,problemfloat);
end

disp(['norm of the residue is ',num2str(residue)]);
problem.semigroup=problemfloat.semigroup;

end

function semigroup=naivesemigroup(x,problem)
% semigroup data for naive semigroup
    D=length(x);
    P=problem.pde.order;
    semigroup.N=0;
    semigroup.Q=ones([1,1,D]);
    semigroup.Qinv=ones([1,1,D]);
    semigroup.Lambda=zeros([1,D]);
    semigroup.ToepV=zeros([1,D,P]); 
    semigroup.Residue=zeros([1,1,D]);
end

