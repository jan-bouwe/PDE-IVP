function [success,rsol,bounds,errorbound]=ivpproof(x,problem,intervals)
% proof with or without intervalarithmetic
% if x(1) is an interval or intervals==true 
% then interval arithmetic is used if it is available

% decide whether to use intervals or not
useintervals=false;
if exist('intervals','var')
    if intervals==true 
        if exist('intval.m','file')
            useintervals=true;
        else
            disp(['No interval arithmetic available!! '...
                        'Continuing without intervals'])
        end

    end
elseif exist('intval.m','file') && altisintval(x(1))
    % default is to use intervals if x(1) is an interval
    useintervals=true;
end

if useintervals
    for d=1:length(x)   
        x{d}=intval(x{d});
    end
    disp('Full proof including interval arithmetic')
else
    for d=1:length(x) 
        x{d}=altmid(x{d});
    end
    disp('Proof attempt without interval arithmetic')
end

for d=1:length(x)
    xd=altmid(x{d});
    if strcmp(problem.symmetry,'cosineseries')
        symmetric=all(imag(xd)==0,'all') && isequal(xd,xd(end:-1:1,:));
    elseif strcmp(problem.symmetry,'sineseries')
        symmetric=all(real(xd)==0,'all') && isequal(xd,-xd(end:-1:1,:));
    else
        symmetric=true;
    end
    if ~symmetric
        error(['x is not in the ',problem.symmetry,' subspace on domain ',num2str(d)]);
    end
end

if ~isfield(problem,'semigroup') 
    problem.semigroup.type='elaborate';
end
% interval arithmetic semigroup data if x are intervals
problem.semigroup=fixsemigroup(x,problem);

% proof
%tic
[success,rsol,bounds,errorbound]=provesol(x,problem);
%toc

end