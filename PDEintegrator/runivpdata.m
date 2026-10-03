function [x,problem,success,rsol,bounds,errorbound] = runivpdata(filename,datadirectory,resultdir)
% run the initial value problem with the name filename
% which has all data precomputed and stored the datadirectory 
% 
% save results (data and figure) in the resultdir (default 'results')

if ~exist('resultdir','var') 
    dirname='results/';
elseif isempty(resultdir)
    dirname=[];
else
    dirname=[resultdir,'/'];
end

intervalsavailable = exist('intval.m','file');
if intervalsavailable
    datatype='intval';
else
    datatype='float'; 
end

% get data
load([datadirectory,'/data',datatype,filename,'.mat'],'x','problem');
plotdata=[];

% floats or intervals or both
if ~isfield(problem,'useintervals') || isempty(problem.useintervals) || ...
                                        ~islogical(problem.useintervals)
    problem.useintervals=[];
    tryboth=true;
else
    tryboth=false;
end
tryintervals = (tryboth || problem.useintervals);
tryfloats = (tryboth || ~problem.useintervals);
if ~intervalsavailable && ~tryboth && tryintervals
    disp('Interval arithmetic is not available!')
    disp('Continue with floats instead.')
    tryintervals=false;
    tryfloats=true;
end

% make a figure of the solution
plotsolution(x,problem,plotdata,[dirname,'figure',filename]);
% save the data
save([dirname,'data',filename,'.mat'],'x','problem');

%start proofs
if tryfloats
    % using floats
    useintervals=false;
    tic
    [success,rsol,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    if success 
        save([dirname,'dataprooffloat',filename,'.mat'],'x','problem',...
              'success','rsol','bounds','errorbound','useintervals','tictoc');
    else
        disp('Proof with floats failed')
        tryintervals=false;
    end
end

if ~intervalsavailable && tryintervals && tryboth
    disp('Interval arithmetic is not available!')
    disp('Cannot continue with intervals.')
    tryintervals=false;
end

% now using interval arithmetic
if tryintervals
    useintervals=true;
    tic
    [success,rsol,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    if success 
        save([dirname,'dataproofintval',filename,'.mat'],'x','problem',...
              'success','rsol','bounds','errorbound','useintervals','tictoc');
    else
        disp('Proof with intervals failed')
    end
end


end