function optimizedchebyshev=generatecheb(problem,griddata)
% chooses numbers of Chebyshev nodes in the solution
% for each subdomain, based on an estimate of the Ybound
%
% In particular:
% optimizedchebyshev.K (problem.sol.K) 
%          the number of nodes in the approximate solution
% optimizedchebyshev.K0 (problem.proof.interpolation.K0)
%          the number of nodes in the interpolation of the output of the map  
% optimizedchebyshev.K1 (problem.proof.integrals.K1)
%          the number of nodes in the computation of an exponential integral
% The latter two are chosen to be equal in this algorithm
% 
% In griddata the user can provide values to override several defaults values
% for parameters used in the algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set algorithmic parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% target size of the Ybound
Ybound1target=setdefault(1e-6,griddata,'Ybound1target');
% target size of the interpolation error, much smaller than Ytarget
Ybound2target=setdefault(1e-8,griddata,'Ybound2target');
% minimal K = number of Chebyshev modes (K+1 really)
Kmin=setdefault(1,griddata,'Kmin');
% maximal K
Kmax=setdefault(20,griddata,'Kmax');
% minimal number of interpolation nodes K0
K0min=setdefault(30,griddata,'K0min');
% maximal number of interpolation nodes K0
K0max=setdefault(500,griddata,'K0max');
% step size for the optimization of K0
K0step=setdefault(20,griddata,'K0step');

%%%%%%%%%%%%%%
% initialize %
%%%%%%%%%%%%%%

D = problem.sol.D;
N = problem.sol.N;

if problem.infinitetime
    D=D-1;   
end

problem.timegrid=altmid(problem.timegrid);
problem.initial.data=altmid(problem.initial.data);
problem.domains=altmid(problem.domains);
problem.nodisplay=true;

disp('Optimizing number of Chebyshev nodes on each domain')

% find numerical solution
x = numericalintegration(problem);   
problem.semigroup = fixsemigroup(x,problem);
x = findzero(x,problem);

% define problem on single domain
pro1=problem;
pro1.infinitetime=false;
pro1.sol.D=1;

newK=zeros(D,1);
newK0=zeros(D,1);

%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over the domains %
%%%%%%%%%%%%%%%%%%%%%%%%%

for d=1:D
    % initialize
    xd{1}=x{d};
    pro1.initial.data=cheb_eval(xd{1},-1);
    if d==1
        pro1.initial.data=problem.initial.data;
    end
    
    % starting values for K and K0
    pro1.proof.interpolation.K0=problem.proof.interpolation.K0(d);
    pro1.sol.K=problem.sol.K(d);
    if d>1 && flagsuccessK 
        % overwrite with data from previous domain
        % if previous step was successful
        pro1.proof.interpolation.K0=newK0(d-1);
        if abs(newK(d-1)-problem.sol.K(d))>0
            pro1.sol.K=newK(d-1);
        end    
    end

    % adapt size 
    if size(xd{1},2)<pro1.sol.K+1
        xd{1}(:,pro1.sol.K+1)=zeros(2*N+1,1);
    elseif size(xd{1},2)>pro1.sol.K+1
        xd{1}(:,pro1.sol.K+2:end)=[];            
    end
    
    % K0 should not be smaller than K
    K0minK=max(K0min,Kmax);
    pro1.proof.interpolation.K0=max(pro1.proof.interpolation.K0,K0minK);
    pro1.proof.integrals.K1=pro1.proof.interpolation.K0;

    pro1.timegrid=problem.timegrid(d:d+1);
    pro1.timegridfloats=problem.timegridfloats(d:d+1);
    pro1.domains=problem.domains(d);

    % recompute (K and initial data have changed)
    pro1.semigroup = fixsemigroup(xd,pro1);
    xd = findzero(xd,pro1);

    % initialize logicals
    flagsuccessK0=false;
    flagnonsuccessK0=false;
    flagsuccessK=false;
    flagnonsuccessK=false;
    stopiteratingK0=false;
    stopiteratingK=false;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iteration process to optimize K and K0 (and K1=K0) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while ~(stopiteratingK0 && stopiteratingK)

        % get bounds
        [~,~,Yextra]=getTYbounds(xd,pro1);
        successK = (Yextra.bounds1 < Ybound1target);
        successK0 = (Yextra.bounds2 < Ybound2target);
        
        reachedK0min = (pro1.proof.interpolation.K0==K0minK);
        reachedK0max = (pro1.proof.interpolation.K0==K0max);
        reachedKmin = (pro1.sol.K==Kmin);
        reachedKmax = (pro1.sol.K==Kmax);

        if ~stopiteratingK
            % we are still iterating K
            if successK0 
                % interpolation can be trusted
                % hence we can adapt K
                if successK
                    % maybe K can be even smaller
                    flagsuccessK=true;
                    Kvalue=pro1.sol.K;
                    K0value=pro1.proof.interpolation.K0;
                    pro1.sol.K=max(pro1.sol.K-1,Kmin);
                    if size(xd{1},2)>pro1.sol.K+1
                        xd{1}(:,pro1.sol.K+2)=[];
                        pro1.semigroup = fixsemigroup(xd,pro1);
                        xd = findzero(xd,pro1);
                    end
                else
                    % K must be larger (unless Kmax already reached)
                    flagnonsuccessK=true;
                    pro1.sol.K=min(pro1.sol.K+1,Kmax);
                    if size(xd{1},2)<pro1.sol.K+1
                        xd{1}(:,pro1.sol.K+1)=zeros(2*N+1,1);
                        pro1.semigroup = fixsemigroup(xd,pro1);
                        xd = findzero(xd,pro1);
                    end
                end
                % we stop iterating K if we have seen both success and nonsuccess
                % or if we reach minimum or maximum values
                stopiteratingK = (flagsuccessK || reachedKmax) && (flagnonsuccessK || reachedKmin);
                
                % if we stop iterating K there are some things to take care
                % of before moving on
                if stopiteratingK
                    flagsuccessK0=true;            
                    if flagsuccessK
                        pro1.sol.K=Kvalue;
                        if size(xd{1},2)<Kvalue+1
                            xd{1}(:,Kvalue+1)=zeros(2*N+1,1); 
                            pro1.semigroup = fixsemigroup(xd,pro1);
                            xd = findzero(xd,pro1);
                        end
                        % K0 should not be smaller than K
                        K0minK=max(K0min,Kvalue);
                        pro1.proof.interpolation.K0=max(K0value-K0step,K0minK);
                        pro1.proof.integrals.K1=pro1.proof.interpolation.K0;
                    else % stop due to reaching Kmax
                        Kvalue=Kmax;
                        K0value=pro1.proof.interpolation.K0;
                        K0minK=max(K0min,Kvalue);
                        pro1.proof.interpolation.K0=max(pro1.proof.interpolation.K0-K0step,K0minK);
                        pro1.proof.integrals.K1=pro1.proof.interpolation.K0;
                    end
                end
            elseif ~reachedK0max
                % cannot trust interpolation, increase K0
                pro1.proof.interpolation.K0=min(pro1.proof.interpolation.K0+K0step,K0max);
                pro1.proof.integrals.K1=pro1.proof.interpolation.K0;
            else
                % still unsuccessful at maximum value of K0
                disp('Please increase K0max to optimize for K');
                K0value=K0max;
                Kvalue=problem.sol.K(d);
                stopiteratingK=true;
                stopiteratingK0=true;
            end
        else % we have already stopped iterating K and are now optimizing K0
            if successK0 && successK
                % maybe K0 can be even smaller
                flagsuccessK0=true;
                K0value=pro1.proof.interpolation.K0;
                pro1.proof.interpolation.K0=max(pro1.proof.interpolation.K0-K0step,K0minK);
                pro1.proof.integrals.K1=pro1.proof.interpolation.K0;
            else
                % K0 is now too small (or K1 is too small)
                flagnonsuccessK0=true;
                pro1.proof.interpolation.K0=min(pro1.proof.interpolation.K0+K0step,K0max);
                pro1.proof.integrals.K1=pro1.proof.interpolation.K0;
            end
            % we stop iterating K0 if we have seen both success and nonsuccess
            % or if we reach minimum or maximum values
            stopiteratingK0 = (flagsuccessK0 || reachedK0max) && (flagnonsuccessK0 || reachedK0min);
        end
                
    end

    newK(d)=Kvalue;
    newK0(d)=K0value;    
    disp(['Optimized number of Chebyshev nodes on domain ',num2str(d)])

end 

% finalize the output
if problem.infinitetime
    newK(D+1)=0;
    newK0(D+1)=1;
end

optimizedchebyshev.K=newK;
optimizedchebyshev.K0=newK0;
optimizedchebyshev.K1=newK0;

disp('Finished optimizing number of Chebyshev nodes')

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value=setdefault(defaultvalue,variablename,fieldname)
% a function to set default values that can be overwritten easily
    if isfield(variablename,fieldname)
        value=variablename.(fieldname);
    else
        value=defaultvalue;
    end
end