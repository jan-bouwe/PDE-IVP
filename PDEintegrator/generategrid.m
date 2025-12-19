function optimizedtimegrid=generategrid(problem,griddata)
% generates a time grid based on the Z-bound
% The number of domains is approximately griddata.D
% and the total time is griddata.T
%
% The default is to start with a uniform grid, but the user may 
% optionally provide a griddata.startinggrid for the search
%
% In griddata the user can provide values to override several defaults values
% for parameters used in the algorithm
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set algorithmic parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we do new numerics once every 3 grids
updatenumericsfrequency=setdefault(3,griddata,'updatenumericsfrequency');
% and we only use the last 3 grids to find a new grid
updategridfrequency=setdefault(3,griddata,'updategridfrequency');
% limit extrapolation to large step sizes after the first 3 grids
limitextrapolationafter=setdefault(3,griddata,'limitextrapolationafter');
% part of the points where the high value bounds should be roughly uniform 
uniformpart=setdefault(0.75,griddata,'uniformpart');

% weight of taking difficulty of diagonalization into account 
difficultyfactorstart=setdefault(0.5,griddata,'difficultyfactorstart');
% increases by factor
difficultyfactorincrease=setdefault(1.3,griddata,'difficultyfactorincrease');  
% after 4 iterations
difficultyfactorincreasestart=setdefault(4,griddata,'difficultyfactorincreasestart');  
% with a maximum value of
difficultyfactormax=setdefault(2,griddata,'difficultyfactormax');
% the part for which difficulty gets capped to take out peaks
difficultyfactortoppart=setdefault(0.02,griddata,'difficultyfactortoppart');     
% maximum number of iterations to find a good grid
maxsteps=setdefault(12,griddata,'maxsteps');
% how uniform you would like the Zbound to be (relative error)
accuracy=setdefault(0.15,griddata,'accuracy');
% how much variation in the number of domains is allowed
D=griddata.D;
Dvar=setdefault(floor(0.05*D),griddata,'Dvar');

%%%%%%%%%%%%%%
% initialize %
%%%%%%%%%%%%%%

K=problem.sol.K;
problem.infinitetime=false;
problem.nodisplay=true;

T=altmid(griddata.T);
problem.initial.data=altmid(problem.initial.data);
chebgrid=(fliplr(cheb_grid(K))+1)/2;

disp('Generating grids to optimize Z-bound')

% generate initial grids
if isfield(griddata,'startinggrid')
    % the user may provide an initial guess for the search
    problem.timegrid=griddata.startinggrid;
else
    % the default is a uniform grid
    problem.timegrid=linspace(0,T,D+1);
end
problem.domains=diff(problem.timegrid);
problem.sol.D=length(problem.domains);
% solve numerically
problem.sol.K=repmat(K,D,1);
[x,alldata]=numericalintegration(problem);
problem.semigroup = fixsemigroup(x,problem);

% factor indicating how large the Q and Qinv matrices are
Qfactor=zeros([problem.sol.D,1]);
for d=1:problem.sol.D
    Qd=problem.semigroup.Q(:,:,d);
    Qinvd=problem.semigroup.Qinv(:,:,d);
    Qfactor(d)=norm(abs(Qd)*abs(Qinvd),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the first bounds for initial grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zaccuracy=accuracy;
Zbound=getZWbounds(x,problem,true);
timegrid{1}=problem.timegrid(1:end-1);
taugrid{1}=problem.domains;
Zgrid{1}=Zbound;
Zaverage=sum(Zbound)/length(Zbound);
Zvariation=max(abs(Zbound-Zaverage));
Zmax=max(Zbound);
difficultyweight=difficultyfactorstart;
disp('Generated testgrid 1');

%%%%%%%%%%%%%%%%%%%%%%%%
% generate other grids %
%%%%%%%%%%%%%%%%%%%%%%%%

ngrids=1;
maxstepsgrid=maxsteps;
lastngrids=min(updategridfrequency,maxstepsgrid);
% at each iteration it tries to find a good value for Zgoal
% that will lead to approximately D domains
% and maxstepsZgoal is the maximum number of steps used for that
maxstepsZgoal=10;

while Zvariation/Zaverage>Zaccuracy && ngrids<maxstepsgrid+1
    % still looking for better grids
    nstoredgrids=min(ngrids,lastngrids);
    tau=zeros(nstoredgrids+1,1);
    Ztau=zeros(nstoredgrids+1,1);
    Dnew=zeros(maxstepsZgoal,1);
    Zgoal=zeros(maxstepsZgoal,1);
    mstep=0;
    % first guess for Zgoal is the average
    Zgoal(1)=Zaverage;
    Dnewest=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adapting Zgoal until we find approximately D grid points %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while abs(Dnewest-D)>Dvar && mstep<maxstepsZgoal && (mstep<1 || Zgoal(mstep)~=Zgoal(mstep+1))
        % we continue changing Zgoal since the number of domains is not right yet
        mstep=mstep+1;
        T0=0;
        newgrid=T0;
        % normalizing the difficulty (due to size of Q and Qinv) factor 
        Qfactorsort=sort(Qfactor);
        minQfactor=Qfactorsort(round(uniformpart*length(Qfactor)));
        maxQfactor=Qfactorsort(round((1-difficultyfactortoppart)*length(Qfactor)));
        Qfactorcapped=min(Qfactor,maxQfactor);
        maxminQfactor=max(Qfactorcapped)-minQfactor;

        while T0<T % generate next grid point

            % take difficulty into account
            QfactorT0=interp1(problem.timegrid(1:end-1),Qfactorcapped,T0,'linear','extrap');
            QfactorT0=min(QfactorT0,maxQfactor);
            Tdifficulty=max(0,QfactorT0-minQfactor)/maxminQfactor;
            % adapt the goal downwards if the Qfactor is large
            ZgoalT0=Zgoal(mstep)/(1+difficultyweight*Tdifficulty);
            
            % use only last lastngrids grids for interpolation 
            jshift=max(0,ngrids-lastngrids);
            for j=1:nstoredgrids
                tau(j+1)=interp1(timegrid{j+jshift},taugrid{j+jshift},T0,'linear','extrap');
                Ztau(j+1)=interp1(timegrid{j+jshift},Zgrid{j+jshift},T0,'linear','extrap');
            end
            [sortZtau,sortI]=sort(Ztau);
            cummaxsorttau=cummax(tau(sortI));
            % guess for the step based on interpolating from previous data
            dT=interp1(sortZtau,cummaxsorttau,ZgoalT0,'linear','extrap');
            if ZgoalT0>sortZtau(end)
                % damping in case of extrapolation (versus interpolation):
                % stepsize not more than linear extrapolation between max values and the origin
                extrapolationslope=cummaxsorttau(end)/sortZtau(end);
                % for later grids we damp extra by a factor
                limitingfactor=1/(1+max(0,mstep-limitextrapolationafter));
                dTmax=extrapolationslope*(sortZtau(end)+limitingfactor*(ZgoalT0-sortZtau(end))); 
                dT=min(dT,dTmax); 
            end
            % add new gridpoint to the grid
            T0=T0+dT;
            newgrid=[newgrid,T0];
        end % all grid points determined

        % determine whether the final grid point falls outside the time domain
        discrepancy=(T0-T)/(T0-newgrid(end-1));
        if discrepancy>0.5
            % remove final point (sticks out too much beyond T)
            newgrid(end)=[];
        end
        % now record the number of grid points
        Dnewest=length(newgrid)-1;
        % include fractional value
        Dnew(mstep)=Dnewest-(newgrid(end)-T)/(newgrid(end)-newgrid(end-1)); 
        if mstep>1
            % new guess for goal for Z based on what we have so far
            Zgoal(mstep+1)=interp1(Dnew(1:mstep),Zgoal(1:mstep),D,'linear','extrap');
            if Zgoal(mstep+1)<=0
                % overwrite the extrapolated negative value
                Zgoal(mstep+1)=Zgoal(mstep)*Dnew(mstep)/D;
            end
        else
            Zgoal(2)=Zgoal(1)*Dnew(mstep)/D;
        end
    end % done with new grid, which has number of domains near D
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % we proceed to find the Zbounds on this new grid %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % rescale grid to cover the interval [0,T]   
    newgrid=newgrid*(T/newgrid(end));
    problem.timegrid=newgrid;
    problem.domains=diff(problem.timegrid);
    problem.sol.D=length(problem.domains);
    
    % update the numerical data for the approximate solution

    if mod(ngrids,updatenumericsfrequency)~=0
        % interpolate from the previously computed approx sol
        x=cell(problem.sol.D,1);
        problem.sol.K=repmat(K,problem.sol.D,1);
        for d=1:problem.sol.D
            times=problem.timegrid(d)+chebgrid*problem.domains(d);
            sol=interp1(alldata.times,alldata.orbit,times,'linear','extrap');
            x{d}=cheb_coeffs(fliplr(sol.'));
        end
    else
        % solve numerically (update the approximate solution)
        problem.sol.K=repmat(K,problem.sol.D,1);
        [x,alldata]=numericalintegration(problem);        
    end

    % update the semigroup
    problem.semigroup = fixsemigroup(x,problem);

    % compute the Zbounds
    Zbound=getZWbounds(x,problem,true);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % finalize some things before we can go to the next iteration step %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % recompute the difficulty = Qfactor
    Qfactor=zeros([problem.sol.D,1]);
    for d=1:problem.sol.D
        Qd=problem.semigroup.Q(:,:,d);
        Qinvd=problem.semigroup.Qinv(:,:,d);
        Qfactor(d)=norm(abs(Qd)*abs(Qinvd),1);
    end
    if ngrids>difficultyfactorincreasestart
        % increase difficultyfactor = how much we take Qfactor into account
        difficultyweight=difficultyfactorincrease*difficultyweight;
        % but never to more than difficultyfactormax
        difficultyweight=min(difficultyweight,difficultyfactormax);
    end

    ngrids=ngrids+1;
    timegrid{ngrids}=problem.timegrid(1:end-1);
    taugrid{ngrids}=problem.domains;
    Zgrid{ngrids}=Zbound;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % determine the relevant variation in the Z bound %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % to see how uniform the Zbounds are 
    % we look at variations in the part where the Qfactor is not large
    % which is also where Z will be large 
    % (since we did not lower the target there)
    Qfactorsort=sort(Qfactor);
    % uniformpart is the fraction of points we take into account
    % to determine the average of Z
    minQfactor=Qfactorsort(round(uniformpart*length(Qfactorsort)));

    ZsmallQfactor=Zbound(Qfactor<=minQfactor);
    Zaverage=sum(ZsmallQfactor)/length(ZsmallQfactor);
    % for the maximum error we do look at all points
    Zvariation=max(Zbound-Zaverage);
    Zmax(ngrids)=max(Zbound);

    disp(['Generated testgrid ',int2str(ngrids)]);

end % move on to next grid

%%%%%%%%%%
% wrapup %
%%%%%%%%%%

% select the best grid rather than the last grid
[~,bestgrid]=min(Zmax);
optimizedtimegrid=[timegrid{bestgrid},T];
Dnew=length(optimizedtimegrid)-1;

maxZbound=Zmax(bestgrid);
disp(['Best generated grid has ',int2str(Dnew), ' domains and approximate Z-bound ',num2str(maxZbound)]);

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