function [success,rsol,bounds,errorbound]=provesol(x,problem)
% Takes the approximate solution x and problem parameters and 
% attempts to prove a solution closeby exists
% If x consists of intervals all computations will be rigorously verified
% using interval arithmetic in Intlab 
% (otherwise it is based on floating point arithmetic)
% Here x is an (2N+1)x(K+1)xD tensor
% and we assume x is already conjugate symmetric in the Fourier component
%
% If the radii polynomial is negative for some positive r, 
% then success=true and rsol are the succesfull radii 
%
% Additionally, bounds contains some additional data about the bounds 
% (for saving) from which the radii polynomials have been built 
%
% If the method is unsuccesfull then success=0 and rsol=NaN

success=false;
rsol=NaN;

if altisintval(x{1}(1)) && ~altisintval(problem.semigroup.ToepV(1))
    % if data of semigroup is not of interval type
    % recompute semigroup to get exact inverses Q and Qinv
    problem.semigroup = fixsemigroup(x,problem);
end

% Now the interval arithmetic computations start
disp('Starting the proof')
problem.type='proof';
% turn all data into intervals if x is intval
if altisintval(x{1}(1))
    disp('Including interval arithmetic')
    problem.timegrid=intval(problem.timegrid);
    problem.initial.tailbound=intval(problem.initial.tailbound);
    problem.proof.nu=intval(problem.proof.nu);
    problem.proof.rstar=intval(problem.proof.rstar);
	if ~isintval(problem.initial.data)
		initialdataradii=1.414214*eps(problem.initial.data);
		problem.initial.data=midrad(problem.initial.data,initialdataradii);     
	end
else
    disp('No intervals yet');
    problem.initial.data=altmid(problem.initial.data);
    problem.timegrid=altmid(problem.timegrid);
    problem.domains=altmid(problem.domains);
end
disp(['The proof is performed with symmetry: ',problem.symmetry]);

% degree of the higher order polynomial which we use to approximate the interpolation error
problem.proof.interpolation.K0 = max(problem.proof.interpolation.K0,max(problem.sol.K));

% % the bound on the preconditioner introduced for the domain decomposition
% interpolationerrorQEQ = estimateG(x,problem);
% problem.interpolationerrorQEQ = interpolationerrorQEQ;

% get the Y bound
[~,Y] = getTYbounds(x,problem);
if isfield(problem,'displaybounds') && problem.displaybounds
    disp('The Y-bound:')
    disp(num2str(altsup(Y)));    
end
% disp(['Rough norm of the Y-bound: ',num2str(altsup(max(Y)))]);


% get the Z bound
[Z,W] = getZWbounds(x,problem);
if isfield(problem,'displaybounds') && problem.displaybounds
    disp('The Z-bound:');
    disp(num2str(altsup(Z)));  
end
% disp(['Dominant eigenvalue of the Z-matrix is roughly ',num2str(eigs(altsup(Z),[],1))]); 
% size of W bound: sum(pagenorm(W,Inf)
sizeW=sum(max(sum(altsup(W),2),[],1,'includenan'));
disp(['Rough size of W: ',num2str(sizeW)]); 

% find a point where the radii polynomials are negative
[rminvec,eta]=polynomialsnegative_ivp(Y,Z,W);
rstarvec=problem.proof.rstar;

if any(isnan(rminvec)) || ~all(rminvec>0)
    % failure :-(
    disp('Failed due not finding a positive rmin');

elseif all(rminvec <= rstarvec)
    % success :-)
    success=true;  
    rsol=rminvec;
else
    % failure :-(
    disp('Failed due to rmin not being smaller than rstar');
    disp(['rmin  = [' ,num2str(rminvec'),' ]']); 
    disp(['rstar = [' ,num2str(altsup(rstarvec')),' ]']); 
end

if success 
    % theta = max(max(problem.proof.eps0,problem.proof.eps1),1);
    % errorbound = altsup(max(theta.*rsol));
    disp('The proof was successful: there is a solution');
    if isfield(problem,'displaybounds') && problem.displaybounds
        disp(['The validation radii are ',num2str(rsol')]);
    end
    errorbound=altsup(max(rsol));
    if problem.infinitetime
        disp(['The uniform error bound on the finite part is ',num2str(altsup(max(rsol(1:end-1))))]);    
        disp(['The uniform error bound on the infinite domain is ',num2str(altsup(max(rsol(end))))]);    
    else        
        disp(['The uniform error is ',num2str(errorbound)]);          
    end
    if altisintval(x{1}(1))
        disp('Full proof including interval arithmetic')
    else
        disp('But this does not include interval arithmetic');
    end
else
    disp('The proof failed');
    rsol=NaN;
    errorbound=NaN;
end 

% storing some bounds
bounds.Y=altsup(Y);
bounds.Z=altsup(Z);
bounds.W=altsup(W);
bounds.eta=eta;

end


