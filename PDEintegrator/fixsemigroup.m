function semigroup = fixsemigroup(x,problem)
% semigroup = fixsemigroup(x,problem)
%
% computes both the semigroup data (eigenvalues,eigenvectors,inverse)
% semigroup.Q
% semigroup.Lambda
% semigroup.Qinv
%
% of the Toeplitz/convolution matrices (one per domain)
% built up from the vectors (one per nonlinearity) 
% semigroup.ToepV
% which is real or purely imaginary if the symmetry imposes it
% 
% as well as the residue of this approximate diagonalization:
% semigroup.Residue=Q*Lambda*Qinv-Toeplitz(ToepV)
%
% and the size of the matrices
% semigroup.N
%
% These may be of interval type or float type
% depending on the type of x
% and the inverse is an outer approximation in the interval case

NQ=problem.semigroup.N;
D=problem.sol.D;
K=problem.sol.K;
P=problem.pde.order;
useintervals=isa(x{1}(1),'intval');

semigroup.N=NQ;
semigroup.Q=altzeros([2*NQ+1,2*NQ+1,D],x{1}(1));
semigroup.Qinv=altzeros([2*NQ+1,2*NQ+1,D],x{1}(1));
semigroup.Lambda=altzeros([2*NQ+1,D],x{1}(1));
semigroup.ToepV=zeros([4*NQ+1,D,P]);
DNQ=(1i*(-NQ:NQ));
for d=1:D
    %the mean value in the domain
    %average=altmid(x(:,1,d));
    % extend x so that we can use convolutions to compute nonlinearities
    xx=altmid([x{d}(:,K(d)+1:-1:2),x{d}]);
    for p=0:problem.pde.order-1
        %for each nonlinearity determine the linearization at the mean value
        %averagenonlinearity=convnonlinearity(average,p,1,problem);
        % for each nonlinearity determine the mean of the derivative 
        [~,nonlinearity]=convnonlinearity(xx,p,1,problem);
        averagenonlinearity=nonlinearity(:,(size(nonlinearity,2)+1)/2);
        averagenonlinearity=settensorsize(averagenonlinearity,2*NQ);
        semigroup.ToepV(:,d,p+1)=averagenonlinearity;
    end
end
% symmetrize
semigroup.ToepV=...
    (semigroup.ToepV+conj(semigroup.ToepV(end:-1:1,:,:)))/2;    
if strcmp(problem.symmetry,'cosineseries')
    semigroup.ToepV=real(semigroup.ToepV);
elseif strcmp(problem.symmetry,'sineseries')
    semigroup.ToepV(:,:,1:2:P)=real(semigroup.ToepV(:,:,1:2:P));
    semigroup.ToepV(:,:,2:2:P)=1i*imag(semigroup.ToepV(:,:,2:2:P));
end
for d=1:D
   % LN is the linearization at the mean value (tilde L in the paper)
   LN=diag(-(abs(-NQ:NQ)').^problem.pde.order);
   for p=1:P
       LN=LN+diag(DNQ.^(p-1))*toeplitzshift(semigroup.ToepV(:,d,p),NQ); 
   end
   
   % approximate diagonalization
   % use some perturbations to make the algorithm more stable
   perturbations=[0,1,100,10000,1000000]*eps;
   sizeQLQinv=zeros(size(perturbations));
   for m=1:length(perturbations)
       [mQ,mL]=symeig(LN,problem.symmetry,perturbations(m),false);
       sizeQLQinv(m)=norm(abs(mQ)*diag(abs(mL))*abs(inv(mQ)),1);
   end
   % pick the best one (with a penalty for the perturbation)
   [~,m0]=min(sizeQLQinv+(2*NQ+1)*abs(perturbations));

   % diagonalize
   [Evectors,Evalues,invEvectors]=symeig(LN,problem.symmetry,perturbations(m0),useintervals);
   semigroup.Q(:,:,d)=Evectors;
   semigroup.Lambda(:,d)=Evalues;
   semigroup.Qinv(:,:,d)=invEvectors;
end

% in case of interval arithmetic, convert to intervals
if useintervals
    semigroup.ToepV=intval(semigroup.ToepV);
end

% now determine the residue of the approximate diagonalization
% with interval arithmetic if appropriate
semigroup.Residue=altzeros([2*NQ+1,2*NQ+1,D],x{1}(1));

for d=1:D
    % LN is the linearization at the mean value (tilde L in the paper)
    LN=diag(-(abs(-NQ:NQ)').^problem.pde.order);
    for p=1:P
        LN=LN+diag(DNQ.^(p-1))*toeplitzshift(semigroup.ToepV(:,d,p),NQ); 
    end
    semigroup.Residue(:,:,d)=...
            semigroup.Q(:,:,d)*diag(semigroup.Lambda(:,d))*semigroup.Qinv(:,:,d)-LN;
end

% Enforce symmetry
if strcmp(problem.symmetry,'cosineseries')
     semigroup.Q(:,NQ+2:2*NQ+1,:)=[];
     semigroup.Lambda(NQ+2:2*NQ+1,:)=[];
     semigroup.Qinv(NQ+2:2*NQ+1,:,:)=[];
elseif strcmp(problem.symmetry,'sineseries')
     semigroup.Q(:,1:NQ+1,:)=[];    
     semigroup.Lambda(1:NQ+1,:)=[];
     semigroup.Qinv(1:NQ+1,:,:)=[];  
end

end

function [Evectors,Evalues,invEvectors]=symeig(LN,symmetry,perturbation,useintervals)
% produces symmetrized eigenvectors 
% and their inverses with interval arithmetic if requested
N=(size(LN,1)-1)/2;
if any(strcmp(symmetry,{'sineseries','cosineseries'}))
    % symmetric case
    % build matrices on symmetry spaces
    Lcos=LN(N+1:2*N+1,N+1:2*N+1);
    Lcos(:,2:N+1)=Lcos(:,2:N+1)+LN(N+1:2*N+1,N:-1:1);
    Lsin=LN(N+2:2*N+1,N+2:2*N+1)-LN(N+2:2*N+1,N:-1:1);
    % compute their eigeninformation
    [Evec_cos,Eval_cos]=symeig(Lcos,'none',perturbation,useintervals);
    [Evec_sin,Eval_sin]=symeig(Lsin,'none',perturbation,useintervals);
    % convert back to full space
    Evec_cos=[Evec_cos(end:-1:2,:);Evec_cos];
    Evec_sin=[-Evec_sin(end:-1:1,:);zeros([1,N]);Evec_sin];
    Evectors=[Evec_cos,Evec_sin];
    Evalues=[Eval_cos;Eval_sin];
else
    % without symmetry
    [Evectors,Evalues]=eig(LN+perturbation);
    Evalues=diag(Evalues);
    if useintervals
        Evalues=intval(Evalues);
        Evectors=intval(Evectors);
    end
end
if nargout==3
   invEvectors=inv(Evectors);
end

end
