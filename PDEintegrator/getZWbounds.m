function [Z,W] = getZWbounds(x,problem,diagonalZonly)
% computes the Z and W bounds
% If diagonalZonly is true then only the Z bound on the diagonal is computed 
% (this is useful for determining a grid, see generategrid.m)

%%%%%%%%%%%%%%%%%%%
%  Initialization %
%%%%%%%%%%%%%%%%%%%

if nargin<3
    diagonalZonly=false;
end

NQ=problem.semigroup.N;
NX=problem.sol.N;
D=problem.sol.D;
symmetry=problem.symmetry;

nu = problem.proof.nu;
tau = problem.domains/2;
P = length(problem.pde.polynomials);
Dphifull = getDphifull(x,problem);
Kmax=zeros(D,1);
for d=1:D
    for p=1:P
        Kmax(d) =max(Kmax(d),(size(Dphifull{p,d},2)-1)/2);
    end
end

Nbar = max(1,nonlinmaxdegree(problem)-1)*NX;

% NG denotes the size of part of G needed/used for the computations
NG = 2*Nbar+NQ; 
indQ = NG+1+(-NQ:NQ);
indnotQ=1:2*NG+1; 
indnotQ(indQ)=[];
NGsym=length(symNrange(NG,symmetry));

der = 1i*(-NG:NG)';
derabs = abs(-NQ:NQ);

ZFinite = altzeros([D,D],x{1}(1));
ZTail = altzeros([D,D],x{1}(1));
W = altzeros([D,D,D],x{1}(1));

weights = nu.^abs(-NG:NG);
rstar=problem.proof.rstar;

% Precomputing real parts of eigenvalues
lambda = computelambdatail(NG,problem,'all');
indL=symNindex(NQ,NG,symmetry);
NL=length(indL);
NLmax=size(lambda,1);
indnotL=1:NLmax;
indnotL(indL)=[];
lambda(indL,:) = problem.semigroup.Lambda;
rlambda = real(lambda);

% The code computes each B^{(m,i)} when needed, 
% which is best memory wise but not speed wise.

for dd = 1:D % dd is i in the paper
    if ~diagonalZonly
        disp(['Computing the Z-bound on domain ',num2str(dd)])
    end

    Qdd=problem.semigroup.Q(:,:,dd);
    Qinvdd=problem.semigroup.Qinv(:,:,dd);
    % Constructing Gamma^{(i)}
    GammaT = altzeros([2*NG+1,2*NG+1,Kmax(dd)+1],x{1}(1));  
    QinvGammaT = altzeros([NGsym,2*NG+1,Kmax(dd)+1],x{1}(1));  
    ZTaildiag = 0;
    for p = 1:P
        if ~isempty(problem.pde.polynomials{p})

            Dphip = Dphifull{p,dd};
            Np =(size(Dphip,1)-1)/2;
            Kp =(size(Dphip,2)-1)/2;
            % full derivative
            gtilde = Dphip(:,Kp+1:2*Kp+1);
            Vp = problem.semigroup.ToepV(:,dd,p);
            % alpha and beta
            indQp = Np+1+(-NQ:NQ);
            alpha = gtilde(indQp,:);
            alpha(:,1) = alpha(:,1) - Vp;
            beta = gtilde;
            beta(Np+1,1) = beta(Np+1,1) - Vp(NQ+1);
            normC0beta = normC0(beta);
    
            % finite part 
            for k=1:Kp+1
                alphak = alpha(:,k);
                betak = settensorsize(beta(:,k),[NG,0]);
                Mpk = toeplitzshift(betak);
                Mpk(indQ,indQ) = toeplitzshift(alphak);
                Mpk = der.^(p-1) .* Mpk;
                GammaT(:,:,k) = GammaT(:,:,k) + Mpk;
            end
    
            % tail part 
            ZTaildiag = ZTaildiag + normL1(normC0beta,nu)*boundchi(p-1,NQ,dd,problem);
        end
    end

    % adding the R_N term
    GammaT(indQ,indQ,1) = GammaT(indQ,indQ,1) + problem.semigroup.Residue(:,:,dd);

    % take into account symmetry, multiply by Qinv
    for k=1:Kmax(dd)+1
        QinvGammaT(:,:,k) = symQinvproject(GammaT(:,:,k),Qinvdd,symmetry);
    end

    % take C0 norm
    QinvGamma=normC02(QinvGammaT);

    % Precomputations for W
    if ~diagonalZonly
        normxrstar = normL1(normC0(x{dd}),nu) + rstar(dd);
        normDDg = altzeros([P,1],x{1}(1));
        chi = altzeros([P,1],x{1}(1));
        for p=1:P 
            % bound on norm of second derivative of nonlinearities
            normDDg(p) = nonlinpoly(normxrstar,p-1,2,problem,'absolute');
            chi(p) = boundchi(p-1,NG-2*Nbar,dd,problem);
        end
    end
    
    % Constructing B^{(m,i)}=B^{(d,dd)}. 
    % We build a proper matrix for the finite part, 
    % and only a vector containing the diagonal elements for the tail part. 
    % The computation is split between three terms, 
    % the middle one being reused as d varies.

    if ~(problem.infinitetime && dd==D)
        BQ_end = diag(expm1div(2*tau(dd)*rlambda(indL,dd))*2*tau(dd));
        BTail_end = expm1div(2*tau(dd)*rlambda(indnotL,dd))*2*tau(dd); 
    else
        % infinite domain
        if all(rlambda(:,D)<0)
            BQ_end = diag(-1./rlambda(indL,dd));
            BTail_end = -1./rlambda(indnotL,dd);
        else
            BQ_end = diag(Inf(length(indL),1));
            BTail_end = Inf(length(indnotL),1);
        end
    end
    if diagonalZonly
        dset=dd;
    else
        dset = dd:D;
    end
    for d = dset % d is m in the draft    
        Qd=problem.semigroup.Q(:,:,d);
        Qinvd=problem.semigroup.Qinv(:,:,d);
        if d == dd
            BQ_start = abs(Qd);
            BTail_start = ones(NLmax-NL,1);    
            BQ_mid = eye(NL);
            BTail_mid = ones(NLmax-NL,1);
        else
            if ~(problem.infinitetime && d==D)
                BQ_start = abs(Qd)*diag(exp(2*tau(d)*max(rlambda(indL,d),0)));
                BTail_start = exp(2*tau(d)*max(rlambda(indnotL,d),0));
            else
                % infinite domain
                if all(rlambda(:,D)<0)
                    BQ_start = abs(Qd);
                    BTail_start = ones(length(indnotL),1);
                else
                    BQ_start = abs(Qd)*diag(Inf(length(indL),1));
                    BTail_start = Inf(length(indnotL),1);
                end
            end

            if d == dd+1
                BQ_mid = Qinvd*Qdd;
                BTail_mid = ones(NLmax-NL,1);
            else
                Qdm1=problem.semigroup.Q(:,:,d-1);
                BQ_mid = (Qinvd*Qdm1)*diag(exp(2*tau(d-1)*lambda(indL,d-1)))*BQ_mid;
                BTail_mid = exp(2*tau(d-1)*lambda(indnotL,d-1)).*BTail_mid;
            end
        end
        % We now compute the actual B^{(m,i)}
        BQ = BQ_start * abs(BQ_mid) * BQ_end;
        BnotL = BTail_start .* abs(BTail_mid) .* BTail_end;    
        
        % Computing B^{(m,i)} * \Gamma^{(i)} 
        BGamma = altzeros([2*NG+1,2*NG+1],x{1}(1)); 
        BGamma(indQ,:) = BQ * QinvGamma(indL,:);

        if any(strcmp(symmetry,{'cosineseries','sineseries'}))
            DBnotL=diag(BnotL);
            BnotQ=[DBnotL(end:-1:1,:);DBnotL];
        else
            BnotQ=diag(BnotL);
        end
        BGamma(indnotQ,:) = BnotQ * QinvGamma(indnotL,:);        

        ZFinite(d,dd) = max((weights*BGamma)./weights);
        mu = exp(boundmu(d,dd,NQ,problem));
        ZTail(d,dd) = mu * ZTaildiag;  

        % W bound
        if ~diagonalZonly
            normBD = altzeros([P,1],x{1}(1));
            BQW=BQ*abs(Qinvdd);
            for p=1:P 
                normBD(p) = max(max((weights(indQ)*(BQW*diag(derabs.^p)))./weights(indQ)),chi(p)*mu);
            end
            W(d,dd,dd) = sum(normBD.*normDDg);
        end
    end
end

% compute total Z
Z = max(ZFinite,ZTail);

if diagonalZonly
    Z=diag(Z);
else
    % display
    disp(['size of Z is ',num2str(altsup(max(diag(Z))))]);
end

end
