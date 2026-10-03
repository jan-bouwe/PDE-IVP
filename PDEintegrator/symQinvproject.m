function symu = symQinvproject(u,Qinv,symmetry)
% projects u onto symmetry variables (left multiplication by Qinv)

Nu=(size(u,1)-1)/2;
NQ=(size(Qinv,2)-1)/2;

Nurange=symNrange(Nu,symmetry);
symu=u(Nu+1+Nurange,:);
% tail: average the n and -n coefficients, as Qinv does,
% so that Z is the operator norm on the symmetric subspace
if strcmp(symmetry,'cosineseries')
    symu(2:end,:)=(symu(2:end,:)+u(Nu:-1:1,:))/2;
elseif strcmp(symmetry,'sineseries')
    symu=(symu-u(Nu:-1:1,:))/2;
end

NLsymindex=symNindex(NQ,Nu,symmetry);
symu(NLsymindex,:)=Qinv*u(Nu+1+(-NQ:NQ),:);

end