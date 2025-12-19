function symu = symQinvproject(u,Qinv,symmetry)
% projects u onto symmetry variables (left multiplication by Qinv)

Nu=(size(u,1)-1)/2;
NQ=(size(Qinv,2)-1)/2;

Nurange=symNrange(Nu,symmetry);
symu=u(Nu+1+Nurange,:);

NLsymindex=symNindex(NQ,Nu,symmetry);
symu(NLsymindex,:)=Qinv*u(Nu+1+(-NQ:NQ),:);

end