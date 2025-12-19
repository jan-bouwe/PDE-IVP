function symu = symrightQproject(u,Q,symmetry)
% right multiplication by Q

Nu=(size(u,2)-1)/2;
NQ=(size(Q,1)-1)/2;

Nurange=symNrange(Nu,symmetry);
symu=u(:,Nu+1+Nurange);

NLsymindex=symNindex(NQ,Nu,symmetry);
symu(:,NLsymindex)=u(:,Nu+1+(-NQ:NQ))*Q;

end

