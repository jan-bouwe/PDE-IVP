function u = symrightQinvproject(symu,Qinv,symmetry)
% right multiplication by Qinv

if strcmp(symmetry,'cosineseries')
    u=[symu(:,end:-1:2),symu];
elseif strcmp(symmetry,'sineseries')
    u=[-symu(:,end:-1:1),zeros([size(symu,1),1]),symu];
else
    u=symu;
end

Nu=(size(u,2)-1)/2;
NQ=(size(Qinv,2)-1)/2;

NLsymindex=symNindex(NQ,Nu,symmetry);
u(:,Nu+1+(-NQ:NQ))=symu(:,NLsymindex)*Qinv;

end

