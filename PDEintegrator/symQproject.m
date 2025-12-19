function u = symQproject(symu,Q,symmetry)
% projects symu onto full variables (left multiplication by Q)

if strcmp(symmetry,'cosineseries')
    u=[symu(end:-1:2,:);symu];
elseif strcmp(symmetry,'sineseries')
    u=[-symu(end:-1:1,:);zeros([1,size(symu,2)]);symu];
else
    u=symu;
end

Nu=(size(u,1)-1)/2;
NQ=(size(Q,1)-1)/2;

NLsymindex=symNindex(NQ,Nu,symmetry);
u(Nu+1+(-NQ:NQ),:)=Q*symu(NLsymindex,:);

end