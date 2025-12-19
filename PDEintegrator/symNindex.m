function Nindex = symNindex(N,Nbase,symmetry)
% returns the range of the counter for different symmetries

if strcmp(symmetry,'cosineseries')
    Nindex = 1:N+1;
elseif strcmp(symmetry,'sineseries')
    Nindex = 1:N;
else
    Nindex = Nbase+1+(-N:N);
end

end