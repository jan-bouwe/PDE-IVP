function Nrange = symNrange(N,symmetry)
% returns the Fourier index set for different symmetries

if strcmp(symmetry,'cosineseries')
    Nrange = 0:N;
elseif strcmp(symmetry,'sineseries')
    Nrange = 1:N;
else
    Nrange = -N:N;
end

end