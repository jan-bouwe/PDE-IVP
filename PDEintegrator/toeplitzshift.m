function matrix = toeplitzshift(v,N)
% toeplitz matrix with v in the middle column
% of size (2N+1)x(2N+1), with entries v_{n-m} for |n|,|m|<=N
% the default for N is (length(v)-1)/2

Nv=(length(v)-1)/2;
if ~exist('N','var')
    N=Nv;
end
M=min(Nv,2*N);
vv=altzeros([4*N+1,1],v(1));
vv(2*N+1+(-M:M))=v(Nv+1+(-M:M));
column=vv(2*N+1:4*N+1);
row=vv(2*N+1:-1:1).';
matrix=toeplitz(column,row);

end

