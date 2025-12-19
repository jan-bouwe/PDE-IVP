function [ y ] = normC02(M,rho)
% normC02(M) computes an upperbound(!) for the C0-norm on [-1,1]
% hence the name of the function is a bit misleading
% Here M is assumed to be in Fourier-Fourier-Chebyshev format
% normC02(M,rho) does the same but for the Bernstein-ellipse of size rho.

if ~exist('rho','var')
    rho=1;
end

S=[size(M),1];
M=reshape(M,[S(1)*S(2),S(3)]);
y=normC0(M,rho);
y=reshape(y,[S(1),S(2)]);

end

