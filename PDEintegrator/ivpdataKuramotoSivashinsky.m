function [problem,plotdata] = ivpdataKuramotoSivashinsky
% data for Kuramoto-Sivashinsky equation
%
% u_t = - u_{xxxx} - u_{xx} - 1/2(u^2)_x
%
% for x in [0,L] and t in [0,tau]
%
% returns the problem definition, including the initial data 
% with the initial data setting the number N of Fourier modes
% the number K of Chebyshev modes (interpolation nodes) in time
% the number D of domains in the (time) domain decomposition
% (note that D is called M in the latex)
%
% the problem definition also includes some constants needed in the proof
%
% additional data for the plot can also be specified

if exist('intval.m','file')
    altpi=intval('pi');
    altone=intval(1);
else
    altpi=pi;
    altone=1;
    epsfactor=[];
end

%% Set problem parameters %%

% define the constants in the equation
% u_t = - u_{xxxx}  - u_{xx} - 1/2(u^2)_x
% for x in [0,L] and t in [0,tau]

alpha = 127*altone/1000;

L=2*altpi/sqrt(alpha);  
tau=(2245*altone/1000)/alpha;

% rescale to spatial domain [0,2*pi] and
% u_t = - u_{xxxx}  - lambda2*u_{xx} - lambda1*(u^2)_x

scale.space=L/(2*altpi);
scale.time=scale.space^4;
integrationtime=tau/scale.time;

% lambda2=scale.time/scale.space^2;
% lambda1=scale.time/scale.space/2;
lambda2=1/alpha;
lambda1=lambda2*scale.space/2;

%% Truncation and grid %%

% number of Fourier modes (2N+1 really)
N=30; 
% number of Chebyshev modes (K+1 really)
K=10; 
% number of time domains (called M in the latex)
D=40; 

% size of matrices in the semigroup
NQ=20;

% total integration time
griddata.T=integrationtime;
griddata.D=D;

 % we let the code pick the subdivision
griddata.variable=true;
% and set some algorithmic parameters for finding this subdivision
griddata.accuracy=0.1;
griddata.Dvar=1;
griddata.maxsteps=10;

% optimize the number K of chebyshev modes and K0 used for approximating
% via interpolation, in each subdomain
griddata.optimizecheb=true;
% and set some algorithmic parameters for this
griddata.Ybound1target=1e-12;
griddata.Ybound2target=1e-14;
griddata.Kmin=1;
griddata.Kmax=100;
griddata.K0min=30;
griddata.K0max=500;
griddata.K0step=20;

%% Initial data %%

% set initial data 
load('data/initialdataKS.mat','xinit')
Ninit = (length(xinit)-1)/2;
initialdata = altzeros([2*N+1,1],altone);
if N>= Ninit
    initialdata(N+1+(-Ninit:Ninit)) = xinit;
else
    initialdata = xinit(Ninit+1+(-N:N));
end
% note that these will be symmetrized later as follows:
% initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

%% Parameters for the proof %%

% weights in the norms for the Banach space
nu=1.0001;
problem.proof.nu=nu;

% an a priori bound on the validation radius
problem.proof.rstar=2e-1;

% computational constants used in computing exponential integrals:
% the number of quadrature nodes used for computations of (exponential) integral
problem.proof.integrals.K1=200; 
% the number of Chebyshev interpolation nodes used to represent the output of the map
problem.proof.interpolation.K0=100; 


%% PDE definition %%

% symbolic variable
syms u;

% the (even) order of the PDE (2J in the latex)
problem.pde.order=4;

if exist('intval.m','file')
  % determine maximum error in rounding, relative to eps
  % for all monomial coefficients
  alllambda=[lambda1,lambda2];
  epsfactor=max(sup(intval(rad(alllambda))./eps(mid(alllambda))));
  % turn relevant coefficient(s) into floats
  lambda2=mid(lambda2);
  lambda1=mid(lambda1);
end

% nonlinearity{j} represents the term g^{(j-1)}(u)
nonlinearity{2} = -lambda1*u^2;
nonlinearity{3} = -lambda2*u;

%% Preprocessing %%
[problem.pde.polynomials,symmetry]=preprocessnonlinearities(nonlinearity,epsfactor);
problem = preprocessdata(problem,initialdata,K,NQ,griddata,symmetry,scale);

%% Final settings %%
problem.useintervals = [];

%% Plot data %%
plotdata=[];

end