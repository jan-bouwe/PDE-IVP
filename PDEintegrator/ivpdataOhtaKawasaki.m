function [problem,plotdata] = ivpdataOhtaKawasaki
% data for Ohta-Kawasaki equation
%
% u_t = - 1/gamma^2 u_{xxxx} - (u-u^3)_{xx} - sigma(u-m)
%
% for x in [0,L] and t in [0,tau]
%
% returns the nonlinearities, 
% the problem definition, including the initial data 
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
% u_t = - 1/gamma^2 u_{xxxx} - (u-u^3)_{xx} - sigma(u-m)
% for x in [0,L] and t in [0,tau]

gamma=sqrt(8*altone); 
sigma=altone/5;
m=altone/10;
L=4*altpi;
tau=30*altone;

% rescale to spatial domain [0,2*pi] and
% u_t = -u_{xxxx} + lambda21*(u-u^3)_{xx} + lambda01*(u-m)

% scale.space=L/(2*altpi);
% scale.time=gamma^2*scale.space^4;
scale.space=2*altone;
scale.time=8*scale.space^4;
integrationtime=tau/scale.time;

% lambda21=-scale.time/scale.space^2;
lambda21=-8*scale.space^2;
lambda23=-lambda21;
lambda01=-sigma*scale.time;
lambda00=-lambda01*m;

%% Truncation and grid %%

% number of Fourier modes (2N+1 really)
N=60; 
% number of Chebyshev modes (K+1 really)
K=4; 
% number of time domains (called M in the latex)
D=80; 

% size of matrices in the semigroup (called N_L in the latex)
NQ=30;

% total integration time
griddata.T=integrationtime;
griddata.D=D;

 % we let the code pick the subdivision
griddata.variable=true;
% and set some algorithmic parameters for finding this subdivision
griddata.accuracy=0.1;
griddata.Dvar=1;
griddata.maxsteps=10;

% we want to integrate all the way to infinity
griddata.infinite=true; 

% optimize the number K of chebyshev modes and K0 used for approximating
% via interpolation, in each subdomain
griddata.optimizecheb=true;
% and set some algorithmic parameters for this
griddata.Ybound1target=1e-8;
griddata.Ybound2target=1e-10;
griddata.Kmin=1;
griddata.Kmax=50;
griddata.K0min=30;
griddata.K0max=600;
griddata.K0step=20;

%% Initial data %%

% set initial data 
initialdata=altzeros([2*N+1,1],m);
initialdata(N+1)=m;
initialdata(N+3)=2*altone/10;
initialdata(N+5)=2*altone/10;

% note that these will be symmetrized later as follows (note factor 2):
% initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

%% Parameters for the proof %%

% weights in the norms for the Banach space
nu=1.0001;
problem.proof.nu=nu;

% an a priori bound on the validation radius
problem.proof.rstar=2e-1;

% computational constants used in computing exponential integrals:
% the number of quadrature nodes used for computations of (exponential) integral
problem.proof.integrals.K1=300; 
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
  alllambda=[lambda21,lambda23,lambda01,lambda00];
  epsfactor=max(sup(intval(rad(alllambda))./eps(mid(alllambda))));
  % turn relevant coefficient(s) into floats
  lambda21=mid(lambda21);
  lambda01=mid(lambda01);
  m=mid(m);
end

% nonlinearity{j} represents the term g^{(j-1)}(u)
nonlinearity{1} = lambda01*(u-m);
nonlinearity{3} = lambda21*(u-u^3);

%% Preprocessing %%
[problem.pde.polynomials,symmetry]=preprocessnonlinearities(nonlinearity,epsfactor);
problem = preprocessdata(problem,initialdata,K,NQ,griddata,symmetry,scale);

%% Final settings %%
problem.useintervals = [];

%% Plot data %%
plotdata=[];

end