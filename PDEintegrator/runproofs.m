% This file runs all proofs of the theorems

initialize

% directory to get the data from
datadir='data'; 
% directory to put the results
resultdir='results';
[~,~]=mkdir(resultdir);

if ~exist('intval.m','file')
    % run all "proofs" with floats
    useintervals=false;
    datatype='float';
else
    % run all proofs rigorously
    useintervals=true;
    datatype='intval';
end

%%%%%%%%%%%%%%%
% Theorem 1.2 %
%%%%%%%%%%%%%%%

problemname='OhtaKawasaki';
datafilename=[datadir,'/data',datatype,problemname];
load(datafilename,'x','problem');

% plot of the solution
plotsolution(x,problem);
% NB: restricting back to the PDE domain 
% (computational domain was 2L because of Neumann BC)
L = problem.scalefloats.space*pi;
xlim([0 L]); 
title('Ohta-Kawasaki')
set(gca,'FontSize',20)
drawnow

% proof
fprintf("\nRunning the proof for the Ohta-Kawasaki equation (Theorem 1.2)\n\n") 
tic
[success,rsol,bounds,errorbound] = ivpproof(x,problem,useintervals);
toc
tictoc=toc;
resultsfilename=[resultdir,'/dataproof',datatype,problemname];
save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','useintervals','tictoc');

%%%%%%%%%%%%%%%
% Theorem 6.2 %
%%%%%%%%%%%%%%%

problemname='SwiftHohenberg';
datafilename=[datadir,'/data',datatype,problemname];
load(datafilename,'x','problem');

% plot of the solution
plotsolution(x,problem);
title('Swift-Hohenberg')
set(gca,'FontSize',20)
drawnow

% proof
fprintf("\n\nRunning the proof for the Swift-Hohenberg equation (Theorem 6.2)\n\n")
tic
[success,rsol,bounds,errorbound] = ivpproof(x,problem,useintervals);
toc
tictoc=toc;
resultsfilename=[resultdir,'/dataproof',datatype,problemname];
save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','useintervals','tictoc');

%%%%%%%%%%%%%%%
% Theorem 6.3 %
%%%%%%%%%%%%%%%

problemname='KuramotoSivashinsky';
datafilename=[datadir,'/data',datatype,problemname];
load(datafilename,'x','problem');

% plot of the solution
plotsolution(x,problem);
title('Kuramoto-Sivashinsky')
set(gca,'FontSize',20)
drawnow

% proof
fprintf("\n\nRunning the proof for the Kuramoto-Sivashinsky equation (Theorem 6.3)\n\n")
tic
[success,rsol,bounds,errorbound] = ivpproof(x,problem,useintervals);
toc
tictoc=toc;
resultsfilename=[resultdir,'/dataproof',datatype,problemname];
save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','useintervals','tictoc');

