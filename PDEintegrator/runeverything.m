% This produces data, theorems and figures for the Ohta-Kawasaki equation
% and for the Kuramoto-Sivashinsky equation. For Swift-Hohenberg there is
% no data to generate, since the data are taken from the previous paper for
% proper comparison.

initialize
datadir='data';
resultdir='results';
[~,~]=mkdir(resultdir);

fprintf("\nGenerating data and running the proof for the Ohta-Kawasaki equation (Theorem 1.2)\n\n")
% This computes the data and runs the proof of Theorem 1.2 and Figures 1 and 3
runivp('OhtaKawasaki',resultdir);
% Because of the homogeneous Neumann BC, we doubled the size of the domain
% and used even symmetry, so the solution of the PDE is half of what is
% shown here

fprintf("\n\nGenerating data and running the proof for the Kuramoto-Sivashinsky equation (Theorem 6.3)\n\n")
% This computes the data and runs the proof of Theorem 6.3 and Figure 4
runivp('KuramotoSivashinsky',resultdir);

fprintf("\n\nRunning the proof for the Swift-Hohenberg equation (Theorem 6.2)\n\n")
% This computes the data and runs the proof of Theorem 6.2 and Figure 2
runivpdata('SwiftHohenberg',datadir,resultdir);

% generate the plots for the paper
fprintf("\nGenerating the figures in the paper\n\n")
runplots

