% This produces all figures (the numbering of the figures corresponds to
% the one used in the paper)

% directory to get the data from
datadir='data'; 
% directory to put the figures in
figuredir='results';
[~,~]=mkdir(figuredir);

%%%%%%%%%%%%
% Figure 1 %
%%%%%%%%%%%%

problemname='OhtaKawasaki';
datafilename=[datadir,'/datafloat',problemname];
load(datafilename,'x','problem');

plotdata.viewangle = [10 60];
plotdata.xfactor = 0.5;
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,plotdata,figurefilename);
plotdata.viewangle = [];
plotdata.xfactor = [];

%%%%%%%%%%%%
% Figure 3 %
%%%%%%%%%%%% 

griddatafilename=[figuredir,'/example',problemname,'_griddata_'];
plotgriddata(problem,griddatafilename);

%%%%%%%%%%%%
% Figure 2 %
%%%%%%%%%%%% 

problemname='SwiftHohenberg';
datafilename=[datadir,'/datafloat',problemname];
load(datafilename,'x','problem');
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,[],figurefilename);

plotdata.viewangle = [-160 30];
figurefilename=[figuredir,'/example',problemname,'_bis'];
plotsolution(x,problem,plotdata,figurefilename);
plotdata.viewangle = [];

%%%%%%%%%%%%
% Figure 4 %
%%%%%%%%%%%% 

problemname='KuramotoSivashinsky';
datafilename=[datadir,'/datafloat',problemname];
load(datafilename,'x','problem');

plotdata.viewangle = [-15 40];
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,plotdata,figurefilename);
plotdata.viewangle = [];


