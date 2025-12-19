function plotgriddata(problem,filename)
% plots the subdomain sizes and the number of Chebyshev nodes per subdomain
% (the numbers used for the solution and for the interpolation)
%
% if filename is provided and nonempty then the plots are saved

plotssaved = (exist('filename','var') && ~isempty(filename));

% figure subdivision
resc = problem.scalefloats.time;
tau = resc * altmid(problem.timegrid(1:end-1));
deltatau = resc * altmid(problem.domains);
fig = figure;
plot(tau, deltatau, '.k', 'MarkerSize', 20)
xlabel('time')
ylabel('length of subdomain')
set(gca,'FontSize',15) 

% save plot
if plotssaved
    saveplot(fig,[filename,'domains']);
end

% figures K 
fig = figure;
plot(tau, problem.sol.K, '.k', 'MarkerSize', 20)
xlabel('time')
ylabel('$K$','Interpreter','latex')
set(gca,'FontSize',15) 

% save plot
if plotssaved
    saveplot(fig,[filename,'K']);
end

% figure K0 (tildeK in the latex)
fig = figure;
plot(tau, problem.proof.interpolation.K0, '.k', 'MarkerSize', 20)
xlabel('time')
ylabel('$\tilde{K}$','Interpreter','latex')
set(gca,'FontSize',15) 

% save plot
if plotssaved
    saveplot(fig,[filename,'K0']);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveplot(fig,filename)
% saves a pdf of figure fig in file filename

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,filename,'-dpdf','-r300')

end



