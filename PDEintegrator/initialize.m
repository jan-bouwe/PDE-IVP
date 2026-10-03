% this is an initialization file
% please adapt the path appropriately

 intlabpath = [];

% If you have intlab provide the path to it below.
% Uncomment the next line (and adapt the path) if you have intlab

% intlabpath = '~/matlab/intlab/Intlab_V12/';

if ~exist('intval','file') && ~isempty(intlabpath)  
    dir = pwd;
    cd(intlabpath)
    startintlab;
    cd(dir)
    setround(0);
end

clearvars;
if exist('intval','file')
    setround(0);
end

