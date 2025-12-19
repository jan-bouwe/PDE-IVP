% this is an initialization file
% please adapt the path appropriately

 intlabpath = [];

% If you have intlab provide the path to it below.
% Comment out the next line if you don't have intlab

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

