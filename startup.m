if ismac
        homedir = '/Users/';
else
        homedir = '/home/';
end

%opengl software
addpath([homedir 'lansdell/matlab/correlationpointprocess/']);

addpath([homedir 'lansdell/projects/connectivity/functions']);
addpath_recurse([homedir 'lansdell/projects/connectivity/preprocess']);
addpath_recurse([homedir 'lansdell/projects/connectivity/eval']);
addpath_recurse([homedir 'lansdell/projects/connectivity/glm/models']);
addpath_recurse([homedir 'lansdell/projects/connectivity/glm/fitting']);
addpath_recurse([homedir 'lansdell/projects/connectivity/ldglm/']);

