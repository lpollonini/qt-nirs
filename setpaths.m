function setpaths(option)

%
% USAGE: 
%
%   setpaths
%   setpaths(1)
%   setpaths(0)
%

if ~exist('option','var')
    option = 1;
end
setpathsPath = which('qtnirs.m');
rootdir = fileparts(setpathsPath);
dir2include = genpath(rootdir);
dirs2include = tokenize(dir2include,';');
for i=1:length(dirs2include)
   disp(dirs2include{i}); 
end

p = { 
    'Snirf';
    'Snirf/Examples';
    'Utils';
    'Utils/Hdf5';
    };
    
for ii = 1:length(p)
    pname = [rootdir, '/', p{ii}];
    if option
        addpath(pname, '-end');
    else
        rmpath(pname);
    end
end


