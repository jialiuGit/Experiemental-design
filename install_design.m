

%Run install_design.m in the root directory of Experimental_design
% If Matlab is started in GPstuff root directory, this startup.m
% adds subfolders automatically to the path
F = mfilename;
S = which(F);
subfolders={'design'  'gp_design'};

for sf=subfolders
  addpath(strrep(S,[F '.m'],sf{:}))
end


% Alternatively copy following lines in your main MATLAB script

%designroot='/...'
%addpath([design 'design'])
%addpath([design 'gp_design'])
 