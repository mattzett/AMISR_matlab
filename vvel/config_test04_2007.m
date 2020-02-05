% Config file for velocity reconstruction through projection along magnetic
% field. (test01)

%------------%
% PARAMETERS %
%------------%

% Approximate declination and dip angles at PFISR location
dec = 22   * pi/180;
dip = 77.5 * pi/180;

%------------------------%
% Reconstruction options %
%------------------------%
% Resolution
Nx = 4; Ny = 4;
% Use Bayesian estimator?
Bayesian = 1;

%---------------------------%
% Data filtering parameters %
%---------------------------%
StartTime = datenum(2007,11,10,09,14,50);
EndTime   = datenum(2007,11,10,09,26,00);
gmin = 310; gmax = 500; beta = 0.1;

DateStr   = datestr(StartTime,'yyyymmdd');

HEIGHT = 240e3;

%----------------%
% HDF5 data file %
%----------------%
% Tint     = '30sec'; % Integration time
% Tint     = '1min';
% Tint     = '2min';
% Tint     = '5min.MJN';
Tint     = '5min.TWB';
h5path   = '/shared/classes/projects/AMISR-experiments/AMISR_fitter_py.20090721/runfiles/20071110_lp/';
% h5subdir = [DateStr,'.001/'];
% h5path   = fullfile(h5path,h5subdir);
h5file   = [DateStr '.001_lp_' Tint '.h5'];
h5file   = fullfile(h5path,h5file);

%-------------------%
% Allsky data files %
%-------------------%
wavelength = 0;
image_path = ['/shared/classes/projects/AMISR-experiments/data/',DateStr,'_Allsky/FITS_',DateStr,'/'];
image_file_list = dir([image_path 'AS' sprintf('%04d',wavelength) '*.FITS']);

%------%
% Misc %
%------%
% Scaling for quiver plots
qscaling = 1.0e-2;
qfac = 0.03;

% Plotting options
output    = 1; % Output results to image files?

outputdir = 'test04_results_2007';
if ( output && ~exist(outputdir,'file') ),
    mkdir(outputdir);
end;
