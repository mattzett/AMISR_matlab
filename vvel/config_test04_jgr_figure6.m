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
% StartTime = datenum(2008,03,26,11,26,00); % ONSET
% EndTime   = datenum(2008,03,26,12,34,00);
StartTime = datenum(2008,03,26,12,00,00); % TRANSLATING BOUNDARY
EndTime   = datenum(2008,03,26,13,00,00);
gmin = 350; gmax = 4000; beta = 0.5;

DateStr   = datestr(StartTime,'yyyymmdd');

HEIGHT = 180e3;

%----------------%
% HDF5 data file %
%----------------%
% Tint     = '30sec'; % Integration time
Tint     = '2min';
% Tint     = '5min';
h5path   = '/shared/classes/projects/AMISR-experiments/data/';
h5subdir = [DateStr,'.001/'];
h5path   = fullfile(h5path,h5subdir);
h5file   = [DateStr '.001_lp_' Tint '.h5'];
h5file   = fullfile(h5path,h5file);

%-------------------%
% Allsky data files %
%-------------------%
wavelength = 0;
image_path = ['/shared/classes/projects/AMISR-experiments/data/',DateStr,'_Allsky/FITS_',DateStr,'/'];
image_file_list = dir([image_path 'PKR_DASC_' sprintf('%04d',wavelength) '*.FITS']);

%------%
% Misc %
%------%
% Scaling for quiver plots
qscaling = 1.0e-2;
qfac = 0.03;

% Plotting options
output    = 1; % Output results to image files?

outputdir = 'test04_results_2008';
if ( output && ~exist(outputdir,'dir') ),
    mkdir(outputdir);
end;
