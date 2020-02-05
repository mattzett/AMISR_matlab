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
StartTime = datenum(2008,03,26,00,00,00);
EndTime   = datenum(2008,03,27,00,00,00);
gmin = 350; gmax = 4000; beta = 0.5;

% % Josh's JGR paper: Figure 9 "at the expansion phase onset"
% StartTime = datenum(2008,03,26,11,33,00);
% EndTime   = datenum(2008,03,26,12,00,00);
% gmin = 330; gmax = 2500; beta = 0.7;

% % Josh's JGR paper: Figure 10 "during a period 2-min prior to onset"
% StartTime = datenum(2008,03,26,11,30,00);
% EndTime   = datenum(2008,03,26,11,40,00);
% gmin = 330; gmax = 800; beta = 0.6;

% % Josh's JGR paper: Figure 11 "an approaching discrete auroral boundary"
% StartTime = datenum(2008,03,26,12,12,00);
% EndTime   = datenum(2008,03,26,12,31,00);
% gmin = 330; gmax = 700; beta = 0.0;

% Josh's JGR paper: use only with test04b_bayesian_just_v_and_Ti.m
% StartTime = datenum(2008,03,26,11,40,00);
% EndTime   = datenum(2008,03,26,16,00,00);

DateStr   = datestr(StartTime,'yyyymmdd');

HEIGHT = 240e3;

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
% wavelength = 558;
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
output    = 0; % Output results to image files?

outputdir = 'test04_results_2008';
if ( output && ~exist(outputdir,'file') ),
    mkdir(outputdir);
end;
