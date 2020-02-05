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
Nx = 1; Ny = 4;
% Reconstruction type
% reconstruction_type = 'individual';
% reconstruction_type = 'holistic'; alpha = 1;
% reconstruction_type = 'holistic_linear'; alpha = 4e-1;
reconstruction_type = 'holistic_incompressible'; alpha = 0.001;
% reconstruction_type = 'holistic_incompressible_altBC'; alpha = 5;

%---------------------------%
% Data filtering parameters %
%---------------------------%
% StartTime = datenum(2009,03,24,07,49,57);
% EndTime   = datenum(2009,03,24,08,36,01);

% StartTime = datenum(2009,03,24,07,49,57);
% EndTime   = datenum(2009,03,24,09,00,03);

StartTime = datenum(2009,03,24,07,49,57);
% EndTime   = datenum(2009,03,24,10,00,03);
EndTime = datenum(2009,03,24,08,40,00);    % <-- Focusing on the period with suspicious flow reversals on the far-most bin (lowest elevation)
gmin = 420; gmax = 1800; beta = 0.5;
cmin = 500; cmax = 3000;

% StartTime = datenum(2009,03,24,08,40,00);
% EndTime   = datenum(2009,03,24,10,00,03);
% gmin = 400; gmax = 600; beta = 0.0;
% cmin = 500; cmax = 3000;

DateStr   = datestr(StartTime,'yyyymmdd');

HEIGHT = 240e3;

%----------------%
% HDF5 data file %
%----------------%
Tint     = '30sec'; % Integration time
% Tint     = '2min';
% Tint     = '5min';
h5path   = '/shared/classes/projects/AMISR-experiments/data/';
h5subdir = [DateStr,'.001/'];
h5path   = fullfile(h5path,h5subdir);
h5file   = [DateStr '.001_lp_' Tint '.h5'];
h5file   = fullfile(h5path,h5file);

%-------------------%
% Allsky data files %
%-------------------%
% wavelength = 0;
wavelength = 558;  % nm
% wavelength = 630;  % nm
image_path = ['/shared/classes/projects/AMISR-experiments/data/',DateStr,'_Allsky/FITS_',DateStr,'/'];
image_file_list = dir([image_path 'PKR_DASC_' sprintf('%04d',wavelength) '*.FITS']);

%------%
% Misc %
%------%
% Scaling for quiver plots
qscaling = 1.0e-2;
qfac = 0.06;

% Plotting options
output    = 1; % Output results to image files?

outputdir = fullfile('main01_results',DateStr,sprintf('%dx%d',Nx,Ny),Tint,reconstruction_type,num2str(alpha));
if ( output && ~exist(outputdir,'file') ),
    mkdir(outputdir);
end;
