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
% Reconstruction type
% reconstruction_type = 'individual';
% reconstruction_type = 'holistic';
% reconstruction_type = 'holistic_linear'; alpha = 3e-1;
reconstruction_type = 'holistic_incompressible'; alpha = 2;
% reconstruction_type = 'holistic_incompressible_altBC'; alpha = 2;

%---------------------------%
% Data filtering parameters %
%---------------------------%

% % Josh's JGR paper: Figure 10 "growth phase"
% StartTime = datenum(2008,03,26,11,38,00);
% EndTime   = datenum(2008,03,26,11,41,00);
% gmin = 350; gmax = 500; beta = 0.4;
% plot_temp = false;
% fig_id = 'Fig10a_growth';

% % Josh's JGR paper: Figure 10 "expansion phase"
% StartTime = datenum(2008,03,26,11,46,30);
% EndTime   = datenum(2008,03,26,11,49,00);
% gmin = 450; gmax = 2500; beta = 0.5;
% plot_temp = false;
% fig_id = 'Fig10b_expansion';

% % Josh's JGR paper: Figure 10 "recovery phase"
% StartTime = datenum(2008,03,26,11,53,00);
% EndTime   = datenum(2008,03,26,11,55,30);
% gmin = 400; gmax = 600; beta = 0.0;
% plot_temp = false;
% fig_id = 'Fig10c_recovery';

% Josh's JGR paper: Figure 12 "arc element activation"
StartTime = datenum(2008,03,26,11,30,00);
EndTime   = datenum(2008,03,26,11,36,33);
gmin = 350; gmax = 700; beta = 0.6;
cmin = 500; cmax = 3000;
plot_temp = true;
fig_id = 'Fig12';

% % Josh's JGR paper: Figure 13 "substorm recovery"
% StartTime = datenum(2008,03,26,12,12,00);
% EndTime   = datenum(2008,03,26,12,31,00);
% gmin = 330; gmax = 700; beta = 0.0;
% cmin = 500; cmax = 3000;
% plot_temp = true;
% fig_id = 'Fig13';

% % Josh's JGR paper: MSP and velocity time history
% StartTime = datenum(2008,03,26,11,00,00);
% EndTime   = datenum(2008,03,26,12,30,00);
% gmin = 350; gmax = 500; beta = 0.4;
% plot_temp = false;
% fig_id = '';

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
if isempty(fig_id),
    output    = 0; % Output results to image files?
else
    output    = 1;
end;

outputdir = fullfile('main01b_results',fig_id);
if ( output && ~exist(outputdir,'file') ),
    mkdir(outputdir);
end;
