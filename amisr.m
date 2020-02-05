%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ONE PROGRAM TO RULE THEM ALL %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script reads in the amisr hdf5 files that are posted on SRIs website:
% http://amisr.com/amisr/links/data-access/
% is doesn't work directly for the madrigal data without some minor modification


clear; close all;
addpath ./geom ./vvel ./plot_imgfiles;


% PICK THE DATA 
%filename='20170222.002_lp_2min-fitcal_ISINGLASS_A.h5';  % ISINGLASS A, 02/22/2017, 10:14:00-10:23:55 UT (subpayload), 10:14:00-10:24:52 UT (main)
% filename='20170302.002_lp_2min-fitcal_ISINGLASS_B.h5';  % ISINGLASS B, 03/02/2017, ???
 filename='20130207.001_lp_1min-cal_VISIONS.h5';         % VISIONS, 2/7/2013
%filename='20131207.003_lp_3min-cal.h5';                 % Scintilation comparison
%filename='20131208.001_lp_1min-cal.h5';                 % Scintilation comparison
%filename='20150318.001_lp_1min-cal.h5';                   %St. Patrick's Day storm


% NEED WINDOWS MOVIE OPTION?
opflag=0; % 1 - Windows, 0 - Not Windows


% CONVERTS FROM HDF5 FORMAT TO MAT FORMAT
saveflag=1;  % 1 = make mat file, 0 = no mat file  (make the .mat files for now)
filelab=hdf5_extract(filename,saveflag);  % produces _rawdata.mat
fprintf(['AMISR --> File selected: ',filelab,'\n'])


% SETUP IMAGE STORAGE LOCATION
fprintf('AMISR --> Setting up storage locations\n');
outdir=['./plot_imgfiles/',filelab];
if ~exist(outdir)
    mkdir(outdir);
else
    delete([outdir,'/*.png'])
end


% SELECT THE TIME STEPS TO PLOT
load(['./datasets/data_mat/',filelab,'_rawdata.mat']);  % loads up data to determine total time frames and relevant UT times
itstart   = 1;  % start time frame
itfin     = size(isne,3); % end time frame
% itstart = 100; % can pick a specific time frame window if desired
% itfin   = 130; 
starting=exp_date(itstart,4)+exp_date(itstart,5)/60+exp_date(itstart,6)/3600;
fining=exp_date(itfin,4)+exp_date(itfin,5)/60+exp_date(itfin,6)/3600;
fprintf('AMISR --> Selected start and stop times\n');
fprintf(['AMISR --> ',num2str(starting),' UT to ',num2str(fining),' UT\n']);


% PLOTS BEAMS AND MAG FIELD 
fprintf('AMISR --> Plotting beams and magnetic field lines\n');
payloadflag=3;    % 0 = no rocket trajectory, 1 = ISINGLASS A, 2 = ISINGLASS B, 3 = VISIONS
plot_grid(filelab,payloadflag)


% CALCULATES THE BEAM ENU COORDINATES
fprintf('AMISR --> Prepping fieldgrid data\n');
grid_ISR(filelab)


% PLOTS 3D PFISR DATA (AND MAKES A MOVIE)
fprintf('AMISR --> Plotting altitudes of PFISR data\n');
makemovie_AMISR_vlos(filelab,itstart,itfin,opflag)


% ESTIMATES FLOW FIELDS (BUTLER 2010 METHOD)
fprintf('AMISR --> Calculating velocities\n');
estimate_flowfield(filelab,itstart,itfin)


% PLOTS THE VELOCITY FLOWS (AND MAKES A MOVIE)
fprintf('AMISR --> Quiver plots\n');
plot_ISRflows(filelab,itstart,itfin,opflag)

fprintf('All Done!\n')

rmpath  ./geom ./vvel ./plot_imgfiles;



