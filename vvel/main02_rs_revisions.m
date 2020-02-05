%-------------------------------------------------------------------------%
% Main 01: Velocity field estimator.  Estimate 3-component vecotr velocity
% from line-of-sight ISR measurements of F-region velocities.
%
% Modified and modularized from test04_bayesian.m
%
% Allow flexibility for future additions:
%   Varying types of reconstruction (whole-A2 vs. individual A2 inversion vs. potential estimation?)
%   For whole-A2 inversion, use different types of regularization (linear variation by bin, adjust kernel by optical intensity, etc.)
%   Bayesian / non-Bayesian
%   ...
% 
% Thomas Butler
% 2 July 2009
%
%-------------------------------------------------------------------------%

clear all; clc;

config_main01_2009;

output=0;
outputdir = ['/shared/classes/2008-2009/Summer2009/rs-paper-velocity-field/figures/20090324_final/',reconstruction_type,'/',Tint];

masochism=false; % Very time- and disk space-consuming plots. Not usually necessary.

%% Gather data from HDF5 file
bco   = hdf5read(h5file,'BeamCodes');
fits  = hdf5read(h5file,'FittedParams/Fits');
errs  = hdf5read(h5file,'FittedParams/Errors');
utime = hdf5read(h5file,'Time/UnixTime');
Alt   = hdf5read(h5file,'Geomag/Altitude');
Range = hdf5read(h5file,'Geomag/Range');

%% Format the data as we need it...

% Convert az & el from degrees to radians.
az = bco(2,:) * pi/180;
el = bco(3,:) * pi/180;
Nr = length(az);

RANGE_GATE = zeros(size(az));
rr = zeros(size(az));
vlos_all  = zeros(size(fits,3), size(fits,4), size(fits,5));
Ti_all    = zeros(size(fits,4), size(fits,5));
dvlos_all = zeros(size(vlos_all));
dTi_all   = zeros(size(Ti_all));

for beam = 1:length(az),
    % Range gate corresponding to HEIGHT
    RANGE_GATE(beam) = find(abs(Alt(:,beam) - HEIGHT) == min(abs(Alt(:,beam) - HEIGHT)));
    
    
    % Actual range (km) for each radar beam corresponding to approximately HEIGHT altitude
    rr(beam) = Range(RANGE_GATE(beam),beam); 
%     % Actual range (km) for each radar beam in the fifth range gate
%     rr(beam) = Range(5,beam);
    
    
    % Gather vlos & Ti (and errors) as a function of range, beam #, and time.
    vlos_all(:,beam,:)  = fits(4,3,:               ,beam,:); % Fit 4 = vlos, Species 3 = electrons
    Ti_all(beam,:)      = fits(2,1,RANGE_GATE(beam),beam,:); % Fit 2 = Temp, Species 1 = hydrogen (same temp for all ions)
    dvlos_all(:,beam,:) = errs(4,1,:               ,beam,:); % Species 1, fitter sets all vlos equal, so only one error
    dTi_all(beam,:)     = errs(2,1,RANGE_GATE(beam),beam,:);

end;


% Compute unit vectors for each beam.
el1 = repmat(el,size(Alt,1),1);            % Replicate az & el, making each the same size as Alt (or Range)
az1 = repmat(az,size(Alt,1),1);
range_idx = find( Alt>200e3 & Alt<350e3 ); % Find just the points within these altitudes.
el2 = el1(range_idx);                      % Extract only those az & el
az2 = az1(range_idx);

% Unit vectors of the beam at each sample point.
kx = cos(el2) .* sin(az2);
ky = cos(el2) .* cos(az2);
kz = sin(el2);

% Cartesian coordinates of each sample point.
xr = Range(range_idx).*kx;
yr = Range(range_idx).*ky;
zr = Range(range_idx).*kz;

% Rotation matrix from geo -> gmag
Rgmag = [cos(dec),         -sin(dec),          0;
         sin(dip)*sin(dec), cos(dec)*sin(dip), cos(dip);
        -cos(dip)*sin(dec),-cos(dec)*cos(dip), sin(dip)];

% Rotate direction vectors from geo -> gmag
direction_vectors = [kx ky kz] * Rgmag';

% Rotate coordinates to geomagnetic
xyzgmag = [xr,yr,zr] * Rgmag';
xgmag = xyzgmag(:,1);
ygmag = xyzgmag(:,2);
zgmag = xyzgmag(:,3);


% Convert Unix time to MATLAB time
mtime = unixtime2matlab(utime,0);


clear bco fits errs utime;


% Determine the main loop limits (indexes into mtime)
tstart = find(mtime(1,:) >= StartTime,1,'first');
tend   = find(mtime(2,:) <= EndTime,  1,'last');

%% Gather allsky data

%-----------------------------&
% LOAD UNWARPED ALLSKY IMAGES %
%-----------------------------&
allsky_dir = ['/shared/classes/projects/AMISR-experiments/data/',DateStr,'_Allsky/'];
allsky_file = fullfile(allsky_dir,sprintf('UnwarpedImages.%04d.mat',wavelength));
load(allsky_file);

%% Plot up-B velocity versus time.

% Grab ranges (i) corresponding to beam 13 (j). Note, just one j index needed.
[i,j] = find(Alt > 200e3 & Alt < 350e3);
i = i(j==13); % Beam 13 = "Up B"
j = j(j==13);
vlos_altrange = squeeze(vlos_all(i,j(1),:));
dvlos_altrange = squeeze(dvlos_all(i,j(1),:));
r_altrange = squeeze(Range(i,j(1)));

W = std(vlos_altrange,[],2).^-1;
W = W / sum(W)  * size(vlos_altrange,1);
W = repmat(W(:),[1,size(vlos_altrange,2)]);
v_upB = mean(vlos_altrange)';
v_upB_w = mean(W.*vlos_altrange)';

dv_upB = mean(dvlos_altrange)';

figure(1);
subplot(2,1,1);

% plot(mtime(1,:),vlos_altrange(:,:)')
hold on;
plot(mtime(1,:),v_upB,  'r','LineWidth',3)
plot(mtime(1,:),v_upB_w,'k','LineWidth',3)
hold off;
datetick('x','HH:MM','keeplimits')
ylim([-1000,1000])
xlabel('Time (UT)');
ylabel('v_{los} ("up B")');
title('Red: Unweighted average  Black: Weighted average')

%% Display all-sky images with a up-B position marked;
if masochism
    
x = lims(1):lims(2);
y = lims(3):lims(4);

x_upB = Range(i,13) * cos(el(13)) * sin(az(13)) / 1e3;
y_upB = Range(i,13) * cos(el(13)) * cos(az(13)) / 1e3;

allsky_dir = '/shared/classes/projects/AMISR-experiments/data/20090324_Allsky/FITS_20090324/';
Az = fitsread('/shared/CLASSES/projects/AMISR-experiments/data/Allsky_Starfield_Fits/PKR_DASC_AZ_Nov2007.FITS');
El = fitsread('/shared/CLASSES/projects/AMISR-experiments/data/Allsky_Starfield_Fits/PKR_DASC_EL_Nov2007.FITS');
% Az = Az(1:2:end,1:2:end); El = El(1:2:end,1:2:end);
Az = resample(resample(Az,1,2)',1,2)';
El = resample(resample(El,1,2)',1,2)';

az360 = 180/pi*az;
el360 = 180/pi*el;
az360(az360<0) = az360(az360<0) + 360;

for beam = 1:26,
    [el_as(beam),az_as(beam)] = find(          abs(Az - az360(beam)) + abs(El - el360(beam)) ==       ... 
                                      min(min( abs(Az - az360(beam)) + abs(El - el360(beam)) )) , 1 );
end;

figure(1);
for t_idx = 1:length(itime),
    v_time = find(abs(mtime(1,:) - itime(t_idx)) == min(abs(mtime(1,:) - itime(t_idx))),1);
    
    subplot(2,1,1);
    hold on;
    if t_idx > 1,
        delete(h_marker);
    end;
    h_marker = vline(mtime(1,v_time),'b--');
    set(h_marker,'LineWidth',2);
    
    subplot(2,2,3);
    imagesc(x,y,duw_all(:,:,t_idx));
    axis xy;
    colormap(gray(256));
    title(sprintf('%s  v_{los}=%6.1f m/s (%6.1f m/s, weighted)',datestr(itime(t_idx),'HH:MM:SS'),v_upB(v_time),v_upB_w(v_time)));
    hold on;
    plot(x_upB,y_upB,'rx');
    hold off;

    subplot(2,2,4);
    allsky_file = fullfile(allsky_dir,['PKR_DASC_0558_090324_',datestr(itime(t_idx),'HHMMSS'),'.000.FITS']);
    allsky_img = fitsread(allsky_file);
    imagesc(allsky_img,[0,1000]);
    colormap(gray(256)); colorbar;
    title(sprintf('%s  v_{los}=%6.1f m/s (%6.1f m/s, weighted)',datestr(itime(t_idx),'HH:MM:SS'),v_upB(v_time),v_upB_w(v_time)));
    hold on;
    plot(az_as,el_as,'rx');
    hold off;
    axis xy;
    axis square;
    
    print('-f1','-dpng','-r100',sprintf('main02_rs_revisions/vlos_upB/%03d',t_idx));
    
end;

end; % masochism

%% Plot histograms of up-B velocity versus range

figure(2);
for r_idx = 1:length(r_altrange),
    subplot(length(r_altrange),1,r_idx);
    hist(vlos_altrange(r_idx,:), -1000:50:1000);
    text(-900,10,sprintf('%3.0f km',r_altrange(r_idx)/1e3));
    xlim([-1000,1000])
end;