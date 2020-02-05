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

%% Prepare containers for velocity estimates and errors.

%----------------------------------------%5
% PREPARE TO GATHER ESTIMATES AND ERRORS %
%----------------------------------------%
VX = zeros(Ny,Nx,tend-tstart+1);
VY = VX;
VZ = VX;

DVX = VX; % i.e. zeros
DVY = VY;
DVZ = VZ;

%% Loop over the radar (H5) time indices
figure(1);

for tr = tstart:tend
    %% Grab the current ISR parameters for time tr
    vlos  = vlos_all(:,:,tr); % Just the current time
    Ti    = Ti_all(:,tr);
    
    dvlos = dvlos_all(:,:,tr);
    dTi   = dTi_all(:,tr);

    vlos  = vlos(range_idx);  % Just the range gates represented by xgmag, ygmag, zgmag
    dvlos = dvlos(range_idx);

    %% Velocity reconstruction

    switch lower(reconstruction_type)
        case 'individual'
            [vest,Pvest,x2,y2] = vfield_individual(vlos,dvlos,xgmag,ygmag,direction_vectors,Nx,Ny);
        case 'holistic'
            [vest,Pvest,x2,y2] = vfield_holistic(vlos,dvlos,xgmag,ygmag,direction_vectors,Nx,Ny,alpha);
        case 'holistic_linear'
            [vest,Pvest,x2,y2] = vfield_holistic_linear(vlos,dvlos,xgmag,ygmag,direction_vectors,Nx,Ny,alpha);
        case 'holistic_incompressible'
            [vest,Pvest,x2,y2] = vfield_holistic_incompressible(vlos,dvlos,xgmag,ygmag,direction_vectors,Nx,Ny,alpha);
        case 'holistic_incompressible_altbc'
            [vest,Pvest,x2,y2] = vfield_holistic_incompressible_altBC(vlos,dvlos,xgmag,ygmag,direction_vectors,Nx,Ny,alpha);
        otherwise
            error('reconstruction_type not recognized');
    end;
    
    % Rotate estimate to geographic coordinates
    rotate = true;
    vest_geo = zeros(size(vest));
    if rotate;
        for x = 1:Nx,
            for y = 1:Ny,
                vest_geo(y,x,:) = Rgmag' * squeeze(vest(y,x,:));
            end;
        end;
        vest = vest_geo;
    end;
    
    VX(:,:,tr) = vest(:,:,1);
    VY(:,:,tr) = vest(:,:,2);
    VZ(:,:,tr) = vest(:,:,3);

    DVX(:,:,tr) = sqrt(Pvest(:,:,1,1));
    DVY(:,:,tr) = sqrt(Pvest(:,:,2,2));
    DVZ(:,:,tr) = sqrt(Pvest(:,:,3,3));
    
    %% Interpolate Ti onto X & Y

    xp = (lims(1):lims(2))*1e3;
    yp = (lims(3):lims(4))*1e3;
    
    Ti_I = griddata(rr .* cos(el) .* sin(az), ...
                    rr .* cos(el) .* cos(az), ...
                    Ti, ...
                    xp',yp);
    
    %% Load the unwarped allsky image
    this_itime = find( (itime >= mtime(1,tr)) & (itime < mtime(2,tr)) );

%     figure(1);

    for ti = this_itime,
%         image_file = image_file_list(ti).name;

        z0 = 120; % kilometers

%         lims = [-100,200,-50,200];

        duw = duw_all(:,:,ti);

        %-----------------%
        % DISPLAY OPTICAL %
        %-----------------%
        display_optical;
        
        hold on;
        
        %---------------------%
        % DISPLAY TEMPERATURE %
        %---------------------%
%         display_Ti;
        
        %------------------------%
        % DISPLAY VELOCITY FIELD %
        %------------------------%
%         plot(xr/1e3,yr/1e3,'bx');
        
%         display_dvlos_overlay;
        
        display_vfield_noTi;
        
        set(gca,'FontSize',6);
        
        hold off;

%         axis equal
        xlabel('Ground distance east (km)','FontSize',6);
        ylabel('Ground distance north (km)','FontSize',6);
        TITLE = sprintf('%s--%s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'));
        title(TITLE,'FontSize',6);
        text(120,290,datestr(itime(ti),'HH:MM:SS'),'FontSize',6,'Color','w');
        
        if output,
            set(1,'PaperUnits','centimeters');
            papersize = 16.9/3.1*[1.0,0.95];
            set(1,'PaperSize',papersize);
            set(1,'PaperPosition',[0 0 papersize]);
            outfilename = sprintf('%s/vest_%03d',outputdir,ti);
            print('-f1','-dpng','-r300',outfilename);
        end;

    end; %ti (image time indices loop)

%     %% Display velocity error ellipses
% 
%     figure(5);
%     clf;
%     display_error_ellipses;
%     
%     annotation('textbox',[.25 .9 .8 .1],...
%                'String',sprintf('Velocity error, %s--%s',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS') ),...
%                'FontSize',12,...
%                'LineStyle','None');
% 
%     outfilename = sprintf('%s/error_ellipses_t%03d',outputdir,tr);
%     print('-f5','-depsc2','-r200','-painters',outfilename);
% %     print('-f5','-dpng','-r600',outfilename);
    
    %% Display temperature variances as an image
    
%     figure(6);
    
end; % tr (radar time indices loop)



%% Display the history as a 1-D plot

figure(2); % vx
plot_idx = 1;
for y = Ny:-1:1,
    for x = 1:Nx,
        subplot(Ny,Nx,plot_idx),
        errorbar( mtime(1,tstart:tend), squeeze(VX(y,x,tstart:tend)), squeeze(DVX(y,x,tstart:tend)), 'b-' );
        datetick('x','HHMM')
        if y == 1,
            xlabel('UT');
        end;
        if x == 1,
            ylabel('v_x');
        end;
        title(sprintf('(%d, %d)', x, y));
%         axis([mtime(1,tstart),mtime(1,tend),-2000,500]);
        axis([mtime(1,tstart),mtime(1,tend),-2500,1100]);
        
        hline(0,'k-');
        
        plot_idx = plot_idx+1;
    end;
end;
set(gcf,'PaperSize',[14,10]);
set(gcf,'PaperPosition',[0,0,14,10]);

figure(3); % vy
plot_idx = 1;
for y = Ny:-1:1,
    for x = 1:Nx,
        subplot(Ny,Nx,plot_idx),
        errorbar( mtime(1,tstart:tend), squeeze(VY(y,x,tstart:tend)), squeeze(DVY(y,x,tstart:tend)), 'b-' );
        datetick('x','HHMM')
        if y == 1,
            xlabel('UT');
        end;
        if x == 1,
            ylabel('v_y');
        end;
        title(sprintf('(%d, %d)', x, y));
        axis([mtime(1,tstart),mtime(1,tend),-1000,1500]);
%         axis([mtime(1,tstart),mtime(1,tend),-1500,500]);
        
        hline(0,'k-');
        
        plot_idx = plot_idx+1;
    end;
end;
set(gcf,'PaperSize',[14,10]);
set(gcf,'PaperPosition',[0,0,14,10]);

figure(4); % vz
plot_idx = 1;
for y = Ny:-1:1,
    for x = 1:Nx,
        subplot(Ny,Nx,plot_idx),
        errorbar( mtime(1,tstart:tend), squeeze(VZ(y,x,tstart:tend)), squeeze(DVZ(y,x,tstart:tend)), 'b-' );
        datetick('x','HHMM')
        if y == 1,
            xlabel('UT');
        end;
        if x == 1,
            ylabel('v_z');
        end;
        title(sprintf('(%d, %d)', x, y));
        axis([mtime(1,tstart),mtime(1,tend),-100,100]);
%         axis([mtime(1,tstart),mtime(1,tend),-200,200]);
        
        hline(0,'k-');
        
        plot_idx = plot_idx + 1;
    end;
end;
set(gcf,'PaperSize',[14,10]);
set(gcf,'PaperPosition',[0,0,14,10]);

print('-f2','-depsc2','-r200',fullfile(outputdir,'history_vx.eps'));
print('-f3','-depsc2','-r200',fullfile(outputdir,'history_vy.eps'));
print('-f4','-depsc2','-r200',fullfile(outputdir,'history_vz.eps'));
