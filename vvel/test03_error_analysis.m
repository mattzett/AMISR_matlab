%-------------------------------------------------------------------------%
% Test 3: Propagate error through estimation process of test02.m.  Data from
% 24 March 2009.
% 
% Thomas Butler
% 18 May 2009
%
% NOTE: (26 May 2009) Still not complete. Only temperature error
% information has been collected.
%-------------------------------------------------------------------------%

clear all; clc;

config_test03;

%% Gather data from HDF5 file
bco   = hdf5read(h5file,'BeamCodes');
fits  = hdf5read(h5file,'FittedParams/Fits');
err   = hdf5read(h5file,'FittedParams/Errors');
utime = hdf5read(h5file,'Time/UnixTime');
Alt   = hdf5read(h5file,'Geomag/Altitude');
Range = hdf5read(h5file,'Geomag/Range');

%% Format the data as we need it...

az = bco(2,:) * pi/180;
el = bco(3,:) * pi/180;
Nr = length(az);

vlos_all = zeros(size(fits,3), size(fits,4), size(fits,5));
Ti_all   = zeros(size(fits,4), size(fits,5));
dvlos_all = zeros(size(vlos_all));
dTi_all   = zeros(size(Ti_all));

RANGE_GATE = zeros(1,length(az));
for beam = 1:length(az),
    % Range gate corresponding to HEIGHT
    RANGE_GATE(beam) = find(abs(Alt(:,beam) - HEIGHT) == min(abs(Alt(:,beam) - HEIGHT)));
    vlos_all(:,beam,:)  = fits(4,3,:               ,beam,:); % Fit 4 = vlos, Species 3 = electrons
    Ti_all(beam,:)      = fits(2,1,RANGE_GATE(beam),beam,:); % Fit 2 = Temp, Species 1 = hydrogen (same temp for all ions)
    dvlos_all(:,beam,:) =  err(4,1,:               ,beam,:); % Species 1, fitter sets all vlos equal, only one error
    dTi_all(beam,:)     =  err(2,1,RANGE_GATE(beam),beam,:);
end;

el1 = repmat(el,size(Alt,1),1);
az1 = repmat(az,size(Alt,1),1);
range_idx = find( Alt>150e3 & Alt<245e3 );
el2 = el1(range_idx);
az2 = az1(range_idx);
kx = cos(el2) .* sin(az2);
ky = cos(el2) .* cos(az2);
kz = sin(el2);
% direction_vectors = [kx ky kz];

xr = Range(range_idx).*kx;
yr = Range(range_idx).*ky;
zr = Range(range_idx).*kz;

Rgmag = [cos(dec),         -sin(dec),          0;
         sin(dip)*sin(dec), cos(dec)*sin(dip), cos(dip);
        -cos(dip)*sin(dec),-cos(dec)*cos(dip), sin(dip)];
direction_vectors = [kx ky kz] * Rgmag';
xyzgmag = [xr,yr,zr] * Rgmag';
xgmag = xyzgmag(:,1);
ygmag = xyzgmag(:,2);
zgmag = xyzgmag(:,3);

mtime = unixtime2matlab(utime,0);

clear bco fits err utime;

% Determine the main loop limits (indexes into mtime)
tstart = find(mtime(1,:) >= StartTime,1,'first');
tend   = find(mtime(2,:) <= EndTime,  1,'last');

%-------------------------------%
% TIME VECTOR FOR ALLSKY IMAGES %
%-------------------------------%
ScanStr = ['PKR_DASC_' sprintf('%04d', wavelength) '_%2d%2d%2d_%2d%2d%2d.000.FITS'];
itime = zeros(1,length(image_file_list));
for t = 1:length(image_file_list),
    A = sscanf( image_file_list(t).name, ScanStr );
    A(1) = A(1) + 2000;
    itime(t) = datenum(A');
end;

%-----------------------------&
% LOAD UNWARPED ALLSKY IMAGES %
%-----------------------------&
load ~butler/classes/projects/AMISR-experiments/data/20090324_Allsky/UnwarpedImages

figure(1);

%% Loop over the radar (H5) time indices
for tr = tstart:tend
%% Grab the current ISR parameters for time tr
    vlos = vlos_all(:,:,tr); % Just the current time
    vlos = vlos(range_idx);  % Just the range gates represented by xgmag, ygmag, zgmag
    Ti   = Ti_all(:,tr);
    dTi  = dTi_all(:,tr);

%% Velocity reconstruction

    x2 = linspace(min(xgmag), max(xgmag), Nx+2);
    y2 = linspace(min(ygmag), max(ygmag), Ny+2);
    [X2,Y2] = meshgrid(x2,y2);
    vest = zeros(Ny,Nx,2);

    for x = 1:Nx,
        for y = 1:Ny,

            % Grab the measurements corresponding to this "pixel"
            vlos_idx = find(xgmag >= x2(x) & xgmag <= x2(x+2) & ygmag >= y2(y) & ygmag <= y2(y+2));
            vlos_current = vlos(vlos_idx);

            A2 = direction_vectors(vlos_idx,1:2);

%             vest(y,x,:) = A2\vlos_current;
            vest(y,x,:) = pinv(A2)*vlos_current;

        end; % for y
    end; % for x
    
%     vest = zeros([size(VERTICES2,1),2]);
%     vest(:,1) = vest_array(1:2:end);     % x-component
%     vest(:,2) = vest_array(2:2:end);     % y-component
% 
    
%% Interpolate Ti onto X & Y

    lims = [-100,200,-50,200]*1e3; % NB: This was set in test06a, so you can't change it here
    xp = lims(1):1e3:lims(2);
    yp = lims(3):1e3:lims(4);
    
    Ti_I = griddata(HEIGHT * cos(el) .* sin(az), ...
                    HEIGHT * cos(el) .* cos(az), ...
                    Ti, ...
                    xp',yp);
    dTi_I= griddata(HEIGHT * cos(el) .* sin(az), ...
                    HEIGHT * cos(el) .* cos(az), ...
                    dTi, ...
                    xp',yp);
    
%% Load the unwarped allsky image
    this_itime = find( (itime >= mtime(1,tr)) & (itime <= mtime(2,tr)) );
    for ti = this_itime,
        image_file = image_file_list(ti).name;

        z0 = 120; % kilometers

%         lims = [-100,200,-50,200];

        duw = duw_all(:,:,ti);

%         xp = lims(1):1:lims(2);
%         yp = lims(3):1:lims(4);
% 
%         xr = Range .* cos(el) .* sin(az);
%         yr = Range .* cos(el) .* cos(az);
%         zr = Range .* sin(el);


        figure(1);

        %-----------------%
        % DISPLAY OPTICAL %
        %-----------------%
%         gmin = min(min(min(images(0.75*xmin:1.2*xmax,0.75*ymin:1.2*ymax,:))));
%         gmax = max(max(max(images(0.75*xmin:1.2*xmax,0.75*ymin:1.2*ymax,:))));

        gmin = 350;
        gmax = 1200;
        
        image_disp = repmat((duw-gmin)/(gmax-gmin),[1,1,3]);
        imshow(xp/1e3,yp/1e3,image_disp);

%         imagesc(xp,yp,duw,[350,600]);
%         colormap(gray(256));
        axis xy
        axis on
        
        hold on;
        
        %---------------------%
        % DISPLAY TEMPERATURE %
        %---------------------%
        a = 15;
        contourf(xp/1e3, yp/1e3,...
                filter2(1/a^2*ones(a),Ti_I,'same'),...
                linspace(500,2000,10),...
                'LineWidth',1);
        caxis([500,2000]);
        colormap(hot(10));
        cbar_handle = colorbar;
        set(get(cbar_handle,'ylabel'),'String','Temperature (K)');

        %------------------------%
        % DISPLAY VELOCITY FIELD %
        %------------------------%

        % Show the reconstructed velocity field
        xc = x2(2:end-1)';
        yc = y2(2:end-1)';
        [XC,YC] = meshgrid(xc,yc);
        temp = [XC(:), YC(:), 240e3*ones(numel(XC),1)] * Rgmag;
        xc = temp(:,1);
        yc = temp(:,2);
        vx = vest(:,:,1);
        vy = vest(:,:,2);
        quiver(xc/1e3,     yc/1e3, ...
               vx(:)*qfac, vy(:)*qfac, ...
               0,...
               'Color','c',...
               'LineWidth',2);
%         xlim([xp(1),xp(end)]/1e3);
%         ylim([yp(1),yp(end)]/1e3);
 
        
        quiver(-90, -40, 1e3*qfac, 0, 0, 'Color', 'c', 'LineWidth', 2);
        text(-90,-30,'1 km/s','color','c')

        hold off;

%         axis equal
        xlabel('Ground distance east (km)');
        ylabel('Ground distance north (km)');
        TITLE = sprintf('Radar: %s--%s UT / Allsky: %s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'),datestr(itime(ti),'HH:MM:SS'));
        title(TITLE);
        
        figure(2);

        %--------------------%
        % DISPLAY TEMP ERROR %
        %--------------------%
        a = 15;
        contourf(xp/1e3, yp/1e3,...
                filter2(1/a^2*ones(a),dTi_I,'same'),...
                linspace(0,1000,20),...
                'LineStyle','none');
        caxis([0,1000]);
        colormap(jet(20));
        cbar_handle = colorbar;
        set(get(cbar_handle,'ylabel'),'String','\Delta Temperature (K)');
               
        xlabel('Ground distance east (km)');
        ylabel('Ground distance north (km)');
        TITLE = sprintf('Radar: %s--%s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'));
        title(TITLE);
        
        if output,
            outfilename = sprintf('%s/Ti_err_%03d',outputdir,tr);
            print('-f2','-dpng','-r100',outfilename);
        end;

        figure(3)
        
        %-------------------------------%
        % DISPLAY FRACTIONAL TEMP ERROR %
        %-------------------------------%
        a = 15;
        contourf(xp/1e3, yp/1e3,...
                filter2(1/a^2*ones(a),dTi_I./Ti_I,'same'),...
                linspace(0,1,20),...
                'LineStyle','none');
        caxis([0,1]);
        colormap(jet(20));
        cbar_handle = colorbar;
        set(get(cbar_handle,'ylabel'),'String','\DeltaT_i / T_i');
               
        xlabel('Ground distance east (km)');
        ylabel('Ground distance north (km)');
        TITLE = sprintf('Radar: %s--%s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'));
        title(TITLE);
        
        if output,
            outfilename = sprintf('%s_2/Ti_err_%03d',outputdir,tr);
            print('-f3','-dpng','-r100',outfilename);
        end;

    end; %ti (image time indices loop)

%%
end; % tr (radar time indices loop)

