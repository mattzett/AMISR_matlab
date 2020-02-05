%-------------------------------------------------------------------------%
% Test 1: Generate and display estimates of velocity field from
% line-of-sight projections by first rotating measurements into geomagnetic
% coordinates, thengrouping measurements into a 2D grid, inverting the
% measurements to get 2D velocity vectors at each grid point, then rotating
% back to the original geodetic coordinates.
% 
% Thomas Butler
% 1 May 2009
%
%-------------------------------------------------------------------------%

clear all; clc;

config_test01;

%% Gather data from HDF5 file
bco   = hdf5read(h5file,'BeamCodes');
fits  = hdf5read(h5file,'FittedParams/Fits');
utime = hdf5read(h5file,'Time/UnixTime');
Alt   = hdf5read(h5file,'Geomag/Altitude');
Range = hdf5read(h5file,'Geomag/Range');

%% Format the data as we need it...

az = bco(2,:) * pi/180;
el = bco(3,:) * pi/180;
Nr = length(az);

vlos_all = zeros(size(fits,3), size(fits,4), size(fits,5));
Ti_all   = zeros(size(fits,4), size(fits,5));

RANGE_GATE = zeros(1,length(az));
for beam = 1:length(az),
    % Range gate corresponding to HEIGHT
    RANGE_GATE(beam)   = find(abs(Alt(:,beam) - HEIGHT) == min(abs(Alt(:,beam) - HEIGHT)));
    vlos_all(:,beam,:) = fits(4,3,:               ,beam,:); % Fit 4 = vlos, Species 3 = electrons
    Ti_all(beam,:)     = fits(2,1,RANGE_GATE(beam),beam,:); % Fit 2 = Temp, Species 1 = hydrogen (same temp for all ions)
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

clear bco fits utime;

% Determine the main loop limits (indexes into mtime)
tstart = find(mtime(1,:) >= StartTime,1,'first');
tend   = find(mtime(2,:) <= EndTime,  1,'last');

%-------------------------------%
% TIME VECTOR FOR ALLSKY IMAGES %
%-------------------------------%
% ScanStr = ['PKR_DASC_' sprintf('%04d', wavelength) '_%2d%2d%2d_%2d%2d%2d.000.FITS'];
% itime = zeros(1,length(image_file_list));
% for t = 1:length(image_file_list),
%     A = sscanf( image_file_list(t).name, ScanStr );
%     A(1) = A(1) + 2000;
%     itime(t) = datenum(A');
% end;

%-----------------------------&
% LOAD UNWARPED ALLSKY IMAGES %
%-----------------------------&
% load ~butler/classes/projects/AMISR-experiments/data/20090324_Allsky/UnwarpedImages20090324

figure(1);

%% Loop over the radar (H5) time indices
for tr = tstart:tend
%% Grab the current ISR parameters for time tr
    vlos = vlos_all(:,:,tr); % Just the current time
    vlos = vlos(range_idx);  % Just the range gates represented by xgmag, ygmag, zgmag
    Ti   = Ti_all(:,tr);

%% Velocity reconstruction
%     [VERTICES2,X2,Y2,Z2] = reconstruction_geometry(xr,yr,zr,reconstruction_type);
    

    x2 = linspace(min(xgmag), max(xgmag), Nx);
    y2 = linspace(min(ygmag), max(ygmag), Ny);
    [X2,Y2] = meshgrid(x2,y2);
    vest = zeros(size(X2,1)-1,size(X2,2)-1,2);
    
    for x = 1:length(x2)-1,
        for y = 1:length(y2)-1,
            
            % Grab the measurements corresponding to this "pixel"
            vlos_idx = find(xgmag >= x2(x) & xgmag <= x2(x+1) & ygmag >= y2(y) & ygmag <= y2(y+1));
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
    
%% Interpolate Ti onto xp & yp (allsky image coordinates)

    lims = [-100,200,-50,200]*1e3; % NB: This was set in test06a, so you can't change it here
    xp = lims(1):1e3:lims(2);
    yp = lims(3):1e3:lims(4);
    
%     Ti_I = griddata(xr,yr,Ti,xp',yp);
    
% % Load the unwarped allsky image
%     this_itime = find( (itime >= mtime(1,tr)) & (itime <= mtime(2,tr)) );
%     for ti = this_itime,
%         image_file = image_file_list(ti).name;
% 
%         z0 = 120; % kilometers
% 
% %         lims = [-100,200,-50,200];
% 
%         duw = duw_all(:,:,ti);
% 
%         %-----------------%
%         % DISPLAY OPTICAL %
%         %-----------------%
% %         gmin = min(min(min(images(0.75*xmin:1.2*xmax,0.75*ymin:1.2*ymax,:))));
% %         gmax = max(max(max(images(0.75*xmin:1.2*xmax,0.75*ymin:1.2*ymax,:))));
% 
%         gmin = 350;
%         gmax = 600;
%         
%         image_disp = repmat((duw-gmin)/(gmax-gmin),[1,1,3]);
%         imshow(xp/1e3,yp/1e3,image_disp);
% 
% %         imagesc(xp,yp,duw,[350,600]);
% %         colormap(gray(256));
%         axis xy
%         axis on
%         
%         hold on;
        
%         %---------------------%
%         % DISPLAY TEMPERATURE %
%         %---------------------%
%         a = 15;
%         contour(xp/1e3, yp/1e3,...
%                 filter2(1/a^2*ones(a),Ti_I,'same'),...
%                 linspace(500,2000,10),...
%                 'LineWidth',2);
%         caxis([500,2000]);
%         colormap(hot(10));
%         cbar_handle = colorbar;
%         set(get(cbar_handle,'ylabel'),'String','Temperature (K)');
% 
%         hold on;

        %------------------------%
        % DISPLAY VELOCITY FIELD %
        %------------------------%

        % Show the reconstructed velocity field
        xc = ( x2(1:end-1) + x2(2:end) )' / 2;
        yc = ( y2(1:end-1) + y2(2:end) )' / 2;
        [XC,YC] = meshgrid(xc,yc);
        temp = [XC(:), YC(:), 240e3*ones(numel(XC),1)] * Rgmag;
        xc = temp(:,1);
        yc = temp(:,2);
        vx = vest(:,:,1);
        vy = vest(:,:,2);
        quiver(xc/1e3,     yc/1e3,...
               vx(:)*qfac, vy(:)*qfac,...
               0,...
               'Color','k',...
               'LineWidth',2);
        xlim([xp(1),xp(end)]/1e3);
        ylim([yp(1),yp(end)]/1e3);
           
        hold on;
        
        quiver(-90, -40, 1e3*qfac, 0, 0, 'Color', 'k', 'LineWidth', 2);
        text(-90,-30,'1 km/s','color','k')

        hold off;

%         axis equal
        xlabel('Ground distance east (km)');
        ylabel('Ground distance north (km)');
        TITLE = sprintf('%s--%s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'));
        title(TITLE);

    if output,
        outfilename = sprintf('%s/arc_%03d',outputdir,tr);
        print('-f1','-dpng','-r80',outfilename);
    end;

%     end; %ti (image time indices loop)

%%
end; % tr (radar time indices loop)

