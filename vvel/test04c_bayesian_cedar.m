%-------------------------------------------------------------------------%
% Test 4c: Bayesian velocity estimator. Figures for CEDAR poster
% 
% Thomas Butler
% 26 May 2009
%
%-------------------------------------------------------------------------%

clear all; clc;

config_test04_cedar;

%% Gather data from HDF5 file
bco   = hdf5read(h5file,'BeamCodes');
fits  = hdf5read(h5file,'FittedParams/Fits');
errs  = hdf5read(h5file,'FittedParams/Errors');
utime = hdf5read(h5file,'Time/UnixTime');
Alt   = hdf5read(h5file,'Geomag/Altitude');
Range = hdf5read(h5file,'Geomag/Range');

%% Format the data as we need it...

az = bco(2,:) * pi/180;
el = bco(3,:) * pi/180;
Nr = length(az);

vlos_all  = zeros(size(fits,3), size(fits,4), size(fits,5));
Ti_all    = zeros(size(fits,4), size(fits,5));
dvlos_all = zeros(size(vlos_all));

RANGE_GATE = zeros(size(az));
rr = zeros(size(az));
for beam = 1:length(az),
    % Range gate corresponding to HEIGHT
    RANGE_GATE(beam) = find(abs(Alt(:,beam) - HEIGHT) == min(abs(Alt(:,beam) - HEIGHT)));
    rr(beam) = Range(RANGE_GATE(beam),beam); % Actual range (km) for each radar beam corresponding to approximately HEIGHT altitude
%     rr(beam) = Range(5,beam); % Actual range (km) for each radar beam corresponding to approximately HEIGHT altitude
    vlos_all(:,beam,:)  = fits(4,3,:               ,beam,:); % Fit 4 = vlos, Species 3 = electrons
    Ti_all(beam,:)      = fits(2,1,RANGE_GATE(beam),beam,:); % Fit 2 = Temp, Species 1 = hydrogen (same temp for all ions)
    dvlos_all(:,beam,:) = errs(4,1,:               ,beam,:); % Species 1, fitter sets all vlos equal, so only one error
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

clear bco fits errs utime;

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
allsky_dir = ['/shared/classes/projects/AMISR-experiments/data/',DateStr,'_Allsky/UnwarpedImages'];
eval(['load ',allsky_dir]);

%----------------------------------------%
% PREPARE TO GATHER ESTIMATES AND ERRORS %
%----------------------------------------%
VX = zeros(Ny,Nx,tend-tstart+1);
VY = VX;
VZ = VX;

DVX = VX; DVY = VY; DVZ = VZ;

figure(1);
set(gcf,'Color','w');
set(gcf,'InvertHardCopy','off');

%% Loop over the radar (H5) time indices
for tr = tstart:tend
%% Grab the current ISR parameters for time tr
    vlos = vlos_all(:,:,tr); % Just the current time
    vlos = vlos(range_idx);  % Just the range gates represented by xgmag, ygmag, zgmag
    Ti   = Ti_all(:,tr);
    
    dvlos = dvlos_all(:,:,tr);
    dvlos = dvlos(range_idx);

%% Velocity reconstruction

    x2 = linspace(min(xgmag), max(xgmag), Nx+2);
    y2 = linspace(min(ygmag), max(ygmag), Ny+2);
    [X2,Y2] = meshgrid(x2,y2);
    if Bayesian,
        vest = zeros(Ny,Nx,3);
    else
        vest = zeros(Ny,Nx,2);
    end;

    for x = 1:Nx,
        for y = 1:Ny,

            % Grab the measurements corresponding to this "pixel"
            vlos_idx = find(xgmag >= x2(x) & xgmag <= x2(x+2) & ygmag >= y2(y) & ygmag <= y2(y+2));
            vlos_current = vlos(vlos_idx);
            dvlos_current = dvlos(vlos_idx);

            if Bayesian,
                A2 = direction_vectors(vlos_idx,1:3);
                Pe = diag(dvlos_current.^2);
                Pv = diag([500, 500, 15].^2);
               
                vest(y,x,:) = Pv*A2'*inv(A2*Pv*A2' + Pe) * vlos_current;
                
                Pvest = Pv - Pv*A2'*inv(A2*Pv*A2' + Pe)*A2*Pv;
                
                VX(y,x,tr) = vest(y,x,1);
                VY(y,x,tr) = vest(y,x,2);
                VZ(y,x,tr) = vest(y,x,3);
                
                DVX(y,x,tr) = sqrt(Pvest(1,1));
                DVY(y,x,tr) = sqrt(Pvest(2,2));
                DVZ(y,x,tr) = sqrt(Pvest(3,3));
            else
                A2 = direction_vectors(vlos_idx,1:2);
                
%                 vest(y,x,:) = A2\vlos_current;
                vest(y,x,:) = pinv(A2)*vlos_current;
                
                VX(y,x,tr) = vest(y,x,1);
                VY(y,x,tr) = vest(y,x,2);
                
                Pe = diag(dvlos_current.^2);
                Pvest = inv(A2' / Pe * A2);             %    ???
                
                DVX(y,x,tr) = sqrt(Pvest(1,1));
                DVY(y,x,tr) = sqrt(Pvest(2,2));
            end;
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
    
    Ti_I = griddata(rr .* cos(el) .* sin(az), ...
                    rr .* cos(el) .* cos(az), ...
                    Ti, ...
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

        %-----------------%
        % DISPLAY OPTICAL %
        %-----------------%

        image_disp = duw;
        image_disp = repmat((image_disp-gmin)/(gmax-gmin),[1,1,3]);
        image_disp = image_disp.^(1 - beta); % Brighten the image
        imshow(image_disp,'XData',xp/1e3,'YData',yp/1e3);
        set(gca,'Color',[0,0,0]);


        colormap(gray(256));
        
        axis xy
        axis on
        
        hold on;
        
        %---------------------%
        % DISPLAY TEMPERATURE %
        %---------------------%
        a = 15;
        contour(xp/1e3, yp/1e3,...
                filter2(1/a^2*ones(a),Ti_I,'same'),...
                logspace(log10(500),log10(3000),10),...
                'LineWidth',2);
        caxis([500,3000]);
        colormap(hot(20));
        cbar_handle = colorbar('Location','East');
        set(get(cbar_handle,'ylabel'),'String','Ion Temperature (K)','Color','w');

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
        xlabel('Ground distance east (km)','FontSize',12);
        ylabel('Ground distance north (km)','FontSize',12);
        TITLE = sprintf('Radar: %s--%s UT / Allsky: %s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'),datestr(itime(ti),'HH:MM:SS'));
        title(TITLE,'FontSize',14);
        
        if output,
            outfilename = sprintf('%s/vest_%03d',outputdir,ti);
            print('-f1','-depsc2','-r100',outfilename);
            print('-f1','-dpng','-r100',outfilename);
        end;

    end; %ti (image time indices loop)

%%
end; % tr (radar time indices loop)

%% Display the history as a 1-D plot
figure(2);
for y = 1:Ny,
    for x = 1:Nx,

        figure(2);
        subplot(Ny,Nx,(y-1)*Ny+x),
        errorbar( mtime(1,tstart:tend), squeeze(VX(y,x,:)), squeeze(DVX(y,x,:)), '.' );
        datetick('x','HHMM')
        if y == Ny,
            xlabel('UT');
        end;
        if x == 1,
            ylabel('v_x');
        end;
        title(sprintf('(%d, %d)', x, y));
        axis([mtime(1,tstart),mtime(1,tend),-2000,500]);
    end;
end;
set(gcf,'PaperSize',[14,10]);
set(gcf,'PaperPosition',[0,0,14,10]);

figure(3); % vy
for y = 1:Ny,
    for x = 1:Nx,
        subplot(Ny,Nx,(y-1)*Ny+x),
        errorbar( mtime(1,tstart:tend), squeeze(VY(y,x,:)), squeeze(DVY(y,x,:)), '.' );
        datetick('x','HHMM')
        if y == Ny,
            xlabel('UT');
        end;
        if x == 1,
            ylabel('v_y');
        end;
        title(sprintf('(%d, %d)', x, y));
        axis([mtime(1,tstart),mtime(1,tend),-500,1500]);
    end;
end;
set(gcf,'PaperSize',[14,10]);
set(gcf,'PaperPosition',[0,0,14,10]);

if Bayesian,
figure(4); % vz
for y = 1:Ny,
    for x = 1:Nx,
        subplot(Ny,Nx,(y-1)*Ny+x),
        errorbar( mtime(1,tstart:tend), squeeze(VZ(y,x,:)), squeeze(DVZ(y,x,:)), '.' );
        datetick('x','HHMM')
        if y == Ny,
            xlabel('UT');
        end;
        if x == 1,
            ylabel('v_z');
        end;
        title(sprintf('(%d, %d)', x, y));
        axis([mtime(1,tstart),mtime(1,tend),-100,100]);
        
    end;
end;
set(gcf,'PaperSize',[14,10]);
set(gcf,'PaperPosition',[0,0,14,10]);
end;

if Bayesian,
    print('-f2','-depsc2','-r200',fullfile(outputdir,'history_vx.eps'));
    print('-f3','-depsc2','-r200',fullfile(outputdir,'history_vy.eps'));
    print('-f4','-depsc2','-r200',fullfile(outputdir,'history_vz.eps'));
else
    print('-f2','-depsc2','-r200',fullfile(outputdir,'history_vx_nobayesian.eps'));
    print('-f3','-depsc2','-r200',fullfile(outputdir,'history_vy_nobayesian.eps'));
end;