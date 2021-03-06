%-------------------------------------------------------------------------%
% Test 4b: Bayesian velocity estimator.  Plot just v and Ti (the latter as
% a filled contour plot) for Josh's JGR(?) paper.
% 
% Thomas Butler
% 16 June 2009
%
%-------------------------------------------------------------------------%

clear all; clc;

config_test04_2008;
outputdir = 'test04b_results_2008';
if ( output && ~exist(outputdir,'file') ),
    mkdir(outputdir);
end;

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

%----------------------------------------%
% PREPARE TO GATHER ESTIMATES AND ERRORS %
%----------------------------------------%
VX = zeros(Ny,Nx,tend-tstart+1);
VY = VX;
VZ = VX;

DVX = VX; DVY = VY; DVZ = VZ;

figure(1);
set(gcf,'Color','w');

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
                Pvest = inv(A2' / Pe * A2); % ???
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
        
        %---------------------%
        % DISPLAY TEMPERATURE %
        %---------------------%
        a = 15;
        contourf(xp/1e3, yp/1e3,...
                 filter2(1/a^2*ones(a),Ti_I,'same'),...
                 'LineWidth',1);
        caxis([0,3000]);
        colormap(jet);
        cbar_handle = colorbar('East');
%         set(get(cbar_handle,'xlabel'),'String','Temperature (K)','Color','w');

        axis xy;
        hold on;
        
        set(gca,'Color',[0,0,0]);

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
               'Color','w',...
               'LineWidth',2);
%         xlim([xp(1),xp(end)]/1e3);
%         ylim([yp(1),yp(end)]/1e3);
 
        
        quiver(-90, -40, 1e3*qfac, 0, 0, 'Color', 'w', 'LineWidth', 2);
        text(-90,-30,'1 km/s','color','w')
        
        hold off;

%         axis equal
        xlabel('E-W distance (km)','FontSize',14);
        ylabel('N-S distance (km)','FontSize',14);
        TITLE = sprintf('%s--%s UT',datestr(mtime(1,tr),'HH:MM:SS'),datestr(mtime(2,tr),'HH:MM:SS'));
        text(-90,195,TITLE,'Color','w','FontSize',18);
        text(125,185,'T_i (K)','Color','w','FontSize',20);
        set(gcf,'InvertHardCopy','off');
        
        if output,
            outfilename = sprintf('%s/vest_%s',outputdir,datestr(mtime(1,tr),'HHMMSS'));
            print('-f1','-dpng','-r100',outfilename);
        end;

end; % tr (radar time indices loop)
