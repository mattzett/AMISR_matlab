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
% 3 February 2010
%
%-------------------------------------------------------------------------%

clear all; clc;

config_main01b_2008;

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
range_idx = find( Alt>150e3 & Alt<250e3 ); % Find just the points within these altitudes.
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
allsky_file = fullfile(allsky_dir,sprintf('UnwarpedImages.%04d(OLD).mat',wavelength));
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
                vest_geo(x,y,:) = Rgmag' * squeeze(vest(x,y,:));
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
        if plot_temp,
            display_Ti;
        end;
        
        %------------------------%
        % DISPLAY VELOCITY FIELD %
        %------------------------%
        display_vfield_noTi;
        
        set(gca,'FontSize',6);
        
        hold off;

%         axis equal
        xlabel('Ground distance east (km)','FontSize',6);
        ylabel('Ground distance north (km)','FontSize',6);
        TITLE = sprintf('Radar: %s--%s UT \n Allsky: %s UT',datestr(mtime(1,tr)),datestr(mtime(2,tr),'HH:MM:SS'),datestr(itime(ti),'HH:MM:SS'));
        title(TITLE,'FontSize',6);
        
        if output,
            set(1,'PaperUnits','centimeters');
            papersize = 16.9/1.75*[1.0, 0.85];
            set(1,'PaperSize',papersize);
            set(1,'PaperPosition',[0 0 papersize]);
            outfilename = fullfile(outputdir,datestr(itime(ti),'HHMMSS'));
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
    % Hmm... never got around to doing this.
    
%     figure(6);
    
end; % tr (radar time indices loop)



% %% Display the history as a 1-D plot
% 
% figure(2); % vx
% plot_idx = 1;
% for y = Ny:-1:1,
%     for x = 1:Nx,
%         subplot(Ny,Nx,plot_idx),
%         errorbar( mtime(1,tstart:tend), squeeze(VX(y,x,tstart:tend)), squeeze(DVX(y,x,tstart:tend)) );
%         datetick('x','HHMM')
%         if y == 1,
%             xlabel('UT');
%         end;
%         if x == 1,
%             ylabel('v_x');
%         end;
%         title(sprintf('(%d, %d)', x, y));
% %         axis([mtime(1,tstart),mtime(1,tend),-2000,500]);
%         axis([mtime(1,tstart),mtime(1,tend),-2000,1500]);
%         
%         hline(0,'k-');
%         
%         plot_idx = plot_idx+1;
%     end;
% end;
% set(gcf,'PaperSize',[14,10]);
% set(gcf,'PaperPosition',[0,0,14,10]);
% 
% figure(3); % vy
% plot_idx = 1;
% for y = Ny:-1:1,
%     for x = 1:Nx,
%         subplot(Ny,Nx,plot_idx),
%         errorbar( mtime(1,tstart:tend), squeeze(VY(y,x,tstart:tend)), squeeze(DVY(y,x,tstart:tend)) );
%         datetick('x','HHMM')
%         if y == 1,
%             xlabel('UT');
%         end;
%         if x == 1,
%             ylabel('v_y');
%         end;
%         title(sprintf('(%d, %d)', x, y));
% %         axis([mtime(1,tstart),mtime(1,tend),-500,1500]);
%         axis([mtime(1,tstart),mtime(1,tend),-1500,500]);
%         
%         hline(0,'k-');
%         
%         plot_idx = plot_idx+1;
%     end;
% end;
% set(gcf,'PaperSize',[14,10]);
% set(gcf,'PaperPosition',[0,0,14,10]);
% 
% figure(4); % vz
% plot_idx = 1;
% for y = Ny:-1:1,
%     for x = 1:Nx,
%         subplot(Ny,Nx,plot_idx),
%         errorbar( mtime(1,tstart:tend), squeeze(VZ(y,x,tstart:tend)), squeeze(DVZ(y,x,tstart:tend)) );
%         datetick('x','HHMM')
%         if y == 1,
%             xlabel('UT');
%         end;
%         if x == 1,
%             ylabel('v_z');
%         end;
%         title(sprintf('(%d, %d)', x, y));
%         axis([mtime(1,tstart),mtime(1,tend),-20,20]);
% %         axis([mtime(1,tstart),mtime(1,tend),-200,200]);
%         
%         hline(0,'k-');
%         
%         plot_idx = plot_idx + 1;
%     end;
% end;
% set(gcf,'PaperSize',[14,10]);
% set(gcf,'PaperPosition',[0,0,14,10]);
% 
% print('-f2','-depsc2','-r200',fullfile(outputdir,'history_vx.eps'));
% print('-f3','-depsc2','-r200',fullfile(outputdir,'history_vy.eps'));
% print('-f4','-depsc2','-r200',fullfile(outputdir,'history_vz.eps'));

%% Plot MSP data along with vectors versus time
if isempty(fig_id),
    
    %% Get data from NetCDF file
    MSP_File = '/shared/classes/projects/msp/MSP_2008086.PF';
    % Open the NetCDF file
    ncid = netcdf.open(MSP_File,'NC_NOWRITE');
    % Get time, peak, and background
    time = netcdf.getVar(ncid,0);
    peak = netcdf.getVar(ncid,1);
    back = netcdf.getVar(ncid,2);
    d_start = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Start_Day_UT');
    % Close the NetCDF file.
    netcdf.close(ncid);

    %% Change time to MATLAB format
    time = time';
    hour = double(time) / 3600.0;
    DOY = str2double(d_start(5:7));
    Year = str2double(d_start(9:13));
    Day_MATLAB = DOY + datenum([Year,0,0,0,0,0]);
    [year,month,day] = datevec(Day_MATLAB);
    msptime = Day_MATLAB + double(time)/(3600*24);
    
    %% Eliminate the background
    diff = double(peak(:,1:4,:) - back(:,1:4,:));

    %% Plot MSP data

    TITLES = {'5577','4278','4861','6300?'};
    
    T1 = find(msptime<=StartTime,1,'last');
    T2 = find(msptime>=EndTime,1,'first');

%     % Plot the MSP data.
%     figure(2);
%     plot_idx = 1;
%     ElRange = 1:90;
%     imagesc(msptime(T1:T2),...
%             ElRange,...
%             squeeze(diff(ElRange,plot_idx,T1:T2)),...
%             [0,10000]);
%     datetick('x','HHMM','keeplimits');
%     colormap(jet(256)); colorbar;
%     xlabel('Time (UT)');
%     ylabel('Elevation (deg)');
%     title(TITLES{plot_idx});
%     
%     % Plot vectors
%     
%     hold on;
%     
%     X = (mtime(1,tstart:tend)'*[1 1  1 1])';
%     NewEl = atan2(120e3,...
%                   sqrt(mean(XC(1:4,2:3),2).^2 ...
%                      + mean(YC(1:4,2:3),2).^2)...
%                  ) * 180/pi;
%     Y = (ones(tend-tstart+1,1)*NewEl')';
% %     Y = (ones(tend-tstart+1,1)*[80,90,100])';
%     U = squeeze(mean(VX(1:4,2:3,tstart:tend),2)) * 3e-6;
%     V = squeeze(mean(VY(1:4,2:3,tstart:tend),2)) * -1e-2;
%     quiver(X,Y,U,V,0,'.','color','w');
%     
%     hold off;
    
    Vmag = sqrt(mean(VX(1:4,2:3,tstart:tend),2).^2 + ...
                mean(VY(1:4,2:3,tstart:tend),2).^2 );
    Vdir = atan2(-mean(VY(1:4,2:3,tstart:tend),2).^2,...
                 mean(VX(1:4,2:3,tstart:tend),2).^2 ) * 180 / pi;

    figure(2);
    lambda_idx = 1;
    NewEl = atan2(120e3,...
                  sqrt(mean(XC(1:4,2:3),2).^2 ...
                     + mean(YC(1:4,2:3),2).^2)...
                 ) * 180/pi;
    for El_idx = 4:-1:1,
        
        subplot(4,1,El_idx);
        
          El = round(NewEl(El_idx));
          
          imagesc(msptime(T1:T2),...
                  linspace(0,1000,5),...
                  squeeze(diff(El-2:El+2,lambda_idx,T1:T2)),...
                  [0,10000]);
          colormap(gray(256));
          datetick('x','HHMM','keeplimits');
          axis xy
          
          hold on;
          
          [AX,Hmag,Hdir] = ...
              plotyy(mtime(1,tstart:tend), squeeze(Vmag(El_idx,:,:)),...
                     mtime(1,tstart:tend), squeeze(Vdir(El_idx,:,:)),...
                     'plot','plot');
          hold off;
          
          set(AX(1),'FontSize',6)
          set(AX(2),'FontSize',6)
          datetick(AX(1),'x','HHMM','keeplimits');
          datetick(AX(2),'x','HHMM','keeplimits');
          set(AX(1),'YLim',[0,1000],'YTickMode','auto');
          set(get(AX(1),'YLabel'),'String','Speed (m/s)','FontSize',8)
          set(AX(2),'YLim',[-180,90],...
                    'YTick',[-180,-90,0,90],...
                    'YTickLabel',{'W','S','E','N'});
          set(get(AX(2),'YLabel'),'String','Direction','FontSize',8);
          set(Hdir,'LineStyle','--');
          xlabel('Time (UT)');
          
          title(sprintf('El: %d deg',El));
    end;
    %% Export to image file.
    papersize = [8.5,11];
    set(2,'PaperUnits','centimeters');
    set(2,'PaperSize',papersize);
    set(2,'PaperPosition',[0,0,papersize]);
    print(2,'-dpng','-r600',fullfile(outputdir,'msp_vi'));
    
    
    %% Line plot of brightness and speed.
    
    figure(3);
    lambda_idx = 1;
    El_idx = 3;
    El = round(NewEl(El_idx));
    
    [AX,HL,Hv] = ...
        plotyy(msptime(T1:T2),       squeeze(mean(diff(El,lambda_idx,T1:T2),1))/1000,...
               mtime(1,tstart:tend), squeeze(Vmag(El_idx,:,:)));
        axis(AX,'tight');
        
    set(AX(2),'XTick',[]);
    datetick(AX(1),'x','HHMM','keeplimits')
    set(AX(1),'FontSize',10);
    set(AX(2),'FontSize',10);
    
    xlabel(AX(1),'Time (UT)');
    ylabel(AX(1),'MSP brightness (counts \times1000)');
    ylabel(AX(2),'Speed (m/s)');
    
    set(AX(1),'YLim',[0,20],'YTickMode','auto');
    set(AX(2),'YLim',[-1000,1500],...
              'YTick',[-1000,0,500,1000],...
              'YTickLabel',{'','0','500','1000'});
    
    set(Hv,'Marker','s');
    
    papersize = [8.5, 5];
    set(3,'PaperUnits','inches');
    set(3,'Papersize',papersize);
    set(3,'PaperPosition',[0 0 papersize]);
    print('-f3','-dpdf','-r600',fullfile(outputdir,'msp_vi2.pdf'));
end