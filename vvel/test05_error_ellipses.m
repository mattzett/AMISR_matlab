%-------------------------------------------------------------------------%
% Test 5: Plot error ellipses for the experimental setup of March 24, 2009
% and the reconstruction technique used so far (binning in both magnetic
% latitude and longitude with 50% overlap).
%
% Plot error ellipse for each bin.
%
% Also do this for non-overlapping bins?
%
% Thomas Butler
% 22 June 2009
%
%-------------------------------------------------------------------------%

clear all; clc;

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
% Use Bayesian estimator?
Bayesian = 1;

%---------------------------%
% Data filtering parameters %
%---------------------------%
StartTime = datenum(2009,03,24,07,49,57);
EndTime   = datenum(2009,03,24,08,36,01);
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

%----------------%
% Output options %
%----------------%
outputdir = 'test05_results';
if ~exist(outputdir,'file'),
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
dvlos_all = zeros(size(vlos_all));

RANGE_GATE = zeros(size(az));
for beam = 1:length(az),
    % Range gate corresponding to HEIGHT
    RANGE_GATE(beam) = find(abs(Alt(:,beam) - HEIGHT) == min(abs(Alt(:,beam) - HEIGHT)));
    vlos_all(:,beam,:)  = fits(4,3,:               ,beam,:); % Fit 4 = vlos, Species 3 = electrons
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

for tr = tstart:tend
%% Grab the current ISR parameters for time tr
    vlos = vlos_all(:,:,tr); % Just the current time
    vlos = vlos(range_idx);  % Just the range gates represented by xgmag, ygmag, zgmag
    
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
                Pv = diag([3000, 3000, 15].^2);
               
                vest(y,x,:) = Pv*A2'*( (A2*Pv*A2' + Pe) \ vlos_current );
                
                Pvest = Pv - Pv*A2'*((A2*Pv*A2' + Pe)\(A2*Pv));
                
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
            
            % DISPLAY ERROR ELLIPSE    
            subplot(Ny,Nx,(y-1)*Ny+x),
            error_ellipse(Pvest(1:2,1:2),'conf',0.394);
            xlim([-200,200]);
            ylim([-150,150]);
            set(gca,'XTick',[-200,-100,0,100,200]);
            set(gca,'XTickLabel',{});
            set(gca,'YTick',[-150,-75,0,75,150]);
            set(gca,'YTickLabel',{});
            axis equal;
%             set(gca,'XTick',[],'YTick',[]);
            if y == Ny,
                xlabel('v_{pe} (m/s)');
                set(gca,'XTickLabel',[-200,-100,0,100,200]);
            end;
            if x == 1,
                ylabel('v_{pn} (m/s)');
                set(gca,'YTickLabel',[-150,-75,0,75,150]);
            end;

        end; % for y
    end; % for x
    
    outfilename = fullfile(outputdir,sprintf('err_%03d',tr));
    print('-f1','-dpng','-r80',outfilename);
    print('-f1','-depsc2','-r80',outfilename);

end; % tr (radar time indices loop)