% Test 06: Plot line-of-sight velocities versus time and beam position.
%
% Thomas Butler
% 12 May 2010

clear all; clc;

StartTime = datenum(2009,03,24,07,49,57);
EndTime   = datenum(2009,03,24,10,00,03);

DateStr   = datestr(StartTime,'yyyymmdd');

Tint     = '30sec'; % Integration time
h5path   = '/shared/classes/projects/AMISR-experiments/data/';
h5subdir = [DateStr,'.001/'];
h5path   = fullfile(h5path,h5subdir);
h5file   = [DateStr '.001_lp_' Tint '.h5'];
h5file   = fullfile(h5path,h5file);

%%
% PlotStartTime = datenum(2009,03,24,08,03,00);
% PlotEndTime = datenum(2009,03,24,08,06,00);
% PlotName = 'Fig11';

PlotStartTime = datenum(2009,03,24,08,15,30);
PlotEndTime = datenum(2009,03,24,08,19,30);
PlotName = 'Fig12';

% PlotStartTime = datenum(2009,03,24,08,29,00);
% PlotEndTime = datenum(2009,03,24,08,45,30);
% PlotName = 'Fig13';

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
%     RANGE_GATE(beam) = find(abs(Alt(:,beam) - HEIGHT) == min(abs(Alt(:,beam) - HEIGHT)));
    
    
    % Actual range (km) for each radar beam corresponding to approximately HEIGHT altitude
%     rr(beam) = Range(RANGE_GATE(beam),beam); 
%     % Actual range (km) for each radar beam in the fifth range gate
%     rr(beam) = Range(5,beam);
    
    
    % Gather vlos & Ti (and errors) as a function of range, beam #, and time.
    vlos_all(:,beam,:)  = fits(4,3,:               ,beam,:); % Fit 4 = vlos, Species 3 = electrons
%     Ti_all(beam,:)      = fits(2,1,RANGE_GATE(beam),beam,:); % Fit 2 = Temp, Species 1 = hydrogen (same temp for all ions)
    dvlos_all(:,beam,:) = errs(4,1,:               ,beam,:); % Species 1, fitter sets all vlos equal, so only one error
%     dTi_all(beam,:)     = errs(2,1,RANGE_GATE(beam),beam,:);

end;

% Convert Unix time to MATLAB time
mtime = unixtime2matlab(utime,0);

%%
beam = [4 14 12 15 7 16 17 11 18 19 3 20 10 21 6 2 22 9 23 24 1 25 8 26 5 13];

figure(1);
clf;

facecolors = jet(5); % 5 == maximum length of the vector i for the given altitude range (see below)

for b = 1:26,
    [i,j] = find(Alt > 200e3 & Alt < 350e3);
    i = i(j==beam(b));
    j = j(j==beam(b));

    if b<26,
        subplot(6,5,b);
    else
        subplot(6,5,28);
    end;
%     plot(mtime(1,:),squeeze(vlos_all(i,beam(b),:)));
    for i_idx = length(i):-1:1,
        this_i = i(i_idx);
        h = area(mtime(1,:),[squeeze(vlos_all(this_i,beam(b),:)-dvlos_all(this_i,beam(b),:)),...
                             2 * squeeze(dvlos_all(this_i,beam(b),:))],...
                'BaseValue',-10000);
        set(h(1),'FaceColor','none','LineStyle','none');  % Make the base of the area object invisible
        set(h(2),'FaceColor',facecolors(i_idx,:));        % Area plot of the 1 sigma confidence interval
        hold on;
    end;
    hold off;
    xlim([PlotStartTime,PlotEndTime]);
    datetick('x','MM','keeplimits');
    
    hline(0,'k--');
    ylim([-1000,1500]);
end;

subplot(6,5,30);
set(gca,'NextPlot','replacechildren','ColorOrder',facecolors)
handle = plot(mtime(1,:),squeeze(vlos_all(i,beam(b),:)));
axis off;
legend(handle,{'215km','250km','287km','323km','359km'});

papersize=[8.5,11];
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',papersize);
set(gcf,'PaperPosition',[0,0,papersize]);

if ~exist('test06_results','dir'),
    mkdir('test06_results');
end;
print('-f1','-deps2c','-r150',fullfile('test06_results',['vlos_',PlotName]));