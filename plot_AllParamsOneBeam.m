clear%, close all

%Note: 15 beams, 24 altitudes.

%Establish variables
BNum = 1; % BEAM NUMBER OF INTEREST

%%%%%%%%%% Plot Altitude Range %%%%%%%%%%%%%
altrange  = [100 500];
cne       = [10.0 11.5];
cTi       = [100 2000];
%cTsperp  = [800 3300];
cv        = [-700 700];
cTe       = [100 3500];
sized     = 18;

%%%%%%% Plot Altitude Error Range %%%%%%%%%%
cdne       = [0 5.0];
cdTi       = [0 750];
cdv        = [0 800];
cdTe       = [0 3500];

%%%%%%%Introduce filelab to this code%%%%%%%
filelab='22Feb2017_2.1min';
% filelab='02Mar2017_2.1min';
fprintf(strcat('AMISR --> File selected: ',filelab,'\n'))

outdir=['./plot_imgfiles/',filelab];
if ~exist(outdir)
    mkdir(outdir);
end

%%%%%%%%%%% Chosen Data %%%%%%%%%%%%%%%%%%%
load(['./datasets/data_mat/',filelab,'_rawdata.mat']);

%REMOVE NEGATIVE DENSITY VALUES AND RELATED DATA
Neg = find(isne<1E-100);
isne(Neg)=NaN;
isti(Neg)=NaN;
isvi(Neg)=NaN;

%Convert experiment date to time in seconds
t=exp_date(:,4)*3600+exp_date(:,5)*60+exp_date(:,6);
time=t/3600/24;

figure(5)
set(figure(5),'PaperUnits','inches','PaperPosition',[0,0,12,20])

subplot(4,1,1)
% h=pcolor(time,isz(:,BNum),(log10(squeeze(isne(:,BNum,:)))));
% shading flat
h=imagesc(time,isz(:,BNum),log10(squeeze(isne(:,BNum,:))));
axis xy
set(h,'AlphaData',log10(squeeze(isdne(:,BNum,:)))<3);
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% day=datenum(dateyymmdd,'yymmdd');
% title(sprintf( datestr(day, 'mm/dd/yyyy')))
%ylim(altrange);
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
clb=colorbar('location','EastOutside');
caxis(cne); 
ylabel(clb,'n_e (m^{-3})');
title(['az: ',num2str(az(BNum)),', el: ',num2str(el(BNum)),', BeamCode: ',num2str(beamcodes(BNum))]);

subplot(4,1,2)
% h=pcolor(time,isz(:,BNum),squeeze(isti(:,BNum,:)));
% shading flat
h=imagesc(time,isz(:,BNum),squeeze(isti(:,BNum,:)));
axis xy
set(h,'AlphaData',squeeze(isdti(:,BNum,:))<1000);
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% ylim(altrange); 
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
caxis(cTi); 
clb=colorbar('location','EastOutside');
ylabel(clb,'Ion Temp (K)');

subplot(4,1,3)
% h=pcolor(time,isz(:,BNum),squeeze(iste(:,BNum,:)));
% shading flat
h=imagesc(time,isz(:,BNum),squeeze(iste(:,BNum,:)));
axis xy
set(h,'AlphaData',squeeze(isdte(:,BNum,:))<500);
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% ylim(altrange); 
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
caxis(cTe); 
clb=colorbar('location','EastOutside');
ylabel(clb,'Electron Temp (K)');

subplot(4,1,4)
% h=pcolor(time,isz(:,BNum),squeeze(isvi(:,BNum,:)));
% shading flat
h=imagesc(time,isz(:,BNum),squeeze(isvi(:,BNum,:)));
axis xy
set(h,'AlphaData',squeeze(isdvi(:,BNum,:))<200);
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% ylim(altrange); 
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
caxis(cv); 
clb=colorbar('location','EastOutside');
ylabel(clb,'v_{los} (m/s)');
xlabel('Time (UT)');

% print -dpng AllParamsOneBeam.png;

%%%%%%%%%%%%%%Plot error values%%%%%%%%%%%%%%%%%

figure(6)
set(figure(6),'PaperUnits','inches','PaperPosition',[0,0,12,20])
subplot(4,1,1)
pcolor(time,isz(:,BNum),real(log10(squeeze(isdne(:,BNum,:)))));
shading flat
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% day=datenum(dateyymmdd,'yymmdd');
% title(sprintf( datestr(day, 'mm/dd/yyyy')))
%ylim(altrange);
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
clb=colorbar('location','EastOutside');
caxis(cdne); 
ylabel(clb,'n_e (m^{-3})');
title(['Error Bars, ','az: ',num2str(az(BNum)),', el: ',num2str(el(BNum)),', BeamCode: ',num2str(beamcodes(BNum))]);

subplot(4,1,2)
pcolor(time,isz(:,BNum),squeeze(isdti(:,BNum,:)));
shading flat
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% ylim(altrange); 
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
caxis(cdTi); 
clb=colorbar('location','EastOutside');
ylabel(clb,'Ion Temp (K)');

subplot(4,1,3)
pcolor(time,isz(:,BNum),squeeze(isdte(:,BNum,:)));
shading flat
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% ylim(altrange); 
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
caxis(cdTe); 
clb=colorbar('location','EastOutside');
ylabel(clb,'Electron Temp (K)');

subplot(4,1,4)
pcolor(time,isz(:,BNum),squeeze(isdvi(:,BNum,:)));
shading flat
set(gca,'FontSize',sized,'XMinorTick','on','YMinorTick','on');
% ylim(altrange); 
xlim([6.5/24 9.5/24])
ylim([min(isz(:,BNum)) 400]);
datetick('x','keeplimits');
ylabel('Altitude (km)');
caxis(cdv); 
clb=colorbar('location','EastOutside');
ylabel(clb,'v_{los} (m/s)');
xlabel('Time (UT)');

% print -dpng AllParamsOneBeam_errors.png;
