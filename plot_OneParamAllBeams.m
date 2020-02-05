%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Reads PFISR input data and plots each paramter for each beam%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

addpath ./geom ./vvel;

%Pick the data to use
filelab='22Feb2017_2.1min';
 
fprintf(strcat('AMISR --> File selected: ',filelab,'\n'))


%Load input PFISR data
load(['./datasets/data_mat/',filelab,'_rawdata.mat']);
load(['./datasets/data_mat/',filelab,'_fieldgrid.mat']);
load(['./datasets/data_mat/',filelab,'_vvelscans.mat']);


%Create variable for number of beams
NumberofBeams=size(isne,2);


%Convert time
t=(exp_date(:,4)*3600+exp_date(:,5)*60+exp_date(:,6))/3600;


%REMOVE NEGATIVE VALUES FROM ISNE, ISTI, ISVI TO REDUCE POTENTIAL ERROR
Neg = find(isne<1E-100);
isne(Neg)=NaN;
isti(Neg)=NaN;
isvi(Neg)=NaN;


%Create figure window
figure(1);
set(gcf,'PaperPosition',[0 0 11 11]);


%Create variable for number of subplots, depends on number of beams
NRowsAndCols=(ceil(sqrt(NumberofBeams+1)));


%Create loop of subplots
for in=1:NumberofBeams
    subplot(NRowsAndCols,NRowsAndCols,in);
    imagesc(t,isz(:,in),(squeeze(isti(:,in,:)))); %isne, isp, iste, isti, isvi Different parameter choices
    set(h,'AlphaData',squeeze(isdti(:,in,:))<500);
%     caxis([1E10 4E11]) %for isne
%     caxis ([0 1]) %for isp
%     caxis ([0 4000]) %for iste
    caxis([0 4000]) %for isti
%     caxis ([-200 200]) %for isvi
    clb = colorbar;
    xlim(timeframe);
    ylim([min(isz(:,in)) 400]);
    title(['az: ',num2str(az(in)),', el: ',num2str(el(in)),', BeamCode: ',num2str(beamcodes(in))]);
    xlabel('Time (UT)');
    ylabel('Altitude (km)');
    ylabel(clb,'isti'); %isne, isp, iste, isti, isvi Different parameter choices
    axis xy;
    %axis square;
end

%%Add plot beams and mag field
%plot_grid(filelab)


