function makemovie_AMISR_vlos(filelab,itstart,itfin,opflag)

%This script makes plots of the usual plasma parameters in 3D slice
%configuration.  Also makes a movie, be sure to set the operating system flag...


load(['./datasets/data_mat/',filelab,'_rawdata.mat']); 
load(['./datasets/data_mat/',filelab,'_fieldgrid.mat']);


outdir=['./plot_imgfiles/',filelab,'/PFISR/'];
if ~exist(outdir)
    mkdir(outdir);
else
    delete([outdir,'/*.png'])
end

% MAKE THE PLOTS AND MOVIE
hat=figure(1);
set(gcf,'PaperPosition',[0 0 8 4.5])
if opflag==1
    vidObj = VideoWriter([outdir,'/PFISR.avi']);
    set(vidObj,'FrameRate',5);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

for i1=itstart:itfin; % runs for the time window specified
    
    data1=isne(:,:,i1);
    data2=isti(:,:,i1);
    data3=iste(:,:,i1);
    data4=isvi(:,:,i1);
    
    inds=find(data1<=0); % get rid of negative density, updated: 2017/03/31
    data1(inds)=NaN;
    data2(inds)=NaN;
    data3(inds)=NaN;
    data4(inds)=NaN;

    if (min(isnan(data1))==0)   % check that we have at least one good data point in the record otherwise pointless to plot
        title1=datestr(exp_date(i1,:));
        sliceplot_params_vlos(Rx,Ry,Rz,data1,data2,data3,data4,title1);
        
        % Save the plot
        filestr=datelab(exp_date(i1,4)+exp_date(i1,5)/60+exp_date(i1,6)/3600,[exp_date(i1,3),exp_date(i1,2),exp_date(i1,1)]);
        print([outdir,filestr,'.png'],'-dpng','-r300');

if opflag==1        
         % Write each frame to the file.
         currFrame = getframe(hat);
         writeVideo(vidObj,currFrame);
end        
    end
end

if opflag==1
    close(vidObj);
elseif opflag~=1
    system(['sh vid_encode.sh ',outdir,' ',outdir,'/PFISR.avi 2']);    %only works in linux or mac os with mplayer installed
end

end



