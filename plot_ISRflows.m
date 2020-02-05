function plot_ISRflows(filelab,itstart,itfin,opflag)

%This script plots the ISR data at a set altitude along with the flow fields


load(['./datasets/data_mat/',filelab,'_rawdata.mat']);
load(['./datasets/data_mat/',filelab,'_fieldgrid.mat']);
load(['./datasets/data_mat/',filelab,'_vvelscans.mat']);


%SIZE OF DATA SET
[lz,lb,lt]=size(isne);


%BEAMS WITH FAVORABLE FIELD ALIGNMENT;
ibg=1:size(Rx,2);    %for now assume that we are using all of the beams, might not be used


%CLEAN UP OUTPUT DIRECTORY
outdir=['./plot_imgfiles/',filelab,'/ISRflows/'];
if ~exist(outdir)
    mkdir(outdir);
else
    delete([outdir,'/*.png'])
end


% SETUP THE PLOTS AND MOVIE
hat=figure;
set(gcf,'PaperPosition',[0 0 11 4.5])
if opflag==1
    vidObj = VideoWriter([outdir,'/ISRflows.avi']);
    set(vidObj,'FrameRate',5);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end


%DEFINE ALTITUDES OF INTEREST AND CORRESPONDING GRID
ngrid=100;
ix2=zeros(1,lb);
for ib=1:lb
    [~,iz2now]=min(abs(Rz(:,ib)-300));   %range index where the altitude hits 300 km for this beam
    iz2(ib)=iz2now;
    Rxset(ib)=Rx(iz2(ib),ib);            %corresponding value for x at this location
    Ryset(ib)=Ry(iz2(ib),ib);            %ditto for y
end
Rxset=Rxset(:); Ryset=Ryset(:);
x=linspace(min(Rxset),max(Rxset),ngrid);
y=linspace(min(Ryset),max(Ryset),ngrid);
[X,Y]=meshgrid(x,y);


%REMOVE NEGATIVE DENSITIES AND RELATED DATA
Neg = find(isne<1E-100);
isne(Neg)=NaN;
isti(Neg)=NaN;
isvi(Neg)=NaN;


%LOOP OVER SCANS AND MAKE PLOTS
ndigits=ceil(log10(lt));
zcisr=zeros(lb,lt);
Tiavg=zeros(1,lt);
vmagavg=zeros(1,lt);

for it=1:itfin-itstart  % runs for the time window specified
    %ELECTRON DENSITY
    for ib=1:lb
        nenow(ib)=isne(iz2(ib),ib,it);
    end
    nenow=nenow(:);
    ne2FN=TriScatteredInterp(Rxset,Ryset,nenow);
    ne2i=ne2FN(X(:),Y(:));
    ne2i=reshape(ne2i,size(X));
    
    
    %ION TEMPERATURE
    for ib=1:lb
        Tinow(ib)=isti(iz2(ib),ib,it);
    end
    Tinow=Tinow(:);
    Ti2FN=TriScatteredInterp(Rxset,Ryset,Tinow);
    Ti2i=Ti2FN(X(:),Y(:));
    Ti2i=reshape(Ti2i,size(X));
    
    
    %PLOTTING
    clf;

    subplot(131)
    h=imagesc(x,y,Ti2i);
    set(h,'alphadata',~isnan(Ti2i));
    axis xy;
    axis square;
    xlabel('dist. E of PFISR (km)');
    ylabel('dist. N of PFISR (km)');
    %     caxis([0, 3000])
    c=colorbar;
    ylabel(c,'T_i (K) at 300 km')
    
    ctriple=[0 0 0];
    
    t=exp_date(itstart+it-1,4)*3600+exp_date(itstart+it-1,5)*60+exp_date(itstart+it-1,6);
    title1=datestr(exp_date(itstart+it-1,:));
    
    %ADD VVELS TO TEMPERATURE PLOT
    hold on;
    %    qfac = 0.03;
    qfac=0.01;
    veast=squeeze(vest_geog(:,:,1,it));
    vnorth=squeeze(vest_geog(:,:,2,it));
    h=quiver(Xvg(:),Yvg(:),veast(:)*qfac,vnorth(:)*qfac,'Color','w','LineWidth',1);
    %h.AutoScale='off';
    h2=quiver(40, -15, 1e3*qfac, 0, 0,'Color',ctriple,'LineWidth',1);
    %h2.AutoScale='off';
    text(40, -10, '1 km/s','color',ctriple,'FontWeight','bold');
    %text(x(floor(ngrid/10)),y(floor(ngrid-ngrid/10)),[num2str(t/3600),' UT'], ...
    %    'FontSize',12,'Color',ctriple,'FontWeight','bold');
    title(title1);    
    hold off;
    
    
    subplot(132)
    h=imagesc(x,y,log10(ne2i));
    set(h,'alphadata',~isnan(ne2i));
    axis xy;
    axis square;
    %    xlabel('dist. E of PFISR (km)');
    %    ylabel('dist. N of PFISR (km)');
    %     caxis([10.9 11.7]);
    c=colorbar;
    ylabel(c,'log_{10} n_e at 300 km')
    
    
    subplot(133)
    vmag=sqrt(veast.^2+vnorth.^2);
    imagesc(vmag);
    %     caxis([0 3000])
    colorbar
    title('|v| debug')
    axis xy;
    axis square;
    
    
    %PRINT THE OUTPUT TO GRAPHICS FILE
    udigits=ceil(log10(it+1));
    lzeros=[];
    for id=1:ndigits-udigits
        lzeros=[lzeros,'0'];
    end
    print([outdir,lzeros,num2str(itstart+it-1),'.png'],'-dpng','-r300');
    
if opflag==1 
    % Write each frame to the file.
    currFrame = getframe(hat);
    writeVideo(vidObj,currFrame);
end    
    
    
    %STORE AVERAGE VALUES FOR EACH FRAME
    Titmp=Ti2i(~isnan(Ti2i));
    Tiavg(it)=mean(Titmp(:));
    Titmp=vmag(~isnan(vmag));
    vmagavg(it)=mean(vmag(:));

end

if opflag==1
    close(vidObj);
elseif opflag~=1
    system(['sh vid_encode.sh ',outdir,' ',outdir,'/ISRflows.avi 2']);    %only works in linux or mac os with mplayer installed
end


%SOME SUMMARY STATISTICS
figure;
plot(vmagavg,Tiavg,'o');
xlabel('|v| @ 300 km');
ylabel('Ti @ 300 km');

print([outdir,'circles.png'],'-dpng','-r300');


