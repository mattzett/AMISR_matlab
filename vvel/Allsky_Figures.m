%Produce Allsky summary image figure for 2009 JGR paper:
%Read in set of allsky images, and el/az information for both images and PFISR.
%produce figures with PFISR beams inset, include time and fiducials.
%
%written by: Joshua Semeter
%last modified: 9 June 2009

path_r='/Users/jls/Documents/PROJECTS/Butler/V_inversion/PFISR_20080326/';
path_i='/Users/jls/Documents/PROJECTS/Butler/V_inversion/Allsky_20080326/';
fn_i=dir([path_i '*000.FITS']);     %make list of all images in directory

% Read in geometry info:  beam positions and el/az for images
fn_r='20080326.001_lp_2min.h5';
BeamCodes=cast(hdf5read([path_r fn_r],'BeamCodes'),'double');
az_r=BeamCodes(2,:);
el_r=BeamCodes(3,:);
%adjust AZ so it goes 0 to 360.
az_r360=az_r;
neg=find(az_r<0);
az_r360(neg)=360+az_r(neg);

% read el/az map for allsky image
fn_az='PKR_DASC_AZ_Nov2007.FITS';
az=fitsread([path_i fn_az]);
azf=rot90(rot90(az));
fn_el='PKR_DASC_EL_Nov2007.FITS';
el=fitsread([path_i fn_el]);
elf=rot90(rot90(el));
%Make North exactly in vertical direction by rotating images
r=-azf(100,1024/2);  %angle needed to accomplish this.
r=-27.;
azf=imrotate(azf,r);
elf=imrotate(elf,r);

%compute pixel positions for pfisr beams
x=[]; y=[];
for i=1:length(az_r)
    cost=abs(azf-az_r360(i))+abs(elf-el_r(i));
    [xtemp,ytemp]=find(cost==min(min(cost)));
    x=[x xtemp(1)]; y=[y ytemp(1)];
end
 
%COLOR MAP EXPERIMENTS
%funky gray/color colormap
cg=colormap(gray);
cgn=imresize(cg,[32,3]);
for i=1:3 cgn(:,i)=decimate(cg(:,i),2); end
cj=colormap(jet);
cjn=imresize(cj,[32,3]);
for i=1:3 cjn(:,i)=decimate(cj(:,i),2); end
cn=[cgn;cjn];
cn(find(cn<=0))=0;
%jet with grayscale lower end.
cj=colormap(jet);
N=6;
for i=1:3 cj(1:N,i)=0:1/(N-1):1; end
%alas, brightened grayscale map seems best
colormap(gray); brighten(0.4)

%Make panels for figure
hf=figure(1); 
set(gcf,'InvertHardcopy','off','color','white')
hold off; 
j=0;
label=['a','b','c','d','e','f','g','h'];
%for i=[6,7,8,9,10,11,12,13]   %growth phase period
for i=[48,49,50,51,52,53,54,55]  %expansion phase period
    clf;
    fn=fn_i(i,1).name;
    d=fitsread([path_i fn]);
    df=rot90(rot90(d));  %put in same orientation as Vi's
    df=imrotate(df,r);  %rotate images to put north towards top of page
    imagesc(df);
    j=j+1;
    axis off;
    caxis([310,6000]);
    axis equal;
    axis off;
    hold on;
    plot(y,x,'X','Color','w'); %place a symbol at each pfisr beam
    time=strcat(fn(22:25),':',fn(26:27));
    text(230,200,strcat(label(j),')  ',time,' UT'),'color','white','fontsize',14)
    if j==1 
        annotation(hf,'arrow',[0.6197 0.621],[0.7141 0.768],'HeadWidth',6,...
        'Color',[1 1 1]);
        annotation(hf,'arrow',[0.6197 0.6734],[0.7141 0.7141],'HeadWidth',6,...
        'Color',[1 1 1]);    
        text(1020,210,'N','color','w','fontsize',14);
        text(1070,340,'E','color','w','fontsize',14);
    end
    colorbar;
    outfile=strcat(fn(22:27),'.pdf')
    saveas(hf,outfile,'pdf');
end
pause;
