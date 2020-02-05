%-------------------------------------------------------------------------%
% Unwarp all allsky images for a given date10 Nov 2007.  Since unwarping is by far the
% most time-consuming part, I only want to do it once, and retrieve the
% results when plotting.
% 
% Thomas Butler
% 17 April, 2009
% 23 July, 2009
%  6 January, 2009 (changed x/y limits to correspond to higher altitude
%                   velocity field projection.
%
%-------------------------------------------------------------------------%

function unwarp_all(DateStr,Tint)

% Have to deal with inconsistent filename conventions. :/
if strcmp(DateStr(1:4),'2007'),
    is_2007 = 1;
else
    is_2007 = 0;
end;

% Paths & filenames
h5path = fullfile('/shared/classes/projects/AMISR-experiments/data/',[DateStr '.001']);
h5file = fullfile(h5path,[DateStr '.001_lp_' Tint '.h5']);

if is_2007,
    image_path = ['/shared/classes/projects/AMISR-experiments/data/' DateStr '_Allsky/FITS_' DateStr '/'];
    image_file_list = dir([image_path 'AS*.FITS']);
else
    image_path = ['/shared/classes/projects/AMISR-experiments/data/' DateStr '_Allsky/FITS_' DateStr '/'];
    image_file_list = dir([image_path 'PKR_DASC*.FITS']);
end;

% Are there more than one wavelenth represented by the image files?  If so,
% determine the wavelengths.
if is_2007,
    wavelength_list = 0;
else
    for t = 1:length(image_file_list),
        A(t) = sscanf(image_file_list(t).name,'PKR_DASC_%4d*');
    end;
    if all(A) == A(1),
        wavelength_list = A(1);
    else
        wavelength_list = unique(A);
    end;
end;

%% Gather data from HDF5 file
bco   = hdf5read(h5file,'BeamCodes');
utime = hdf5read(h5file,'Time/UnixTime');

%% Format the data as we need it...

az = bco(2,:) * pi/180;
el = bco(3,:) * pi/180;

mtime = unixtime2matlab(utime,0);
clear bco fits utime;

% Determine the main loop limits (indexes into mtime)
tstart = 1; %find(mtime(1,:) >= StartTime,1,'first');
tend   = size(mtime,2); %find(mtime(2,:) <= EndTime,  1,'last');

% Some other parameters
z0 = 120; % kilometers

if is_2007,
    lims = [-100,100,-50,125]; % Good for 11 November 2007
else
    lims = [-150,200,-50,300]; % Good for 2008 - 2009
end;


for wavelength = wavelength_list,

    wavelength_str = sprintf('%04d',wavelength);
    if numel(wavelength_list) > 1,
        image_file_list = dir([image_path 'PKR_DASC_' wavelength_str '_*.FITS']);
    end;
    
    %-------------------------------%
    % TIME VECTOR FOR ALLSKY IMAGES %
    %-------------------------------%
    itime = zeros(1,length(image_file_list));
    if is_2007,
        for t = 1:length(image_file_list),
            A = sscanf(image_file_list(t).name,'AS%2d%2d%2d.FITS');
            [yyyy,mm,dd] = datevec(DateStr,'yyyymmdd');
            A = [yyyy; mm; dd; A];
            itime(t) = datenum(A');
        end;
    else
        for t = 1:length(image_file_list),
            A = sscanf(image_file_list(t).name,['PKR_DASC_' wavelength_str '_%2d%2d%2d_%2d%2d%2d.*']);
            A(1) = A(1) + 2000;
            itime(t) = datenum(A');
        end;
    end;

    duw_all = zeros(lims(4)-lims(3)+1,...
                    lims(2)-lims(1)+1,...
                    length(image_file_list));
    itime_all = zeros(size(itime));
    index = 0;

    %% Loop over the radar (H5) time indices
    for tr = tstart:tend
        %% Load and unwarp allsky image
        this_itime = find( (itime >= mtime(1,tr)) & (itime <= mtime(2,tr)) );
        for ti = this_itime,
            image_file = image_file_list(ti).name;

            tic
            duw = Warp_Allsky([image_path image_file],z0,lims);
            toc
    %         [duw,lims] = Warp_Allsky([image_path image_file],z0);

            index = index+1;

            duw_all(:,:,index) = duw;
            itime_all(index) = itime(ti);

            disp(datestr(itime(ti)));

        end; %ti (image time indices loop)

    end; % tr (radar time indices loop)

    output_file = ['/shared/classes/projects/AMISR-experiments/data/' DateStr '_Allsky/UnwarpedImages.' wavelength_str '.mat'];

    itime = itime_all;
    
    eval(['save ' output_file ' duw_all itime lims']);
    
end; % for wavelength
