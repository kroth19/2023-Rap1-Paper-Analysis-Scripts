%% read in images and annotations for wound healing experiments and calculate properties of wound over time

clear all;


%% specify image sets to analyze
con = {
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo7/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201002/Embryo1/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201002/Embryo3/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo1/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo2/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo3/annotations_woundedge.mat',...
    };

Rap1DN = {
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo2/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo3/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo4/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo1/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo2/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo3/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo6/annotations_woundedge.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo5/annotations_woundedge.mat',...
    };

dirs1_orig = con;
dirs2_orig = Rap1DN;
legend_str = {'control', 'Rap1DN'};
label488 = 'sqh';
label561 = 'Ecad';
channel488 = 'kr-C488-D';
channel561 = 'kr-C561-D';


%% Parameters for analysis

% Define directories to analyze/plot
dirs1 = dirs1_orig;
dirs2 = dirs2_orig;

% Flags for analysis & plotting
analyze_flag = 0; % Run analysis + save
intensity_flag = 1; % Run intensity section of the analysis (unless analyze_flag == 0)
plot_flag = 1; % Show plots
save_flag = 0; % Save plots
black_bg = 0; % If the pixel value of background pixels is strictly zero, set this to one. This affects the normalization of fluorescence over time in images that have been registered and where black pixels have been introduced.
normalize_intensity_flag = 2; % This is normally 0. Use 1 in case of large bleaching (e.g. dirs_sqhmcherry_moegfp_wt) to normalize each time point to its mean fluorescence. Use 2 to substract the background (image mode) and normalize each time point to the mean fluorescence of the initial time points (1:cut_t).

save_area = 1; % Save areas
extension = 'mat'; % Extension format to save the areas - options are mat or csv

% Define key variables
t_res = 30; % Time resolution in seconds
xyres = 16/(60*1.5); % Microns per pixel: 0.33 for 40X, 0.2 for 63X and 0.13 for 100X (assuming binning of 2 in all cases).
ntp = 91; % Number of time points - should exceed max # points in all data sets
% ntp = 61;
the_cut_t = 2; % Number of time points that occur prior to wounding
timeshifts = -16:16; % Time shifts (in number of time points) used to calculate the best correlation between area and intensity.
t_range = [-5 50]; % Time range for plotting (in minutes)

t = ((1:ntp)-the_cut_t-1).*t_res/60; % Time in minutes.

brush_sz = 3;

% Define plotting features
color_str = 'brk';
marker_str = '^s';
line_str = {'-'; '-'; '-'};
title_label = '';

%% Run annotation analysis in batch mode.
if analyze_flag
    
    hwb = waitbar(0, cat(2, 'Analyzing movie 1/', cat(2, num2str(numel(dirs1)), ' ...')));
    
    for ii = 1:numel(dirs1)
        if exist(dirs1{ii}, 'file') % check whether a file with this name exists
            load(dirs1{ii});
            indsep = strfind(dirs1{ii}, filesep);
            analysis_file = strrep(dirs1{ii}(indsep(end)+1:end), 'annotations', 'analysis');
            if strcmp(analysis_file, dirs1{ii}(indsep(end)+1:end)) % if annotations filename doesn't include "annotations"
                analysis_file = 'analysis_healing.mat';
            end
            dirs1{ii} = dirs1{ii}(1:indsep(end));
            
        elseif exist(cat(2, dirs1{ii}, 'annotations_healing.mat'), 'file') % generic annotations filename
            load(cat(2, dirs1{ii}, 'annotations_healing'));
            analysis_file = 'analysis_healing.mat';
        else % no detectable annotations file
            disp(cat(2, 'Error when loading ', dirs1{ii}));
        end
        
        cd(dirs1{ii}); % make folder containing annotations file the current directory
        
        first = 1; % pre-define first time point
        while isempty(ud.rpolygons{1, 1, first}) % check for first time point with annotations
            first = first + 1;
        end
        
        last = size(ud.rpolygons, 3); % pre-define last time point
        while isempty(ud.rpolygons{1,1,last})
            last = last - 1;
        end
        
        p1 = [];
        p2 = [];
        
        % Find the time point immediately before wounding by the presence of a "unique" marker.
        if isfield(ud,'runiquecoords')
            cut_t = find(ud.runiquecoords(:, 1) ~= -1);
        else
            cut_t = [];
        end
        
        if ~isempty(cut_t)
            cut_t = cut_t(1);
        else
            % This assumes that the same number of time points is always annotated before wounding.
            cut_t = the_cut_t - first + 1;
        end
        
        if cut_t < 1
            cut_t = 1; % WHAT IS CUT_T USED FOR? HOW IS IT DIFFERENT FROM THE_CUT_T?
        end
        
        % An empty image is passed so that, if any of the polygons self-intersects, the area, perimeter and so on can be measured on a mask on the dip_image.
        % THIS RETURNS MEASUREMENTS ONLY FOR THE TIME POINTS ANNOTATED.
        
        [l,a,p,~,lh,lv,lellipseh,lellipsev] = cut_analysis(p1, p2, ud.rpolygons(:, :, first:last), [], cut_t, t_res, 0, 0, newim(ud.imsize(1:2)));
        
        % Here we quantify image intensities.
        % The filenames will be of this format: folder_label_channel.tif
        indsep = strfind(dirs1{ii}, filesep);
        filebase = dirs1{ii}(indsep(end-1)+1:indsep(end)-1);
        
        % Allocate memory to store the measurements.
        i488_peri = nan .* zeros(1, last-first+1);
        i488_center = nan .* zeros(1, last-first+1);
        i488_immean = nan .* zeros(1, last-first+1);
        i488_immode = nan .* zeros(1, last-first+1);
        i568_peri = nan .* zeros(1, last-first+1);
        i568_center = nan .* zeros(1, last-first+1);
        i568_immean = nan .* zeros(1, last-first+1);
        i568_immode = nan .* zeros(1, last-first+1);
        m_corr_i488_peri = nan .* zeros(1, last-first+1);
        s_corr_i488_peri = nan .* zeros(1, last-first+1);
        m_corr_i568_peri = nan .* zeros(1, last-first+1);
        s_corr_i568_peri = nan .* zeros(1, last-first+1);
        
        if intensity_flag
            % Load 488 image.
            if exist(cat(2, dirs1{ii}, filebase,'_', label488, '_', channel488, '.tif'), 'file')
                im = tiffread(cat(2, dirs1{ii}, filebase,'_', label488, '_', channel488, '.tif'));
            elseif exist(cat(2,dirs1{ii},'_', channel488, '.tif'),'file')
                im = tiffread(cat(2,dirs1{ii},'_',channel488, '.tif'));
            else
                fprintf('File not found: %s',cat(2, dirs1{ii}, filebase,'_', label488, '_', channel488, '.tif'))
            end
            
            % Go through each slice.
            for jj = first:last
                % Create a mask with the wound outline.
                msk_peri = trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 0);
                mask_file = cat(2,dirs1{ii},'mask_', channel488, '.tif');
                msk_total = fillholes(msk_peri);
                cast_fn = @uint16;
                imwrite(cast_fn(msk_total),mask_file, 'tif', 'Compression', 'none', 'WriteMode', 'append');
                msk_center = fillholes(msk_peri) - msk_peri; %trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 1);
                
                % If the mask is not black.
                if nnz(double(msk_peri))
                    slice = squeeze(im(:, :, jj-1));
                    % Measure the mask mean intensity and the mean and mode of the slice.
                    i488_peri(jj-first+1) = mean(slice(msk_peri));
                    try
                        i488_center(jj-first+1) = mean(slice(msk_center));
                    catch
                        disp(jj);
                        dipshow(msk_peri);
                    end
                    % This deals with situations in which the images have been registered and black pixels introduced. But this will mess up images in which the background contains black pixels due to high signal to noise ratios. Thus, we introduce a new flag.
                    if ~black_bg
                        if numel(find(slice == 0)) > min(ud.imsize(1:2))
                            i488_immean(jj-first+1) = mean(slice(slice~=0));
                            disp('mean modified to del with zeros')
                        else
                            i488_immean(jj-first+1) = mean(slice);
                        end
                        
                    else
                        i488_immean(jj-first+1) = mean(slice);
                    end
                    themode = mode(reshape(double(slice), [1 prod(size(slice))]));
                    if themode == max(max(double(slice)))
                        h = diphist(slice, [0 65535], 65536);
                        [tmp,themode(1)] = max(h(1:(end-1)));
                        themode = themode(1) - 1;
                    elseif themode(1) == 0
                        h = diphist(slice, [0 65535], 65536);
                        [tmp,themode(1)] = max(h(2:(end)));
                        themode = themode(1);
                    end
                    i488_immode(jj-first+1) = themode;
                    % Use these 3 lines if you want the mode to be calculated by fitting a Gaussian to the histogram and taking the mean of the Gaussian.
                    %h = diphist(slice, [0 65535], 65536);
                    %[mu, sigma, scale, fwhm] = fitgaussian(h,[themode 100 max(h)], 0);
                    %i488_immode(jj-first+1) = round(mu);
                    m_corr_i488_peri(jj-first+1) = mean((slice(msk_peri)-i488_immode(jj-first+1))/i488_immean(jj-first+1));
                    s_corr_i488_peri(jj-first+1) = std((slice(msk_peri)-i488_immode(jj-first+1))/i488_immean(jj-first+1));
                end
            end
            
            % Load 568 image.
            if exist(cat(2, dirs1{ii}, filebase,'_',label561, '_', channel561, '.tif'), 'file')
                im = tiffread(cat(2, dirs1{ii}, filebase,'_',label561,'_',channel561,'.tif'));
            elseif exist(cat(2,dirs1{ii},'_', channel561, '.tif'),'file')
                im = tiffread(cat(2,dirs1{ii},'_',channel561, '.tif'));
                % Go through each slice.
                for jj = first:last
                    % Create a mask with the wound outline.
                    msk_peri = trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 0);
                    msk_center = fillholes(msk_peri) - msk_peri; %trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 1);
                    
                    % If the mask is not black.
                    if nnz(double(msk_peri))
                        slice = squeeze(im(:, :, jj-1));
                        % Measure the mask mean intensity and the mean and mode of the slice.
                        i568_peri(jj-first+1) = mean(slice(msk_peri));
                        i568_center(jj-first+1) = mean(slice(msk_center));
                        % This deals with situations in which the images have been registered and black pixels introduced. But this will mess up images in which the background contains black pixels due to high signal to noise ratios. Thus, we introduce a new flag.
                        if ~black_bg
                            if numel(find(slice == 0)) > min(ud.imsize(1:2))
                                i568_immean(jj-first+1) = mean(slice(slice~=0));
                                disp('mean modified to del with zeros')
                            else
                                i568_immean(jj-first+1) = mean(slice);
                            end
                            
                        else
                            i568_immean(jj-first+1) = mean(slice);
                        end
                        themode = mode(reshape(double(slice), [1 prod(size(slice))]));
                        if themode == max(max(double(slice)))
                            h = diphist(slice, [0 65535], 65536);
                            [tmp,themode(1)] = max(h(1:(end-1)));
                            themode = themode(1) - 1;
                        elseif themode(1) == 0
                            h = diphist(slice, [0 65535], 65536);
                            [tmp,themode(1)] = max(h(2:(end)));
                            themode = themode(1);
                        end
                        i568_immode(jj-first+1) = themode;
                        % Use these 3 lines if you want the mode to be calculated by fitting a Gaussian to the histogram and taking the mean of the Gaussian.
                        %h = diphist(slice, [0 65535], 65536);
                        %[mu, sigma, scale, fwhm] = fitgaussian(h,[themode 100 max(h)], 0);
                        %i568_immode(jj-first+1) = round(mu);
                        m_corr_i568_peri(jj-first+1) = mean((slice(msk_peri)-i568_immode(jj-first+1))/i568_immean(jj-first+1));
                        s_corr_i568_peri(jj-first+1) = std((slice(msk_peri)-i568_immode(jj-first+1))/i568_immean(jj-first+1));
                    end
                end
            end
        end
        
        save(analysis_file, 'l', 'a', 'p', 'lh', 'lv', 'lellipseh', 'lellipsev', 'cut_t', 'first', 'last', 't_res', 'i488_peri', 'i488_center', 'i488_immean', 'i488_immode','m_corr_i488_peri','s_corr_i488_peri', 'i568_peri', 'i568_center', 'i568_immean', 'i568_immode','m_corr_i568_peri','s_corr_i568_peri');
        
        if ii < numel(dirs1)
            waitbar((ii)/numel(dirs1), hwb, cat(2, 'Analyzing movie ', cat(2, cat(2, num2str(ii+1), filesep), cat(2, num2str(numel(dirs1)), ' ...'))));
        else
            waitbar(1, hwb, 'Done!');
        end
    end
    
    if ~ isempty(dirs2)
        hwb = waitbar(0, cat(2, 'Analyzing movie 1/', cat(2, num2str(numel(dirs2)), ' ...')));
        
        for ii = 1:numel(dirs2)
            if exist(dirs2{ii}, 'file')
                load(dirs2{ii});
                indsep = strfind(dirs2{ii}, filesep);
                analysis_file = strrep(dirs2{ii}(indsep(end)+1:end), 'annotations', 'analysis');
                if strcmp(analysis_file, dirs2{ii}(indsep(end)+1:end))
                    analysis_file = 'analysis_polygons.mat';
                end
                dirs2{ii} = dirs2{ii}(1:indsep(end));
                
            elseif exist(cat(2, dirs2{ii}, 'annotations_healing.mat'), 'file')
                load(cat(2, dirs2{ii}, 'annotations_healing'));
                analysis_file = 'analysis_healing.mat';
            else
                disp(cat(2, 'Error when loading ', dirs2{ii}));
            end
            
            cd(dirs2{ii});
            
            first = 1; % pre-define first time point
            while isempty(ud.rpolygons{1, 1, first}) % check for first time point with annotations
                first = first + 1;
            end
            
            last = size(ud.rpolygons, 3); % pre-define last time point
            while isempty(ud.rpolygons{1,1,last})
                last = last - 1;
            end
            
            p1 = [];
            p2 = [];
            
            % Find the time point immediately before wounding by the presence of a "unique" marker.
            if isfield(ud,'runiquecoords')
                cut_t = find(ud.runiquecoords(:, 1) ~= -1);
            else
                cut_t = [];
            end
            
            if ~isempty(cut_t)
                cut_t = cut_t(1);
            else
                % This assumes that the same number of time points is always annotated before wounding.
                cut_t = the_cut_t - first + 1;
            end
            
            if cut_t < 1
                cut_t = 1;
            end
            
            % An empty image is passed so that, if any of the polygons self-intersects, the area, perimeter and so on can be measured on a mask on the dip_image.
            
            [l,a,p,~,lh,lv,lellipseh,lellipsev] = cut_analysis(p1, p2, ud.rpolygons(:, :, first:last), [], cut_t, t_res, 0, 0, newim(ud.imsize(1:2)));
            % Here we quantify image intensities.
            % The filenames will be of this format: folder_label_channel.tif
            indsep = strfind(dirs2{ii}, filesep);
            filebase = dirs2{ii}(indsep(end-1)+1:indsep(end)-1);
            
            % Allocate memory to store the measurements.
            i488_peri = nan .* zeros(1, last-first+1);
            i488_center = nan .* zeros(1, last-first+1);
            i488_immean = nan .* zeros(1, last-first+1);
            i488_immode = nan .* zeros(1, last-first+1);
            i568_peri = nan .* zeros(1, last-first+1);
            i568_center = nan .* zeros(1, last-first+1);
            i568_immean = nan .* zeros(1, last-first+1);
            i568_immode = nan .* zeros(1, last-first+1);
            m_corr_i488_peri = nan .* zeros(1, last-first+1);
            s_corr_i488_peri = nan .* zeros(1, last-first+1);
            m_corr_i568_peri = nan .* zeros(1, last-first+1);
            s_corr_i568_peri = nan .* zeros(1, last-first+1);
            
            if intensity_flag
                % Load 488 image.
                if exist(cat(2, dirs2{ii}, filebase,'_', label488, '_', channel488, '.tif'), 'file')
                    im = tiffread(cat(2, dirs2{ii}, filebase,'_', label488, '_', channel488, '.tif'));
                elseif exist(cat(2,dirs2{ii},'_', channel488, '.tif'),'file')
                    im = tiffread(cat(2,dirs2{ii},'_',channel488, '.tif'));
                else
                    fprintf('File not found: %s',cat(2, dirs2{ii}, filebase,'_', label488, '_', channel488, '.tif'))
                end
                
                % Go through each slice.
                for jj = first:last
                    % Create a mask with the wound outline.
                    msk_peri = trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 0);
                    mask_file = cat(2,dirs2{ii},'mask_', channel488, '.tif');
                    msk_total = fillholes(msk_peri);
                    cast_fn = @uint16;
                    imwrite(cast_fn(msk_total),mask_file, 'tif', 'Compression', 'none', 'WriteMode', 'append');
                    msk_center = fillholes(msk_peri) - msk_peri; %trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 1);
                    
                    % If the mask is not black.
                    if nnz(double(msk_peri))
                        slice = squeeze(im(:, :, jj-1));
                        
                        % Measure the mask mean intensity and the mean and mode of the slice.
                        i488_peri(jj-first+1) = mean(slice(msk_peri));
                        try
                            i488_center(jj-first+1) = mean(slice(msk_center));
                        catch
                            disp(jj);
                            dipshow(msk_peri);
                        end
                        % This deals with situations in which the images have been registered and black pixels introduced. But this will mess up images in which the background contains black pixels due to high signal to noise ratios. Thus, we introduce a new flag.
                        if ~black_bg
                            if numel(find(slice == 0)) > min(ud.imsize(1:2))
                                i488_immean(jj-first+1) = mean(slice(slice~=0));
                                disp('mean modified to deal with zeros')
                            else
                                i488_immean(jj-first+1) = mean(slice);
                            end
                            
                        else
                            i488_immean(jj-first+1) = mean(slice);
                        end
                        %i488_immode(jj) =  mode(reshape(double(slice), [1 prod(size(slice))]));
                        themode = mode(reshape(double(slice), [1 prod(size(slice))]));
                        if themode == max(max(double(slice)))
                            h = diphist(slice, [0 65535], 65536);
                            [tmp,themode(1)] = max(h(1:(end-1)));
                            themode = themode(1) - 1;
                        elseif themode(1) == 0
                            h = diphist(slice, [0 65535], 65536);
                            [tmp,themode(1)] = max(h(2:(end)));
                            themode = themode(1);
                        end
                        i488_immode(jj-first+1) = themode;
                        % Use these 3 lines if you want the mode to be calculated by fitting a Gaussian to the histogram and taking the mean of the Gaussian.
                        %h = diphist(slice, [0 65535], 65536);
                        %[mu, sigma, scale, fwhm] = fitgaussian(h,[themode 100 max(h)], 0);
                        %i488_immode(jj-first+1) = round(mu);
                        m_corr_i488_peri(jj-first+1) = mean((slice(msk_peri)-i488_immode(jj-first+1))/i488_immean(jj-first+1));
                        s_corr_i488_peri(jj-first+1) = std((slice(msk_peri)-i488_immode(jj-first+1))/i488_immean(jj-first+1));
                    end
                end
                
                % Load 568 image.
                if exist(cat(2, dirs2{ii}, filebase,'_',label561, '_', channel561, '.tif'), 'file')
                    im = tiffread(cat(2, dirs2{ii}, filebase,'_',label561,'_',channel561,'.tif'));
                elseif exist(cat(2,dirs2{ii},'_', channel561, '.tif'),'file')
                    im = tiffread(cat(2,dirs2{ii},'_',channel561, '.tif'));
                    % Go through each slice.
                    for jj = first:last
                        % Create a mask with the wound outline.
                        msk_peri = trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 0);
                        msk_center = fillholes(msk_peri) - msk_peri; %trajectories2mask(ud.rpolygons(1, :, jj), ud.imsize(1:2), brush_sz, 1);
                        
                        % If the mask is not black.
                        if nnz(double(msk_peri))
                            slice = squeeze(im(:, :, jj-1));
                            % Measure the mask mean intensity and the mean and mode of the slice.
                            i568_peri(jj-first+1) = mean(slice(msk_peri));
                            i568_center(jj-first+1) = mean(slice(msk_center));
                            
                            % This deals with situations in which the images have been registered and black pixels introduced. But this will mess up images in which the background contains black pixels due to high signal to noise ratios. Thus, we introduce a new flag.
                            if ~black_bg
                                if numel(find(slice == 0)) > min(ud.imsize(1:2))
                                    i568_immean(jj-first+1) = mean(slice(slice~=0));
                                    disp('mean modified to del with zeros')
                                else
                                    i568_immean(jj-first+1) = mean(slice);
                                end
                                
                            else
                                i568_immean(jj-first+1) = mean(slice);
                            end
                            %i568_immode(jj) =  mode(reshape(double(slice), [1 prod(size(slice))]));
                            themode = mode(reshape(double(slice), [1 prod(size(slice))])); %#ok<*PSIZE>
                            if themode == max(max(double(slice)))
                                h = diphist(slice, [0 65535], 65536);
                                [tmp,themode(1)] = max(h(1:(end-1)));
                                themode = themode(1) - 1;
                            elseif themode(1) == 0
                                h = diphist(slice, [0 65535], 65536);
                                [tmp,themode(1)] = max(h(2:(end)));
                                themode = themode(1);
                            end
                            i568_immode(jj-first+1) = themode;
                            % Use these 3 lines if you want the mode to be calculated by fitting a Gaussian to the histogram and taking the mean of the Gaussian.
                            %h = diphist(slice, [0 65535], 65536);
                            %[mu, sigma, scale, fwhm] = fitgaussian(h,[themode 100 max(h)], 0);
                            %i568_immode(jj-first+1) = round(mu);
                            m_corr_i568_peri(jj-first+1) = mean((slice(msk_peri)-i568_immode(jj-first+1))/i568_immean(jj-first+1));
                            s_corr_i568_peri(jj-first+1) = std((slice(msk_peri)-i568_immode(jj-first+1))/i568_immean(jj-first+1));
                        end
                    end
                end
            end
            save(analysis_file, 'l', 'a', 'p', 'lh', 'lv', 'lellipseh', 'lellipsev', 'cut_t', 'first', 'last', 't_res', 'i488_peri', 'i488_center', 'i488_immean', 'i488_immode', 'i568_peri', 'i568_center', 'i568_immean', 'i568_immode','m_corr_i488_peri','s_corr_i488_peri','m_corr_i568_peri','s_corr_i568_peri');
            
            if ii < numel(dirs2)
                waitbar((ii)/numel(dirs2), hwb, cat(2, 'Analyzing movie ', cat(2, cat(2, num2str(ii+1), filesep), cat(2, num2str(numel(dirs2)), ' ...'))));
            else
                waitbar(1, hwb, 'Done!');
            end
        end
    end
    
    cd(pwd);
end

%% Run plots.
dirs1 = dirs1_orig;
dirs2 = dirs2_orig;

if plot_flag
    a_abs_1 = nan.*zeros(numel(dirs1), ntp);
    rate_absarea_contraction_1 = nan.*zeros(numel(dirs1), ntp);
    t_absareacontraction_1 = nan.*zeros(1, numel(dirs1));
    p_abs_1 = nan.*zeros(numel(dirs1), ntp);
    rate_absperim_contraction_1 = nan.*zeros(numel(dirs1), ntp);
    t_absperimcontraction_1 = nan.*zeros(1, numel(dirs1));
    a_norm_1 = nan.*zeros(numel(dirs1), ntp);
    p_norm_1 = nan.*zeros(numel(dirs1), ntp);
    aplength_1 = nan.*zeros(numel(dirs1), ntp);
    rate_aplength_1 = nan.*zeros(numel(dirs1), ntp);
    aplengthellipse_1 = nan.*zeros(numel(dirs1), ntp);
    rate_aplengthellipse_1 = nan.*zeros(numel(dirs1), ntp);
    dvlength_1 = nan.*zeros(numel(dirs1), ntp);
    rate_dvlength_1 = nan.*zeros(numel(dirs1), ntp);
    dvlengthellipse_1 = nan.*zeros(numel(dirs1), ntp);
    rate_dvlengthellipse_1 = nan.*zeros(numel(dirs1), ntp);
    anisotropy_1 = nan.*zeros(numel(dirs1), ntp);
    anisotropy_norm_1 = nan.*zeros(numel(dirs1), ntp);
    anisotropyellipse_norm_1 = nan.*zeros(numel(dirs1), ntp);
    rate_anisotropy_1 = nan.*zeros(numel(dirs1), ntp);
    anisotropyellipse_1 = nan.*zeros(numel(dirs1), ntp);
    rate_anisotropyellipse_1 = nan.*zeros(numel(dirs1), ntp);
    i488peri_norm_1 = nan.*zeros(numel(dirs1), ntp);
    rate_i488perinorm_1 = nan.*zeros(numel(dirs1), ntp);
    i568peri_norm_1 = nan.*zeros(numel(dirs1), ntp);
    rate_i568perinorm_1 = nan.*zeros(numel(dirs1), ntp);
    i488center_norm_1 = nan.*zeros(numel(dirs1), ntp);
    rate_i488centernorm_1 = nan.*zeros(numel(dirs1), ntp);
    i568center_norm_1 = nan.*zeros(numel(dirs1), ntp);
    rate_i568centernorm_1 = nan.*zeros(numel(dirs1), ntp);
    total_i568centernorm_1 = nan.*zeros(numel(dirs1), ntp);
    rate_total_i568centernorm_1 = nan.*zeros(numel(dirs1), ntp);
    total_i568perinorm_1 = nan.*zeros(numel(dirs1), ntp);
    total_i488centernorm_1 = nan.*zeros(numel(dirs1), ntp);
    total_i488perinorm_1 = nan.*zeros(numel(dirs1), ntp);
    
    i488_immean_all_1 = nan.*zeros(numel(dirs1), ntp);
    i568_immean_all_1 = nan.*zeros(numel(dirs1), ntp);
    i488_immode_all_1 = nan.*zeros(numel(dirs1), ntp);
    i568_immode_all_1 = nan.*zeros(numel(dirs1), ntp);
    
    m_i488_corr_all_1 = nan.*zeros(numel(dirs1), ntp);
    s_i488_corr_all_1 = nan.*zeros(numel(dirs1), ntp);
    m_i568_corr_all_1 = nan.*zeros(numel(dirs1), ntp);
    s_i568_corr_all_1 = nan.*zeros(numel(dirs1), ntp);
    
    for ii = 1:numel(dirs1)
        if exist(dirs1{ii}, 'file')
            load(strrep(dirs1{ii}, 'annotations', 'analysis'));
        elseif exist(cat(2, dirs1{ii}, 'analysis_healing.mat'), 'file')
            load(cat(2, dirs1{ii}, 'analysis_healing'));
        else
            disp(cat(2, 'Error when loading ', dirs1{ii}));
        end
        
        startr = cut_t - the_cut_t + 1;
        if startr < 1
            startr = 1;
        end
        endr = min(ntp, numel(a));
        if startr > endr
            startr = endr-numel(a)+1;
        end
        
        endw = endr - startr + 1;
        
        a_abs_1(ii, 1:endw) = power(xyres, 2) .* a(startr:endr); % Convert from pixels to square microns.
        a_norm_1(ii, 1:endw) = 100 .* a(startr:endr) ./ mean(a(startr:cut_t)); % Percentage of the initial value.
        a_max_norm_1(ii,1:endw) = 100 .* a(startr:endr) ./ max(a(startr:endr)); % Percentage of the maximum wound area.
        
        % save the areas of the movies....................................................
        if save_area==1
            area_absolute=a_abs_1(ii, 1:endw);
            area_normalized=a_norm_1(ii, 1:endw);
            indsep = strfind(dirs1{ii}, filesep);
            dirs1{ii} = dirs1{ii}(1:indsep(end));
            cd(dirs1{ii})
            if (strcmp(extension, 'mat'))
                area='area_curve.mat';
                save(area, 'area_absolute','area_normalized');
            elseif (strcmp(extension, 'csv'))
                area='area_curve.csv';
                savecsv(area, 'area_absolute',area_absolute,' area_normalized', area_normalized);
            end
        end
        %..............................................................................................
        a_0 = circshift(a_abs_1(ii, 1:endw), [0 (60/t_res)]);
        % contraction rate in um2/min
        rate_absarea_contraction_1(ii, 1:endw) = a_abs_1(ii, 1:endw) - a_0;
        rate_absarea_contraction_1(ii, 1:(60/t_res)) = nan;
        
        i488_immean_all_1(ii, 1:endw) = i488_immean(startr:endr);
        i568_immean_all_1(ii, 1:endw) = i568_immean(startr:endr);
        i488_immode_all_1(ii, 1:endw) = i488_immode(startr:endr);
        i568_immode_all_1(ii, 1:endw) = i568_immode(startr:endr);
        m_i488_corr_all_1(ii, 1:endw) = m_corr_i488_peri(startr:endr);
        s_i488_corr_all_1(ii, 1:endw) = s_corr_i488_peri(startr:endr);
        m_i568_corr_all_1(ii, 1:endw) = m_corr_i568_peri(startr:endr);
        s_i568_corr_all_1(ii, 1:endw) = s_corr_i568_peri(startr:endr);
        
        
        
        % AP and DV lengths and their rates of change.
        if exist('lh', 'var')
            aplength_1(ii, 1:endw) = xyres .* lh(startr:endr); % convert horizontal length to um
            dvlength_1(ii, 1:endw) = xyres .* lv(startr:endr); % convert vertical length to um
            anisotropy_1(ii, 1:endw) = lh(startr:endr) ./ lv(startr:endr);
            anisotropy_norm_1(ii, 1:endw) = 100 .* (lh(startr:endr) ./ lv(startr:endr)) ./ mean(lh(startr:cut_t)./lv(startr:cut_t));
            
            ap_0 = circshift(aplength_1(ii, 1:endw), [0 (60/t_res)]);
            dv_0 = circshift(dvlength_1(ii, 1:endw), [0 (60/t_res)]);
            an_0 = circshift(anisotropy_1(ii, 1:endw), [0 (60/t_res)]);
            
            % rate of length/aniosotropy change per min
            rate_aplength_1(ii, 1:endw) = aplength_1(ii, 1:endw) - ap_0;
            rate_aplength_1(ii, 1:(60/t_res)) = nan;
            rate_dvlength_1(ii, 1:endw) = dvlength_1(ii, 1:endw) - dv_0;
            rate_dvlength_1(ii, 1:(60/t_res)) = nan;
            rate_anisotropy_1(ii, 1:endw) = anisotropy_1(ii, 1:endw) - an_0;
            rate_anisotropy_1(ii, 1:(60/t_res)) = nan;
        end
        
        % Compute the time point where the contraction rate gets closer to zero while being still positive (the wound has not started to close).
        vect = sign(smooth_curve2(rate_absarea_contraction_1(ii, :), 3));
        vect = circshift(vect, [0 1]) + vect;
        t_min = find(vect < eps);
        t_min = t_min(t_min>the_cut_t);
        if ~ isempty(t_min)
            t_min = t_min(1)-1;
        else
            t_min = 1;
        end
        t_absareacontraction_1(ii) = t(t_min);
        
        % Computer perimeter
        p_abs_1(ii, 1:endw) = xyres .* p(startr:endr); % Convert to microns.
        p_norm_1(ii, 1:endw) = 100 .* p(startr:endr) ./ mean(p(startr:cut_t)); % Normalize to the percentage of the initial perimeter.
        p_0 = circshift(p_abs_1(ii, 1:endw), [0 (60/t_res)]);
        rate_absperim_contraction_1(ii, 1:endw) = p_abs_1(ii, 1:endw) - p_0;
        rate_absperim_contraction_1(ii, 1:(60/t_res)) = nan;
        
        % Compute the time point where the contraction rate gets closer to zero while being still positive (the wound has not started to close).
        vect = sign(smooth_curve2(rate_absperim_contraction_1(ii, :), 3));
        vect = circshift(vect, [0 1]) + vect;
        t_min = find(vect < eps);
        t_min = t_min(t_min>the_cut_t);
        if ~ isempty(t_min)
            t_min = t_min(1)-1;
        else
            t_min = 1;
        end
        t_absperimcontraction_1(ii) = t(t_min);
        
        % Perimeter intensity.
        switch normalize_intensity_flag
            case 2
                i488peri_norm_1(ii, 1:endw) = ((i488_peri(startr:endr)) - i488_immode(startr:endr)) ./ i488_immean(startr:endr); % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                %                 i488peri_norm_1(ii, 1:endw) = i488peri_norm_1(ii, 1:endw)./repmat(rmsnan(i488peri_norm_1(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
            case 1
                i488peri_norm_1(ii, 1:endw) = (i488_peri(startr:endr)) ./ i488_immean(startr:endr); % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work.
            case 0
                i488peri_norm_1(ii, 1:endw) = (i488_peri(startr:endr));
        end
        p_0 = circshift(i488peri_norm_1(ii, 1:endw), [0 (60/t_res)]);
        rate_i488perinorm_1(ii, 1:endw) = i488peri_norm_1(ii, 1:endw) - p_0;
        rate_i488perinorm_1(ii, 1:(60/t_res)) = nan;
        total_i488perinorm_1(ii, 1:endw) = i488peri_norm_1(ii, 1:endw) .* p_abs_1(ii, 1:endw) ./ (xyres);
        
        switch normalize_intensity_flag
            case 2
                i568peri_norm_1(ii, 1:endw) = ((i568_peri(startr:endr)) - i568_immode(startr:endr)) ./ i568_immean(startr:endr);  % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                %                 i568peri_norm_1(ii, 1:endw) = i568peri_norm_1(ii, 1:endw)./repmat(rmsnan(i568peri_norm_1(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
            case 1
                i568peri_norm_1(ii, 1:endw) =  (i568_peri(startr:endr)) ./ i568_immean(startr:endr);  % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work (kind off)
            case 0
                i568peri_norm_1(ii, 1:endw) =  (i568_peri(startr:endr));
        end
        p_0 = circshift(i568peri_norm_1(ii, 1:endw), [0 (60/t_res)]);
        rate_i568perinorm_1(ii, 1:endw) = i568peri_norm_1(ii, 1:endw) - p_0;
        rate_i568perinorm_1(ii, 1:(60/t_res)) = nan;
        total_i568perinorm_1(ii, 1:endw) = i568peri_norm_1(ii, 1:endw) .* p_abs_1(ii, 1:endw) ./ xyres;
        
        % Central intensity.
        switch normalize_intensity_flag
            case 2
                i488center_norm_1(ii, 1:endw) = ((i488_center(startr:endr)) - i488_immode(startr:endr)) ./ i488_immean(startr:endr);  % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                %                 i488center_norm_1(ii, 1:endw) = i488center_norm_1(ii, 1:endw)./repmat(rmsnan(i488center_norm_1(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
            case 1
                i488center_norm_1(ii, 1:endw) = (i488_center(startr:endr)) ./ i488_immean(startr:endr); % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work.
            case 0
                i488center_norm_1(ii, 1:endw) = (i488_center(startr:endr));
        end
        p_0 = circshift(i488center_norm_1(ii, 1:endw), [0 (60/t_res)]);
        rate_i488centernorm_1(ii, 1:endw) = i488center_norm_1(ii, 1:endw) - p_0;
        rate_i488centernorm_1(ii, 1:(60/t_res)) = nan;
        total_i488centernorm_1(ii, 1:endw) = i488center_norm_1(ii, 1:endw) .* a_abs_1(ii, 1:endw) ./ (xyres * xyres);
        switch normalize_intensity_flag
            case 2
                i568center_norm_1(ii, 1:endw) = ((i568_center(startr:endr)) - i568_immode(startr:endr)) ./ i568_immean(startr:endr); % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                %                 i568center_norm_1(ii, 1:endw) = i568center_norm_1(ii, 1:endw)./repmat(rmsnan(i568center_norm_1(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
            case 1
                i568center_norm_1(ii, 1:endw) =  (i568_center(startr:endr)) ./ i568_immean(startr:endr);  % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work.
            case 0
                i568center_norm_1(ii, 1:endw) =  (i568_center(startr:endr));
        end
        p_0 = circshift(i568center_norm_1(ii, 1:endw), [0 (60/t_res)]);
        rate_i568centernorm_1(ii, 1:endw) = i568center_norm_1(ii, 1:endw) - p_0;
        rate_i568centernorm_1(ii, 1:(60/t_res)) = nan;
        
        total_i568centernorm_1(ii, 1:endw) = i568center_norm_1(ii, 1:endw) .* a_abs_1(ii, 1:endw) ./ (xyres * xyres);
        p_0 = circshift(total_i568centernorm_1(ii, 1:endw), [0 (60/t_res)]);
        rate_total_i568centernorm_1(ii, 1:endw) = total_i568centernorm_1(ii, 1:endw) - p_0;
        rate_total_i568centernorm_1(ii, 1:(60/t_res)) = nan;
    end
    
    [m_aabs_1, s_aabs_1, e_aabs_1] = rmsnan(a_abs_1);
    [m_rateabsarea_1, s_rateabsarea_1, e_rateabsarea_1] = rmsnan(rate_absarea_contraction_1);
    [m_pabs_1, s_pabs_1, e_pabs_1] = rmsnan(p_abs_1);
    [m_rateabsperim_1, s_rateabsperim_1, e_rateabsperim_1] = rmsnan(rate_absperim_contraction_1);
    [m_anorm_1, s_anorm_1, e_anorm_1] = rmsnan(a_norm_1);
    [m_max_anorm_1, s_max_anorm_1, e_max_anorm_1] = rmsnan(a_max_norm_1);
    [m_pnorm_1, s_pnorm_1, e_pnorm_1] = rmsnan(p_norm_1);
    [m_anisotropy_1, s_anisotropy_1, e_anisotropy_1] = rmsnan(anisotropy_1);
    [m_anisotropynorm_1, s_anisotropynorm_1, e_anisotropynorm_1] = rmsnan(anisotropy_norm_1);
    [m_anisotropyellipse_1, s_anisotropyellipse_1, e_anisotropyellipse_1] = rmsnan(anisotropyellipse_1);
    [m_anisotropyellipsenorm_1, s_anisotropyellipsenorm_1, e_anisotropyellipsenorm_1] = rmsnan(anisotropyellipse_norm_1);
    [m_i488perinorm_1, s_i488perinorm_1, e_i488perinorm_1] = rmsnan(i488peri_norm_1);
    [m_rate488perinorm_1, s_rate488perinorm_1, e_rate488perinorm_1] = rmsnan(rate_i488perinorm_1);
    [m_i568perinorm_1, s_i568perinorm_1, e_i568perinorm_1] = rmsnan(i568peri_norm_1);
    [m_rate568perinorm_1, s_rate568perinorm_1, e_rate568perinorm_1] = rmsnan(rate_i568perinorm_1);
    [m_i488centernorm_1, s_i488centernorm_1, e_i488centernorm_1] = rmsnan(i488center_norm_1);
    [m_rate488centernorm_1, s_rate488centernorm_1, e_rate488centernorm_1] = rmsnan(rate_i488centernorm_1);
    [m_i568centernorm_1, s_i568centernorm_1, e_i568centernorm_1] = rmsnan(i568center_norm_1);
    [m_rate568centernorm_1, s_rate568centernorm_1, e_rate568centernorm_1] = rmsnan(rate_i568centernorm_1);
    
    if ~ isempty(dirs2)
        clear lh lellipseh; % This is necessary to prevent use of the values for dirs1.
        
        a_abs_2 = nan.*zeros(numel(dirs2), ntp);
        rate_absarea_contraction_2 = nan.*zeros(numel(dirs2), ntp);
        t_absareacontraction_2 = nan.*zeros(1, numel(dirs2));
        p_abs_2 = nan.*zeros(numel(dirs2), ntp);
        rate_absperim_contraction_2 = nan.*zeros(numel(dirs2), ntp);
        t_absperimcontraction_2 = nan.*zeros(1, numel(dirs2));
        a_norm_2 = nan.*zeros(numel(dirs2), ntp);
        p_norm_2 = nan.*zeros(numel(dirs2), ntp);
        aplength_2 = nan.*zeros(numel(dirs2), ntp);
        rate_aplength_2 = nan.*zeros(numel(dirs2), ntp);
        aplengthellipse_2 = nan.*zeros(numel(dirs2), ntp);
        rate_aplengthellipse_2 = nan.*zeros(numel(dirs2), ntp);
        dvlength_2 = nan.*zeros(numel(dirs2), ntp);
        rate_dvlength_2 = nan.*zeros(numel(dirs2), ntp);
        dvlengthellipse_2 = nan.*zeros(numel(dirs2), ntp);
        rate_dvlengthellipse_2 = nan.*zeros(numel(dirs2), ntp);
        anisotropy_2 = nan.*zeros(numel(dirs2), ntp);
        anisotropy_norm_2 = nan.*zeros(numel(dirs2), ntp);
        anisotropyellipse_norm_2 = nan.*zeros(numel(dirs2), ntp);rate_anisotropy_2 = nan.*zeros(numel(dirs2), ntp);
        anisotropyellipse_2 = nan.*zeros(numel(dirs2), ntp);
        rate_anisotropyellipse_2 = nan.*zeros(numel(dirs2), ntp);
        i488peri_norm_2 = nan.*zeros(numel(dirs2), ntp);
        rate_i488perinorm_2 = nan.*zeros(numel(dirs2), ntp);
        i568peri_norm_2 = nan.*zeros(numel(dirs2), ntp);
        rate_i568perinorm_2 = nan.*zeros(numel(dirs2), ntp);
        i488center_norm_2 = nan.*zeros(numel(dirs2), ntp);
        rate_i488centernorm_2 = nan.*zeros(numel(dirs2), ntp);
        i568center_norm_2 = nan.*zeros(numel(dirs2), ntp);
        rate_i568centernorm_2 = nan.*zeros(numel(dirs2), ntp);
        total_i568centernorm_2 = nan.*zeros(numel(dirs2), ntp);
        rate_total_i568centernorm_2 = nan.*zeros(numel(dirs2), ntp);
        total_i568perinorm_2 = nan.*zeros(numel(dirs2), ntp);
        total_i488centernorm_2 = nan.*zeros(numel(dirs2), ntp);
        total_i488perinorm_2 = nan.*zeros(numel(dirs2), ntp);
        
        i488_immean_all_2 = nan.*zeros(numel(dirs2), ntp);
        i568_immean_all_2 = nan.*zeros(numel(dirs2), ntp);
        i488_immode_all_2 = nan.*zeros(numel(dirs2), ntp);
        i568_immode_all_2 = nan.*zeros(numel(dirs2), ntp);
        
        m_i488_corr_all_2 = nan.*zeros(numel(dirs2), ntp);
        s_i488_corr_all_2 = nan.*zeros(numel(dirs2), ntp);
        m_i568_corr_all_2 = nan.*zeros(numel(dirs2), ntp);
        s_i568_corr_all_2 = nan.*zeros(numel(dirs2), ntp);
        
        
        for ii = 1:numel(dirs2)
            if exist(dirs2{ii}, 'file')
                load(strrep(dirs2{ii}, 'annotations', 'analysis'));
            elseif exist(cat(2, dirs2{ii}, 'analysis_healing.mat'), 'file')
                load(cat(2, dirs2{ii}, 'analysis_healing'));
            else
                disp(cat(2, 'Error when loading ', dirs2{ii}));
            end
            
            startr = cut_t - the_cut_t + 1;
            if startr < 1
                startr = 1;
            end
            endr = min(ntp, numel(a));
            if startr > endr
                startr = endr-numel(a)+1;
            end
            endw = endr - startr + 1;
            
            a_abs_2(ii, 1:endw) = power(xyres, 2) .* a(startr:endr); % Convert from pixels to square microns.
            a_norm_2(ii, 1:endw) = 100 .* a(startr:endr) ./ mean(a(startr:cut_t)); % Percentage of the initial value.
            a_max_norm_2(ii,1:endw) = 100 .* a(startr:endr) ./ max(a(startr:endr)); % Percentage of the maximum wound area.
            
            if save_area==1
                area_absolute=a_abs_2(ii, 1:endw);
                area_normalized=a_norm_2(ii, 1:endw);
                indsep = strfind(dirs2{ii}, filesep);
                dirs2{ii} = dirs2{ii}(1:indsep(end));
                cd(dirs2{ii})
                if (strcmp(extension, 'mat'))
                    area='area_curve.mat';
                    save(area, 'area_absolute','area_normalized');
                elseif (strcmp(extension, 'csv'))
                    area='area_curve.csv';
                    savecsv(area, 'area_absolute',area_absolute,' area_normalized', area_normalized);
                end
            end
            %..............................................................................................
            a_0 = circshift(a_abs_2(ii, 1:endw), [0 (60/t_res)]);
            rate_absarea_contraction_2(ii, 1:endw) = a_abs_2(ii, 1:endw) - a_0;
            rate_absarea_contraction_2(ii, 1:(60/t_res)) = nan;
            
            i488_immean_all_2(ii, 1:endw) = i488_immean(startr:endr);
            i568_immean_all_2(ii, 1:endw) = i568_immean(startr:endr);
            i488_immode_all_2(ii, 1:endw) = i488_immode(startr:endr);
            i568_immode_all_2(ii, 1:endw) = i568_immode(startr:endr);
            
            m_i488_corr_all_2(ii, 1:endw) = m_corr_i488_peri(startr:endr);
            s_i488_corr_all_2(ii, 1:endw) = s_corr_i488_peri(startr:endr);
            m_i568_corr_all_2(ii, 1:endw) = m_corr_i568_peri(startr:endr);
            s_i568_corr_all_2(ii, 1:endw) = s_corr_i568_peri(startr:endr);
            
            % AP and DV lengths and their rates of change.
            if exist('lh', 'var')
                aplength_2(ii, 1:endw) = xyres .* lh(startr:endr);
                dvlength_2(ii, 1:endw) = xyres .* lv(startr:endr);
                anisotropy_2(ii, 1:endw) = lh(startr:endr) ./ lv(startr:endr);
                anisotropy_norm_2(ii, 1:endw) = 100 .* (lh(startr:endr) ./ lv(startr:endr)) ./ mean(lh(startr:cut_t)./lv(startr:cut_t));
                
                ap_0 = circshift(aplength_2(ii, 1:endw), [0 (60/t_res)]);
                dv_0 = circshift(dvlength_2(ii, 1:endw), [0 (60/t_res)]);
                an_0 = circshift(anisotropy_2(ii, 1:endw), [0 (60/t_res)]);
                
                rate_aplength_2(ii, 1:endw) = aplength_2(ii, 1:endw) - ap_0;
                rate_aplength_2(ii, 1:(60/t_res)) = nan;
                rate_dvlength_2(ii, 1:endw) = dvlength_2(ii, 1:endw) - dv_0;
                rate_dvlength_2(ii, 1:(60/t_res)) = nan;
                rate_anisotropy_2(ii, 1:endw) = anisotropy_2(ii, 1:endw) - an_0;
                rate_anisotropy_2(ii, 1:(60/t_res)) = nan;
            end
            
            % AP and DV lengths and their rates of change (based on ellipses).
            if exist('lellipseh', 'var')
                aplengthellipse_2(ii, 1:endw) = xyres .* lellipseh(startr:endr);
                dvlengthellipse_2(ii, 1:endw) = xyres .* lellipsev(startr:endr);
                anisotropyellipse_2(ii, 1:endw) = lellipseh(startr:endr) ./ lellipsev(startr:endr);
                anisotropyellipse_norm_2(ii, 1:endw) = 100 .* (lellipseh(startr:endr) ./ lellipsev(startr:endr)) ./ mean(lellipseh(startr:cut_t)./lellipsev(startr:cut_t));
                
                ap_0 = circshift(aplengthellipse_2(ii, 1:endw), [0 (60/t_res)]);
                dv_0 = circshift(dvlengthellipse_2(ii, 1:endw), [0 (60/t_res)]);
                an_0 = circshift(anisotropyellipse_2(ii, 1:endw), [0 (60/t_res)]);
                
                rate_aplengthellipse_2(ii, 1:endw) = aplengthellipse_2(ii, 1:endw) - ap_0;
                rate_aplengthellipse_2(ii, 1:(60/t_res)) = nan;
                rate_dvlengthellipse_2(ii, 1:endw) = dvlengthellipse_2(ii, 1:endw) - dv_0;
                rate_dvlengthellipse_2(ii, 1:(60/t_res)) = nan;
                rate_anisotropyellipse_2(ii, 1:endw) = anisotropyellipse_2(ii, 1:endw) - an_0;
                rate_anisotropyellipse_2(ii, 1:(60/t_res)) = nan;
            end
            
            % Compute the time point where the contraction rate gets closer to zero while being still positive (the wound has not started to close).
            vect = sign(smooth_curve2(rate_absarea_contraction_2(ii, :), 3));
            vect = circshift(vect, [0 1]) + vect;
            t_min = find(vect < eps);
            t_min = t_min(t_min>the_cut_t);
            if ~ isempty(t_min)
                t_min = t_min(1)-1;
            else
                t_min = 1;
            end
            t_absareacontraction_2(ii) = t(t_min);
            
            
            p_abs_2(ii, 1:endw) = xyres .* p(startr:endr); % Convert to microns.
            p_norm_2(ii, 1:endw) = 100 .* p(startr:endr) ./ mean(p(startr:cut_t)); % Normalize to the percentage of the initial perimeter.
            p_0 = circshift(p_abs_2(ii, 1:endw), [0 (60/t_res)]);
            rate_absperim_contraction_2(ii, 1:endw) = p_abs_2(ii, 1:endw) - p_0;
            rate_absperim_contraction_2(ii, 1:(60/t_res)) = nan;
            
            % Compute the time point where the contraction rate gets closer to zero while being still positive (the wound has not started to close).
            vect = sign(smooth_curve2(rate_absperim_contraction_2(ii, :), 3));
            vect = circshift(vect, [0 1]) + vect;
            t_min = find(vect < eps);
            t_min = t_min(t_min>the_cut_t);
            if ~ isempty(t_min)
                t_min = t_min(1)-1;
            else
                t_min = 1;
            end
            t_absperimcontraction_2(ii) = t(t_min);
            
            % Perimeter intensity.
            switch normalize_intensity_flag
                case 2
                    i488peri_norm_2(ii, 1:endw) = ((i488_peri(startr:endr)) - i488_immode(startr:endr)) ./ i488_immean(startr:endr);  % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                    %                     i488peri_norm_2(ii, 1:endw) = i488peri_norm_2(ii, 1:endw)./repmat(rmsnan(i488peri_norm_2(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
                case 1
                    i488peri_norm_2(ii, 1:endw) = (i488_peri(startr:endr)) ./ i488_immean(startr:endr); % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work.
                case 0
                    i488peri_norm_2(ii, 1:endw) = (i488_peri(startr:endr));
            end
            p_0 = circshift(i488peri_norm_2(ii, 1:endw), [0 (60/t_res)]);
            rate_i488perinorm_2(ii, 1:endw) = i488peri_norm_2(ii, 1:endw) - p_0;
            rate_i488perinorm_2(ii, 1:(60/t_res)) = nan;
            total_i488perinorm_2(ii, 1:endw) = i488peri_norm_2(ii, 1:endw) .* p_abs_2(ii, 1:endw) ./ (xyres);
            
            %*   volver
            switch normalize_intensity_flag
                case 2
                    i568peri_norm_2(ii, 1:endw) = ((i568_peri(startr:endr)) - i568_immode(startr:endr)) ./ i568_immean(startr:endr);  % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                    %                     i568peri_norm_2(ii, 1:endw) = i568peri_norm_2(ii, 1:endw)./repmat(rmsnan(i568peri_norm_2(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
                case 1
                    i568peri_norm_2(ii, 1:endw) =  (i568_peri(startr:endr)) ./ i568_immean(startr:endr);  % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work (kind off)
                case 0
                    i568peri_norm_2(ii, 1:endw) =  (i568_peri(startr:endr));
            end
            p_0 = circshift(i568peri_norm_2(ii, 1:endw), [0 (60/t_res)]);
            rate_i568perinorm_2(ii, 1:endw) = i568peri_norm_2(ii, 1:endw) - p_0;
            rate_i568perinorm_2(ii, 1:(60/t_res)) = nan;
            total_i568perinorm_2(ii, 1:endw) = i568peri_norm_2(ii, 1:endw) .* p_abs_2(ii, 1:endw) ./ xyres;
            
            % Center intensity.
            switch normalize_intensity_flag
                case 2
                    i488center_norm_2(ii, 1:endw) = ((i488_center(startr:endr)) - i488_immode(startr:endr)) ./ i488_immean(startr:endr);  % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                    %                     i488center_norm_2(ii, 1:endw) = i488center_norm_2(ii, 1:endw)./repmat(rmsnan(i488center_norm_2(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
                case 1
                    i488center_norm_2(ii, 1:endw) = (i488_center(startr:endr)) ./ i488_immean(startr:endr); % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work.
                case 0
                    i488center_norm_2(ii, 1:endw) = (i488_center(startr:endr));
            end
            p_0 = circshift(i488center_norm_2(ii, 1:endw), [0 (60/t_res)]);
            rate_i488centernorm_2(ii, 1:endw) = i488center_norm_2(ii, 1:endw) - p_0;
            rate_i488centernorm_2(ii, 1:(60/t_res)) = nan;
            total_i488centernorm_2(ii, 1:endw) = i488center_norm_2(ii, 1:endw) .* a_abs_2(ii, 1:endw) ./ (xyres * xyres);
            
            switch normalize_intensity_flag
                case 2
                    i568center_norm_2(ii, 1:endw) = ((i568_center(startr:endr)) - i568_immode(startr:endr)) ./ i568_immean(startr:endr);  % This division by im_mean does not do anything, as in the next line you will divide by average value of this same vector, effectively multiplying again by im_mean.
                    %                     i568center_norm_2(ii, 1:endw) = i568center_norm_2(ii, 1:endw)./repmat(rmsnan(i568center_norm_2(ii, the_cut_t+1:end)')', [1, numel(1:endw)]); % Normalize to the average intensity after wounding to be able to compare different wounds.
                case 1
                    i568center_norm_2(ii, 1:endw) =  (i568_center(startr:endr)) ./ i568_immean(startr:endr);  % Dividing by the mean does not work if the image is registered and there are black areas as a result of the shift. The mode, on the other hand, should still work.
                case 0
                    i568center_norm_2(ii, 1:endw) =  (i568_center(startr:endr));
            end
            p_0 = circshift(i568center_norm_2(ii, 1:endw), [0 (60/t_res)]);
            rate_i568centernorm_2(ii, 1:endw) = i568center_norm_2(ii, 1:endw) - p_0;
            rate_i568centernorm_2(ii, 1:(60/t_res)) = nan;
            
            total_i568centernorm_2(ii, 1:endw) = i568center_norm_2(ii, 1:endw) .* a_abs_2(ii, 1:endw) ./ (xyres * xyres);
            p_0 = circshift(total_i568centernorm_2(ii, 1:endw), [0 (60/t_res)]);
            rate_total_i568centernorm_2(ii, 1:endw) = total_i568centernorm_2(ii, 1:endw) - p_0;
            rate_total_i568centernorm_2(ii, 1:(60/t_res)) = nan;
        end
        
        [m_aabs_2, s_aabs_2, e_aabs_2] = rmsnan(a_abs_2);
        [m_rateabsarea_2, s_rateabsarea_2, e_rateabsarea_2] = rmsnan(rate_absarea_contraction_2);
        [m_pabs_2, s_pabs_2, e_pabs_2] = rmsnan(p_abs_2);
        [m_rateabsperim_2, s_rateabsperim_2, e_rateabsperim_2] = rmsnan(rate_absperim_contraction_2);
        [m_anorm_2, s_anorm_2, e_anorm_2] = rmsnan(a_norm_2);
        [m_max_anorm_2,s_max_anorm_2,e_max_anorm_2] = rmsnan(a_max_norm_2);
        [m_pnorm_2, s_pnorm_2, e_pnorm_2] = rmsnan(p_norm_2);                        [m_anisotropy_2, s_anisotropy_2, e_anisotropy_2] = rmsnan(anisotropy_2);
        [m_anisotropynorm_2, s_anisotropynorm_2, e_anisotropynorm_2] = rmsnan(anisotropy_norm_2);
        [m_anisotropyellipse_2, s_anisotropyellipse_2, e_anisotropyellipse_2] = rmsnan(anisotropyellipse_2);
        [m_anisotropyellipsenorm_2, s_anisotropyellipsenorm_2, e_anisotropyellipsenorm_2] = rmsnan(anisotropyellipse_norm_2);
        [m_i488perinorm_2, s_i488perinorm_2, e_i488perinorm_2] = rmsnan(i488peri_norm_2);
        [m_rate488perinorm_2, s_rate488perinorm_2, e_rate488perinorm_2] = rmsnan(rate_i488perinorm_2);
        [m_i568perinorm_2, s_i568perinorm_2, e_i568perinorm_2] = rmsnan(i568peri_norm_2);
        [m_rate568perinorm_2, s_rate568perinorm_2, e_rate568perinorm_2] = rmsnan(rate_i568perinorm_2);
        [m_i488centernorm_2, s_i488centernorm_2, e_i488centernorm_2] = rmsnan(i488center_norm_2);
        [m_rate488centernorm_2, s_rate488centernorm_2, e_rate488centernorm_2] = rmsnan(rate_i488centernorm_2);
        [m_i568centernorm_2, s_i568centernorm_2, e_i568centernorm_2] = rmsnan(i568center_norm_2);
        [m_rate568centernorm_2, s_rate568centernorm_2, e_rate568centernorm_2] = rmsnan(rate_i568centernorm_2);
    end
    
    % PLOTS.
    
    % AREA
    % abs area vs time, each experimental group -----------------------------------
    cm1 = jet(numel(dirs1));
    cm2 = jet(numel(dirs2));
    
    f1=figure; hold on;
    for ii = 1:numel(dirs1)
        plot(t, a_abs_1(ii, :), 'Color', cm1(ii, :), 'LineWidth', 4);
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('area (\mum^2)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:100:y_range(2));
    
    if save_flag
        saveas(f1, cat(2, 'time_area_', legend_str{1}), 'fig');
    end
    
    if ~isempty(dirs2)
        f2=figure; hold on;
        for ii = 1:numel(dirs2)
            plot(t, a_abs_2(ii, :), 'Color', cm2(ii, :), 'LineWidth', 4);
            hold on;
        end
        xlim(t_range);
        y_range = ylim;
        ylim([0 y_range(2)]);
        xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
        ylabel('area (\mum^2)', 'FontSize', 32, 'FontWeight', 'bold');
        set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:100:y_range(2));
        
        if save_flag
            saveas(f2, cat(2, 'time_area_', legend_str{2}), 'fig');
        end
    end
    
    % Absolute area vs time all experiments--------------------------------------------------
    
    f3=figure; hold on;
    h = zeros(1:4);
    ind_max = ntp;
    h(1) = errorbar(t(1:ind_max),m_aabs_1,e_aabs_1);
    set(h(1), 'LineStyle', 'none', 'LineWidth', 2, 'Color', [0.6 0.6 1]);
    if ~isempty(dirs2)
        h(3) = errorbar(t(1:ind_max), m_aabs_2, e_aabs_2);
        set(h(3), 'LineStyle', 'none', 'LineWidth', 2, 'Color', [1 0.6 0.6]);
    end
    h(2) = plot(t(1:ind_max), m_aabs_1, 'Color', 'b', 'LineWidth', 4,'DisplayName',legend_str{1});
    if ~isempty(dirs2)
        h(4) = plot(t(1:ind_max), m_aabs_2, 'Color', 'r', 'LineWidth', 4,'DisplayName',legend_str{2});
    end
    hold off;
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('mean area  (\mum^2)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:100:y_range(2));
    
    if save_flag
        saveas(f3, cat(2,'time_meanarea_',legend_str{1},'_',legend_str{2}),'fig')
    end
    
    % norm area vs time, all experiments -----------------------------------
    f4=figure; hold on;
    for ii = 1:numel(dirs1)
        plot(t, a_max_norm_1(ii, :), 'Color', cm1(ii, :), 'LineWidth', 4);
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('area (%)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:50:y_range(2));
    
    if save_flag
        saveas(f4, cat(2, 'time_normarea_', legend_str{1}), 'fig');
    end
    
    if ~isempty(dirs2)
        f5=figure; hold on;
        for ii = 1:numel(dirs2)
            plot(t, a_max_norm_2(ii, :), 'Color', cm2(ii, :), 'LineWidth', 4);
        end
        xlim(t_range);
        y_range = ylim;
        ylim([0 y_range(2)]);
        xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
        ylabel('area (%)', 'FontSize', 32, 'FontWeight', 'bold');
        set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:50:y_range(2));
        
        if save_flag
            saveas(f5, cat(2, 'time_normarea_', legend_str{2}), 'fig');
        end
    end
    
    % mean norm area vs time
    f6=figure; hold on;
    ind_max =  length(t);
    h = errorbar(t(1:ind_max), m_max_anorm_1(1:ind_max), e_max_anorm_1(1:ind_max));
    set(h, 'LineStyle', 'none', 'LineWidth', 2, 'Color', [0.6 0.6 1]);
    plot(t(1:ind_max), m_max_anorm_1(1:ind_max), 'Color', 'b', 'LineWidth', 4);
    if ~isempty(dirs2)
        h = errorbar(t(1:ind_max), m_max_anorm_2(1:ind_max), e_max_anorm_2(1:ind_max));
        set(h, 'LineStyle', 'none', 'LineWidth', 2, 'Color', [1 0.6 0.6]);
        plot(t(1:ind_max), m_max_anorm_2(1:ind_max), 'Color', 'r', 'LineWidth', 4);
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('area (%)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:50:y_range(2));
    
    if save_flag
        saveas(f6, cat(2, 'time_meannormarea_', legend_str{1}, '_', legend_str{2}), 'fig');
    end
    
    % FLUORESCENCE
    % all 488 perimeter fluorescence -----------------------------------------------------------------------
    f13 = figure; hold on;
    i488peri_norm_norm_1 = zeros(numel(dirs1),ntp);
    for ii = 1:numel(dirs1)
%         i488peri_norm_norm_1(ii,:) = 100.*i488peri_norm_1(ii,:)./mean(i488peri_norm_1(ii, 1:the_cut_t));
        i488peri_norm_norm_1(ii,:) = 100 * i488peri_norm_1(ii,:)./i488peri_norm_1(ii,the_cut_t+1);
        plot(t,i488peri_norm_1(ii,:),'Color',cm1(ii,:), 'LineWidth', 4)
        %         plot(t,i488peri_norm_1(ii,:),'Color',cm1(ii,:), 'LineWidth', 4)
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('488 intensity (%)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:100:y_range(2));
    
    if save_flag
        saveas(f13, cat(2, 'time_norm488_', legend_str{1}), 'fig');
    end
    
    if ~isempty(dirs2)
        f14 = figure; hold on;
        i488peri_norm_norm_2 = zeros(numel(dirs2),ntp);
        for ii = 1:numel(dirs2)
%             i488peri_norm_norm_2(ii,:) = 100.*i488peri_norm_2(ii,:)./mean(i488peri_norm_2(ii, 1:the_cut_t));
            i488peri_norm_norm_2(ii,:) = 100 * i488peri_norm_2(ii,:)./i488peri_norm_2(ii,the_cut_t+1);
            plot(t,i488peri_norm_2(ii,:),'Color',cm2(ii,:), 'LineWidth', 4)
            %             plot(t,i488peri_norm_2(ii,:),'Color',cm2(ii,:), 'LineWidth', 4)
        end
        xlim(t_range);
        y_range = ylim;
        ylim([0 y_range(2)]);
        xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
        ylabel('488 intensity (%)', 'FontSize', 32, 'FontWeight', 'bold');
        set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:100:y_range(2));
        
        if save_flag
            saveas(f14, cat(2, 'time_norm488_', legend_str{2}), 'fig');
        end
    end
    
    
    
    
    
    % Average 488 perimeter fluorescence
    f15=figure; hold on;
    h = zeros(1:4);
%         [m_i488peri_norm_norm_1, s_i488peri_norm_norm_1, e_i488peri_norm_norm_1] = rmsnan(i488peri_norm_norm_1);
%         [m_i488peri_norm_norm_2, s_i488peri_norm_norm_2, e_i488peri_norm_norm_2] = rmsnan(i488peri_norm_norm_2);
    [m_i488peri_norm_norm_1, s_i488peri_norm_norm_1, e_i488peri_norm_norm_1] = rmsnan(i488peri_norm_1);
    [m_i488peri_norm_norm_2, s_i488peri_norm_norm_2, e_i488peri_norm_norm_2] = rmsnan(i488peri_norm_2);
    h(1) = errorbar(t(1:ind_max), m_i488peri_norm_norm_1(1:ind_max), e_i488peri_norm_norm_1(1:ind_max));
    set(h(1), 'LineStyle', 'none', 'LineWidth', 2, 'Color', [.6 .6 1]);
    
    if ~isempty(dirs2)
        h(3) = errorbar(t(1:ind_max), m_i488peri_norm_norm_2(1:ind_max), e_i488peri_norm_norm_2(1:ind_max));
        set(h(3), 'LineStyle', 'none', 'LineWidth', 2, 'Color', [1 0.6 0.6]);
    end
    h(2) = plot(t(1:ind_max), m_i488peri_norm_norm_1(1:ind_max), 'Color', 'b', 'LineWidth', 4,'DisplayName',legend_str{1});
    if ~isempty(dirs2)
        h(4) = plot(t(1:ind_max), m_i488peri_norm_norm_2(1:ind_max), 'Color', 'r', 'LineWidth', 4,'DisplayName',legend_str{2});
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    %     legend(h([2 4]))
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('myosin intensity (a.u.)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:0.2:y_range(2));
    
    if save_flag
        saveas(f15, cat(2, 'time_meannorm488_', legend_str{1},'_',legend_str{2}), 'fig');
    end
    
    % all 568 perimeter fluorescence -----------------------------------------------------------------------
    f16 = figure; hold on;
    i568peri_norm_norm_1 = zeros(numel(dirs1),ntp);
    for ii = 1:numel(dirs1)
        i568peri_norm_norm_1(ii,:) = 100.*i568peri_norm_1(ii,:)./mean(i568peri_norm_1(ii, 1:the_cut_t));
        %         plot(t,i568peri_norm_norm_1(ii,:),'Color',cm1(ii,:), 'LineWidth', 4)
        plot(t,i568peri_norm_1(ii,:),'Color',cm1(ii,:), 'LineWidth', 4)
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('568 intensity (%)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:20:y_range(2));
    
    if save_flag
        saveas(f16, cat(2, 'time_norm568_', legend_str{1}), 'fig');
    end
    
    if ~isempty(dirs2)
        f17 = figure; hold on;
        i568peri_norm_norm_2 = zeros(numel(dirs2),ntp);
        for ii = 1:numel(dirs2)
            i568peri_norm_norm_2(ii,:) = 100.*i568peri_norm_2(ii,:)./mean(i568peri_norm_2(ii, 1:the_cut_t));
            %             plot(t,i568peri_norm_norm_2(ii,:),'Color',cm2(ii,:), 'LineWidth', 4)
            plot(t,i568peri_norm_2(ii,:),'Color',cm2(ii,:), 'LineWidth', 4)
        end
        xlim(t_range);
        y_range = ylim;
        ylim([0 y_range(2)]);
        xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
        ylabel('568 intensity (%)', 'FontSize', 32, 'FontWeight', 'bold');
        set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:20:y_range(2));
        
        if save_flag
            saveas(f17, cat(2, 'time_norm568_', legend_str{2}), 'fig');
        end
    end
    
    % Initial 568 fluorescence
    f18a = figure; hold on;
    i568peri_norm_init_1 = i568peri_norm_1(:,1);
    i568peri_norm_init_2 = i568peri_norm_2(:,2);
    boxplotMWstd([0.5 1.5],i568peri_norm_init_1,2,'bo',1)
    boxplotMWstd([2 3],i568peri_norm_init_2,2,'ro',1)
    set(gca, 'Box', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3);
    set(gca, 'XTick', [1 2.5], 'XTickLabel', legend_str);
    ylabel( 'initial 568 intensity (a.u.)', 'FontSize', 32)
    xlim([0 3.5]);
    
    % Average 568 perimeter fluorescence
    f18=figure; hold on;
    %     [m_i568peri_norm_norm_1, s_i568peri_norm_norm_1, e_i568peri_norm_norm_1] = rmsnan(i568peri_norm_norm_1);
    %     [m_i568peri_norm_norm_2, s_i568peri_norm_norm_2, e_i568peri_norm_norm_2] = rmsnan(i568peri_norm_norm_2);
    [m_i568peri_norm_norm_1, s_i568peri_norm_norm_1, e_i568peri_norm_norm_1] = rmsnan(i568peri_norm_1);
    [m_i568peri_norm_norm_2, s_i568peri_norm_norm_2, e_i568peri_norm_norm_2] = rmsnan(i568peri_norm_2);
    h = errorbar(t(1:ind_max), m_i568peri_norm_norm_1(1:ind_max), e_i568peri_norm_norm_1(1:ind_max));
    set(h, 'LineStyle', 'none', 'LineWidth', 2, 'Color', [0.6 0.6 1]);
    plot(t(1:ind_max), m_i568peri_norm_norm_1(1:ind_max), 'Color', 'b', 'LineWidth', 4);
    if ~isempty(dirs2)
        h = errorbar(t(1:ind_max), m_i568peri_norm_norm_2(1:ind_max), e_i568peri_norm_norm_2(1:ind_max));
        set(h, 'LineStyle', 'none', 'LineWidth', 2, 'Color', [1 0.6 0.6]);
        plot(t(1:ind_max), m_i568peri_norm_norm_2(1:ind_max), 'Color', 'r', 'LineWidth', 4);
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('568 intensity (%)', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:20:y_range(2));
    
    if save_flag
        saveas(f18, cat(2, 'time_meannorm568_', legend_str{1},'_',legend_str{2}), 'fig');
    end
  
    % heterogeneity in signals
    f25 = figure; hold on; box on;
    hetero_i488_peri_1 = zeros(numel(dirs1),ntp);
    norm_hetero_i488_peri_1 = zeros(numel(dirs1),ntp);
    for ii = 1:numel(dirs1)
        hetero_i488_peri_1(ii,:) = s_i488_corr_all_1(ii, :)./m_i488_corr_all_1(ii, :);
        norm_hetero_i488_peri_1(ii,:) = hetero_i488_peri_1(ii,:)./hetero_i488_peri_1(ii,1);
        plot(t, hetero_i488_peri_1(ii,:), 'Color', cm1(ii, :), 'LineWidth', 4);
    end
    xlim(t_range);
    y_range = max(abs(ylim));
    ylim([0 y_range]);
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('i488 heterogeneity', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',-y_range:0.4:y_range);
    
    if ~isempty(dirs2)
        f26 = figure; hold on; box on;
        hetero_i488_peri_2 = zeros(numel(dirs2),ntp);
        norm_hetero_i488_peri_2 = zeros(numel(dirs2),ntp);
        for ii = 1:numel(dirs2)
            hetero_i488_peri_2(ii,:) = s_i488_corr_all_2(ii, :)./m_i488_corr_all_2(ii, :);
            norm_hetero_i488_peri_2(ii,:) = hetero_i488_peri_2(ii,:)./hetero_i488_peri_2(ii,1);
            plot(t, hetero_i488_peri_2(ii,:), 'Color', cm2(ii, :), 'LineWidth', 4);
        end
        xlim(t_range);
        y_range = max(abs(ylim));
        ylim([0 y_range]);
        xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
        ylabel('i488 heterogeneity', 'FontSize', 32, 'FontWeight', 'bold');
        set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',-y_range:0.4:y_range);
    end
    
    % Average 488 perimeter heterogeneity
    f27=figure; hold on;
    h = zeros(1:4);
    [m_hetero_i488_peri_1, s_hetero_i488_peri_1, e_hetero_i488_peri_1] = rmsnan(hetero_i488_peri_1);
    [m_hetero_i488_peri_2, s_hetero_i488_peri_2, e_hetero_i488_peri_2] = rmsnan(hetero_i488_peri_2);
    h(1) = errorbar(t(1:ind_max), m_hetero_i488_peri_1(1:ind_max), e_hetero_i488_peri_1(1:ind_max));
    set(h(1), 'LineStyle', 'none', 'LineWidth', 2, 'Color', [.6 .6 1]);
    if ~isempty(dirs2)
        h(3) = errorbar(t(1:ind_max), m_hetero_i488_peri_2(1:ind_max), e_hetero_i488_peri_2(1:ind_max));
        set(h(3), 'LineStyle', 'none', 'LineWidth', 2, 'Color', [1 0.6 0.6]);
    end
    h(2) = plot(t(1:ind_max), m_hetero_i488_peri_1(1:ind_max), 'Color', 'b', 'LineWidth', 4,'DisplayName',legend_str{1});
    if ~isempty(dirs2)
        h(4) = plot(t(1:ind_max), m_hetero_i488_peri_2(1:ind_max), 'Color', 'r', 'LineWidth', 4,'DisplayName',legend_str{2});
    end
    xlim(t_range);
    y_range = ylim;
    ylim([0 y_range(2)]);
    %     legend(h([2 4]))
    xlabel('time (min)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('488 heterogeneity', 'FontSize', 32, 'FontWeight', 'bold');
    set(gca, 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3,'XTick',0:10:t_range(end),'YTick',0:0.5:y_range(2));
    
   
    % Initial 488 fluorescence
    figure; hold on;
    i488peri_norm_init_1 = mean(i488peri_norm_1(:,1:2),2);
    i488peri_norm_init_2 = mean(i488peri_norm_2(:,1:2),2);
    boxplotMWstd([0.5 1.5],i488peri_norm_init_1,3,'bo',3);
    boxplotMWstd([2 3],i488peri_norm_init_2,3,'ro',3);
    set(gca, 'XTick', [1 2.5], 'XTickLabel', legend_str);
    set(gca, 'Box', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3);
    ylabel( 'initial 488 intensity (a.u.)', 'FontSize', 32)
    xlim([0 3.5]);
    
    % Final 488 fluorescence
    figure; hold on;
    i488peri_norm_final_1 = mean(i488peri_norm_1(:,end-4:end),2);
    i488peri_norm_final_2 = mean(i488peri_norm_2(:,end-4:end),2);
    boxplotMWstd([0.5 1.5],i488peri_norm_final_1,3,'bo',3);
    boxplotMWstd([2 3],i488peri_norm_final_2,3,'ro',3);
    set(gca, 'XTick', [1 2.5], 'XTickLabel', legend_str);
    set(gca, 'Box', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3);
    ylabel( 'final 488 intensity (a.u.)', 'FontSize', 32)
    xlim([0 3.5]);

    % Max 488 fluorescence
    figure; hold on;
    i488peri_norm_max_1 = max(i488peri_norm_1,[],2);
    i488peri_norm_max_2 = max(i488peri_norm_2,[],2);
    boxplotMWstd([0.5 1.5],i488peri_norm_max_1,3,'bo',3);
    boxplotMWstd([2 3],i488peri_norm_max_2,3,'ro',3);
    set(gca, 'XTick', [1 2.5], 'XTickLabel', legend_str);
    set(gca, 'Box', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3);
    ylabel( 'max 488 intensity (a.u.)', 'FontSize', 32)
    xlim([0 3.5]);
    
    % Percent change in 488 fluorescence
    figure; hold on;
    i488peri_norm_per_1 = 100*i488peri_norm_final_1./i488peri_norm_init_1;
    i488peri_norm_per_2 = 100*i488peri_norm_final_2./i488peri_norm_init_2;
    boxplotMWstd([0.5 1.5],i488peri_norm_per_1,3,'bo',3)
    boxplotMWstd([2 3],i488peri_norm_per_2,3,'ro',3)
    set(gca, 'XTick', [1 2.5], 'XTickLabel', legend_str);
    set(gca, 'Box', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 3);
    ylabel( 'percent myosin increase (%)', 'FontSize', 32)
    xlim([0 3.5]);
    
    
end
