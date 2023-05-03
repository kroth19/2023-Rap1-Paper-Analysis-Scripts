%% Read in fiducials marking tricell jxns & calculate intensity over time

%% Define annotation files

con = {
    % yw 30 sec
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo7/annotations_fiducials.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201002/Embryo1/annotations_fiducials.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201002/Embryo3/annotations_fiducials.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo1/annotations_fiducials.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo2/annotations_fiducials.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo3/annotations_fiducials.mat',...
    };

exp = {
    % Rap1-DN 30 sec
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo2/annotations_fiducials.mat';
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo3/annotations_fiducials.mat';
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo4/annotations_fiducials.mat';
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo2/annotations_fiducials.mat';
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo3/annotations_fiducials.mat';
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo6/annotations_fiducials.mat';
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo5/annotations_fiducials.mat';
    };
    
dirs1_orig = con;
dirs2_orig = exp;


%% Define necessary parameters

dirs1 = dirs1_orig;
dirs2 = dirs2_orig;
legend_str = {'control','Rap1DN'};
label561 = 'Ecad';
label488 = 'sqh';
channel561 = 'kr-C561-D';
channel488 = 'kr-C488-D';

analyze_flag = 0; % Run analysis + save
plot_flag = 1; % Show plots

t_res = 30;
xyres = 16/(60*1.5);
ntp = 33;
the_cut_t = 2;
t = (-2*t_res:t_res:(ntp-the_cut_t-1)*t_res)/60;
mask_sz = 2;
brush_sz = 3;
im_sz = 512;

%% Run annotation analysis

if analyze_flag
    im_meani488 = zeros(numel(dirs1),ntp);
    im_modei488 = zeros(numel(dirs1),ntp);
    im_mediani488 = zeros(numel(dirs1),ntp);
    im_meani568 = zeros(numel(dirs1),ntp);
    im_modei568 = zeros(numel(dirs1),ntp);
    im_mediani568 = zeros(numel(dirs1),ntp);
    hwb = waitbar(0, cat(2, 'Analyzing movie 1/', cat(2, num2str(numel(dirs1)), ' ...')));
    for ii = 1:numel(dirs1)
        if exist(dirs1{ii}, 'file') % check whether a file with this name exists
            load(dirs1{ii});
            indsep = strfind(dirs1{ii}, filesep);
            analysis_file = strrep(dirs1{ii}(indsep(end)+1:end), 'annotations', 'analysis');
            dirs1{ii} = dirs1{ii}(1:indsep(end));
        end
        cd(dirs1{ii});
        indsep = strfind(dirs1{ii}, filesep);
        filebase = dirs1{ii}(indsep(end-1)+1:indsep(end)-1);
        exist488 = 0;
        if exist(cat(2, dirs1{ii}, filebase,'_', label488, '_', channel488, '.tif'), 'file')
            im488 = tiffread(cat(2, dirs1{ii}, filebase,'_', label488, '_', channel488, '.tif'));
            exist488 = 1;
        elseif exist(cat(2,dirs1{ii},'_', channel488, '.tif'),'file')
            im488 = tiffread(cat(2,dirs1{ii},'_',channel488,'.tif'));
            exist488 = 1;
        end
        exist568 = 0;
        if exist(cat(2, dirs1{ii}, filebase,'_',label561, '_', channel561, '.tif'), 'file')
            im568 = tiffread(cat(2, dirs1{ii}, filebase,'_',label561,'_',channel561,'.tif'));
            exist568 = 1;
        elseif exist(cat(2,dirs1{ii},'_', channel561, '.tif'),'file')
            im568 = tiffread(cat(2,dirs1{ii},'_',channel561,'.tif'));
            exist568 = 1;
        end
        i488_tricell_norm = zeros(20,ntp);
        i488_tricell = zeros(20,ntp);
        i568_tricell_norm = zeros(20,ntp);
        i568_tricell = zeros(20,ntp);
        for jj = 1:ntp
            % create small mask from fiducials
            fids = ud.rfiducials(:,:,jj);
            fids = fids(fids(:,1)~=-1,:);
            [r,c] = size(fids);
            if exist488
                slice488 = squeeze(im488(:, :, jj-1));
                im_meani488(ii,jj) = mean(slice488);
                im_modei488(ii,jj) = mode(double(slice488(:)));
                im_mediani488(ii,jj) = median(double(slice488(:)));
            end
            if exist568
                slice568 = squeeze(im568(:,:,jj-1));
                im_meani568(ii,jj) = mean(slice568);
                im_modei568(ii,jj) = mode(double(slice568(:)));
                im_mediani568(ii,jj) = median(double(slice568(:)));
            end
            for kk = 1:r
                fidkk = fids(kk,1:2);
                xi = [fidkk(1)-mask_sz, fidkk(1)+mask_sz, fidkk(1)+mask_sz, fidkk(1)-mask_sz, fidkk(1)-mask_sz];
                yi = [fidkk(2)+mask_sz, fidkk(2)+mask_sz, fidkk(2)-mask_sz, fidkk(2)-mask_sz, fidkk(2)+mask_sz];
                fidmask = trajectories2mask({[xi',yi']}, ud.imsize(1:2), brush_sz, 2);
                if exist488
                    i488_tricell(kk,jj) = mean(slice488(fidmask));
                    i488_tricell_norm(kk,jj) = (mean(slice488(fidmask))-im_modei488(ii,jj))./mean(slice488); % subtract bkg & normalize to mean
                end
                if exist568
                    i568_tricell(kk,jj) = mean(slice568(fidmask));
                    i568_tricell_norm(kk,jj) = (mean(slice568(fidmask))-im_modei568(ii,jj))./mean(slice568); % subtract bkg & normalize to mean
                end
            end
        end
        im_i568_tricell_norm = mean(i568_tricell_norm);
        save(analysis_file, 'i488_tricell','i488_tricell_norm','i568_tricell','i568_tricell_norm','im_i568_tricell_norm')
        if ii < numel(dirs1)
            waitbar((ii)/numel(dirs1), hwb, cat(2, 'Analyzing movie ', cat(2, cat(2, num2str(ii+1), '/'), cat(2, num2str(numel(dirs1)), ' ...'))));
        else
            waitbar(1, hwb, 'Done!');
        end
    end
    
    if ~isempty(dirs2)
        im_meani488 = zeros(numel(dirs2),ntp);
        im_modei488 = zeros(numel(dirs2),ntp);
        im_mediani488 = zeros(numel(dirs2),ntp);
        im_meani568 = zeros(numel(dirs2),ntp);
        im_modei568 = zeros(numel(dirs2),ntp);
        im_mediani568 = zeros(numel(dirs2),ntp);
        hwb = waitbar(0, cat(2, 'Analyzing movie 1/', cat(2, num2str(numel(dirs2)), ' ...')));
        for ii = 1:numel(dirs2)
            if exist(dirs2{ii}, 'file') % check whether a file with this name exists
                load(dirs2{ii});
                indsep = strfind(dirs2{ii}, filesep);
                analysis_file = strrep(dirs2{ii}(indsep(end)+1:end), 'annotations', 'analysis');
                dirs2{ii} = dirs2{ii}(1:indsep(end));
            end
            cd(dirs2{ii});
            indsep = strfind(dirs2{ii}, filesep);
            filebase = dirs2{ii}(indsep(end-1)+1:indsep(end)-1);
            exist488 = 0;
            if exist(cat(2, dirs2{ii}, filebase,'_', label488, '_', channel488, '.tif'), 'file')
                im488 = tiffread(cat(2, dirs2{ii}, filebase,'_', label488, '_', channel488, '.tif'));
                exist488 = 1;
            elseif exist(cat(2,dirs2{ii},'_', channel488, '.tif'),'file')
                im488 = tiffread(cat(2,dirs2{ii},'_',channel488,'.tif'));
                exist488 = 1;
            end
            exist568 = 0;
            if exist(cat(2, dirs2{ii}, filebase,'_',label561, '_', channel561, '.tif'), 'file')
                im568 = tiffread(cat(2, dirs2{ii}, filebase,'_',label561,'_',channel561,'.tif'));
                exist568 = 1;
            elseif exist(cat(2,dirs2{ii},'_', channel561, '.tif'),'file')
                im568 = tiffread(cat(2,dirs2{ii},'_',channel561,'.tif'));
                exist568 = 1;
            end
            i488_tricell_norm = zeros(20,ntp);
            i488_tricell = zeros(20,ntp);
            i568_tricell_norm = zeros(20,ntp);
            i568_tricell = zeros(20,ntp);
            for jj = 1:ntp
                % create small mask from fiducials
                fids = ud.rfiducials(:,:,jj);
                fids = fids(fids(:,1)~=-1,:);
                [r,c] = size(fids);
                if exist488
                    slice488 = squeeze(im488(:, :, jj-1));
                    im_meani488(ii,jj) = mean(slice488);
                    im_modei488(ii,jj) = mode(double(slice488(:)));
                    im_mediani488(ii,jj) = median(double(slice488(:)));
                end
                if exist568
                    slice568 = squeeze(im568(:,:,jj-1));
                    im_meani568(ii,jj) = mean(slice568);
                    im_modei568(ii,jj) = mode(double(slice568(:)));
                    im_mediani568(ii,jj) = median(double(slice568(:)));
                end
                for kk = 1:r
                    fidkk = fids(kk,1:2);
                    xi = [fidkk(1)-mask_sz, fidkk(1)+mask_sz, fidkk(1)+mask_sz, fidkk(1)-mask_sz, fidkk(1)-mask_sz];
                    yi = [fidkk(2)+mask_sz, fidkk(2)+mask_sz, fidkk(2)-mask_sz, fidkk(2)-mask_sz, fidkk(2)+mask_sz];
                    xi(xi == 0) = 3;
                    yi(yi == 0) = 3;
                    fidmask = trajectories2mask({[xi',yi']}, ud.imsize(1:2), brush_sz, 2);
                    if exist488
                        i488_tricell(kk,jj) = mean(slice488(fidmask));
                        i488_tricell_norm(kk,jj) = (mean(slice488(fidmask))-im_modei488(ii,jj))./mean(slice488); % subtract bkg & normalize to mean
                    end
                    if exist568
                        i568_tricell(kk,jj) = mean(slice568(fidmask));
                        i568_tricell_norm(kk,jj) = (mean(slice568(fidmask))-im_modei568(ii,jj))./mean(slice568); % subtract bkg & normalize to mean
                    end
                end
            end
            im_i568_tricell_norm = mean(i568_tricell_norm);
            save(analysis_file, 'i488_tricell','i488_tricell_norm','i568_tricell','i568_tricell_norm','im_i568_tricell_norm')
            if ii < numel(dirs2)
                waitbar((ii)/numel(dirs2), hwb, cat(2, 'Analyzing movie ', cat(2, cat(2, num2str(ii+1), '/'), cat(2, num2str(numel(dirs2)), ' ...'))));
            else
                waitbar(1, hwb, 'Done!');
            end
        end
    end
end

%% Read in data for plots

if plot_flag
    dirs1 = dirs1_orig;
    dirs2 = dirs2_orig;
    
    all_i488_tricell_1 = [];
    all_i488_tricell_norm_1 =[];
    all_i568_tricell_1 = [];
    all_i568_tricell_norm_1 =[];
    all_im_i568_tricell_norm_1 = [];
    
    for ii = 1:numel(dirs1)
        if exist(dirs1{ii}, 'file')
            load(strrep(dirs1{ii}, 'annotations', 'analysis'));
        else
            disp(cat(2, 'Error when loading ', dirs1{ii}));
        end
        
        del = find(i488_tricell_norm(:,1)==0);
        i488_tricell_norm(del,:) = [];
        [r,c] = size(i488_tricell_norm);
        all_i488_tricell_1 = [all_i488_tricell_1; i488_tricell];
        all_i488_tricell_norm_1 = [all_i488_tricell_norm_1; i488_tricell_norm];
        
        del = find(i568_tricell_norm(:,1)==0);
        i568_tricell_norm(del,:) = [];
        [r,c] = size(i568_tricell_norm);
        all_i568_tricell_1 = [all_i568_tricell_1; i568_tricell];
        all_i568_tricell_norm_1 = [all_i568_tricell_norm_1; i568_tricell_norm];
        all_im_i568_tricell_norm_1 = [all_im_i568_tricell_norm_1; im_i568_tricell_norm];
    end
    
    all_i488_tricell_norm_1(all_i488_tricell_norm_1(:,1)==0,:) = [];
    all_i568_tricell_norm_1(all_i568_tricell_norm_1(:,1)==0,:) = [];
    
    cm1 = jet(size(all_i488_tricell_1,1));
    all_i488_tricell_norm_norm_1 = zeros(size(all_i488_tricell_norm_1));
    all_i568_tricell_norm_norm_1 = zeros(size(all_i568_tricell_norm_1));
    
    all_i488_tricell_2 = [];
    all_i488_tricell_norm_2 =[];
    all_i568_tricell_2 = [];
    all_i568_tricell_norm_2 =[];
    all_im_i568_tricell_norm_2 = [];
    
    if ~isempty(dirs2)
        for ii = 1:numel(dirs2)
            if exist(dirs2{ii}, 'file')
                load(strrep(dirs2{ii}, 'annotations', 'analysis'));
            else
                disp(cat(2, 'Error when loading ', dirs1{ii}));
            end
            
            del = find(i488_tricell_norm(:,1)==0);
            i488_tricell_norm(del,:) = [];
            [r,c] = size(i488_tricell_norm);
            all_i488_tricell_2 = [all_i488_tricell_2; i488_tricell];
            all_i488_tricell_norm_2 = [all_i488_tricell_norm_2; i488_tricell_norm];
            
            del = find(i568_tricell_norm(:,1)==0);
            i568_tricell_norm(del,:) = [];
            [r,c] = size(i568_tricell_norm);
            all_i568_tricell_2 = [all_i568_tricell_2; i568_tricell];
            all_i568_tricell_norm_2 = [all_i568_tricell_norm_2; i568_tricell_norm];
            all_im_i568_tricell_norm_2 = [all_im_i568_tricell_norm_2; im_i568_tricell_norm];
        end
        
        all_i488_tricell_norm_2(all_i488_tricell_norm_2(:,1)==0,:) = [];
        all_i568_tricell_norm_2(all_i568_tricell_norm_2(:,1)==0,:) = [];
        
        cm2 = jet(size(all_i488_tricell_2,1));
        all_i488_tricell_norm_norm_2 = zeros(size(all_i488_tricell_norm_2));
        all_i568_tricell_norm_norm_2 = zeros(size(all_i568_tricell_norm_2));
    end
    
    %% Make plots
    
    figure; hold on;
    for ii = 1:size(all_i568_tricell_norm_1,1)
        plot(t,all_i568_tricell_norm_1(ii,:),'Color',cm1(ii,:),'LineWidth', 4);
        all_i568_tricell_norm_norm_1(ii,:) = all_i568_tricell_norm_1(ii,:)./all_i568_tricell_norm_1(ii,1);
    end
    
    all_i568_tricell_norm_norm_1(all_i568_tricell_norm_norm_1(:,end)>5 | all_i568_tricell_norm_norm_1(:,end)<-3| all_i568_tricell_norm_norm_1(:,7)>3,:) = [];
    
    figure; hold on;
    for ii = 1:size(all_i568_tricell_norm_norm_1,1)
        plot(t,all_i568_tricell_norm_norm_1(ii,:),'Color',cm1(ii,:),'LineWidth',4);
    end
    
    figure; hold on;
    for ii = 1:size(all_i568_tricell_norm_2,1)
        plot(t,all_i568_tricell_norm_2(ii,:),'Color',cm2(ii,:),'LineWidth', 4);
        all_i568_tricell_norm_norm_2(ii,:) = all_i568_tricell_norm_2(ii,:)./all_i568_tricell_norm_2(ii,1);
    end
    
    all_i568_tricell_norm_norm_2(all_i568_tricell_norm_norm_2(:,end)>3 | all_i568_tricell_norm_norm_2(:,end)<-3 | all_i568_tricell_norm_norm_2(:,7)>3,:) = [];
    
    figure; hold on;
    for ii = 1:size(all_i568_tricell_norm_norm_2,1)
        plot(t,all_i568_tricell_norm_norm_2(ii,:),'Color',cm2(ii,:),'LineWidth',4);
    end
    
    [m_all_i568_tricell_norm_norm_1,~,e_all_i568_tricell_norm_norm_1] = rmsnan(all_i568_tricell_norm_1);
    [m_all_i568_tricell_norm_norm_2,~,e_all_i568_tricell_norm_norm_2] = rmsnan(all_i568_tricell_norm_2);
    
    f4=figure; hold on;
    errorbar(t,m_all_i568_tricell_norm_norm_1,e_all_i568_tricell_norm_norm_1,'b','LineWidth',2)
    plot(t,m_all_i568_tricell_norm_norm_1,'Color', 'b', 'LineWidth', 4)
    errorbar(t,m_all_i568_tricell_norm_norm_2,e_all_i568_tricell_norm_norm_2,'r','LineWidth',2)
    plot(t,m_all_i568_tricell_norm_norm_2,'Color','r','LineWidth',4)
    set(gca, 'XTick', [0:5:25], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlim([-2 16])
    xlabel('time (min)')
    ylabel('intensity (a.u.)')
    
    ratios_1 = 100*(mean(all_i568_tricell_norm_1(:,end-1:end),2)-mean(all_i568_tricell_norm_1(:,1:2),2))./mean(all_i568_tricell_norm_1(:,1:2),2);
    figure; boxplotMWstd(0.5:1.5,ratios_1,2,'bo',1);
    hold on;
    ratios_2 = 100*(mean(all_i568_tricell_norm_2(:,end-1:end),2)-mean(all_i568_tricell_norm_2(:,1:2),2))./mean(all_i568_tricell_norm_2(:,1:2),2);
    boxplotMWstd(2:3,ratios_2,2,'ro',1);
    xlim([0 3.5])
    box off;
    set(gca, 'LineWidth', 2);
    set(gca, 'XTick', [1.0 2.5], 'XTickLabel', legend_str, 'YTickLabel', [-1000:50:1000], 'YTick', [-1000:50:1000], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('% intensity change', 'FontSize', 32, 'FontWeight', 'bold');
    plot([0 4],[0 0],'k--','LineWidth',2)
    ylim([-150 150])

    disp('group 1')
    fprintf('mean = %.2f, ste = %.2f, P = %.3f\n',mean(ratios_1),std(ratios_1)/sqrt(length(ratios_1)), signrank(ratios_1))
    disp('group 2')
    fprintf('mean = %.2f, ste = %.2f, P = %.3f\n',mean(ratios_2),std(ratios_2)/sqrt(length(ratios_2)), signrank(ratios_2))
    disp('between groups')
    fprintf('P = %.3f\n',2*normpdf(rod_mannwhitney(ratios_1,ratios_2)))
end
