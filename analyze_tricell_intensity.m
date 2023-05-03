%% Read in fiducials marking tricell jxns & calculate intensity over time

%% Define annotation files

files = {
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190524/Embryo2/annotations_tricell.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo1/annotations_tricell2.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo2/annotations_tricell.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo3/annotations_tricell2.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo4/annotations_tricell2.mat';
    };

dirs1_orig = files;
dirs2_orig = [];

%% Define necessary parameters

dirs1 = dirs1_orig;
dirs2 = dirs2_orig;
legend_str = {'Ecad','Rap1'};
label561 = 'Ecad';
label488 = 'Rap1';
channel561 = 'kr-C561-D';
channel488 = 'kr-C488-D';

analyze_flag = 1; % Run analysis + save
plot_flag = 1; % Show plots

t_res = 15;
xyres = 16/(60*1.5);
ntp = 61;
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
        if exist(cat(2,dirs1{ii},'localz_Rap1.tif'),'file')
            im488 = tiffread(cat(2,dirs1{ii},'localz_Rap1.tif'));
            exist488 = 1;
        end
        exist568 = 0;
        if exist(cat(2,dirs1{ii},'localz_ecad.tif'),'file')
            im568 = tiffread(cat(2,dirs1{ii},'localz_ecad.tif'));
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
                data = double(slice488(:));
%                 thresh = prctile(data,25);
                im_modei488(ii,jj) = mode(data);
                im_mediani488(ii,jj) = median(data);
            end
            if exist568
                slice568 = squeeze(im568(:,:,jj-1));
                im_meani568(ii,jj) = mean(slice568);
                data = double(slice568(:));
%                 thresh = prctile(data,25);
                im_modei568(ii,jj) = mode(data);
                im_mediani568(ii,jj) = median(data);
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
        im_i488_tricell_norm = mean(i488_tricell_norm);
        im_i568_tricell_norm = mean(i568_tricell_norm);
        save(analysis_file, 'i488_tricell','i488_tricell_norm','i568_tricell','i568_tricell_norm','im_i488_tricell_norm','im_i568_tricell_norm')
        if ii < numel(dirs1)
            waitbar((ii)/numel(dirs1), hwb, cat(2, 'Analyzing movie ', cat(2, cat(2, num2str(ii+1), '/'), cat(2, num2str(numel(dirs1)), ' ...'))));
        else
            waitbar(1, hwb, 'Done!');
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
    all_im_i488_tricell_norm_1 = [];
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
        all_im_i488_tricell_norm_1 = [all_im_i488_tricell_norm_1; im_i488_tricell_norm];
        
        del = find(i568_tricell_norm(:,1)==0);
        i568_tricell_norm(del,:) = [];
        [r,c] = size(i568_tricell_norm);
        %     all_i568_tricell_norm(1:r,:) = i568_tricell_norm;
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
        end
        
        all_i488_tricell_norm_2(all_i488_tricell_norm_2(:,1)==0,:) = [];
        all_i568_tricell_norm_2(all_i568_tricell_norm_2(:,1)==0,:) = [];
        
        cm2 = jet(size(all_i488_tricell_2,1));
        all_i488_tricell_norm_norm_2 = zeros(size(all_i488_tricell_norm_2));
        all_i568_tricell_norm_norm_2 = zeros(size(all_i568_tricell_norm_2));
    end
    
    
    %% Make plots
    
    %        individual 488 intensities vs time
    f1=figure; hold on;
    for ii = 1:size(all_i488_tricell_norm_1,1)
        plot(t,all_i488_tricell_norm_1(ii,:),'Color',cm1(ii,:),'LineWidth', 4);
        all_i488_tricell_norm_norm_1(ii,:) = all_i488_tricell_norm_1(ii,:)./mean(all_i488_tricell_norm_1(ii,1:2),2);
    end
    
    all_i488_tricell_norm_norm_1(all_i488_tricell_norm_norm_1(:,end)>5 | all_i488_tricell_norm_norm_1(:,end)<-5,:) = [];
    
    %        individual 568 intensities vs time
    f2=figure; hold on;
    for ii = 1:size(all_i568_tricell_norm_1,1)
        plot(t,all_i568_tricell_norm_1(ii,:),'Color',cm1(ii,:),'LineWidth', 4);
        all_i568_tricell_norm_norm_1(ii,:) = all_i568_tricell_norm_1(ii,:)./mean(all_i568_tricell_norm_1(ii,1:2),2);
    end
    
    all_i568_tricell_norm_norm_1(all_i568_tricell_norm_norm_1(:,end)>5 | all_i568_tricell_norm_norm_1(:,end)<-5,:) = [];
    
    %        individual norm 568 intensities vs time
    f2=figure; hold on;
    for ii = 1:size(all_i568_tricell_norm_norm_1,1)
        plot(t,all_i568_tricell_norm_norm_1(ii,:),'Color',cm1(ii,:),'LineWidth', 4);
    end

    %        means 488 vs 568
    f3=figure; hold on;
    [m_all_i488_tricell_norm,~,e_all_i488_tricell_norm] = rmsnan(all_i488_tricell_norm_1);
    [m_all_i568_tricell_norm,~,e_all_i568_tricell_norm] = rmsnan(all_i568_tricell_norm_1);
    errorbar(t,m_all_i488_tricell_norm,e_all_i488_tricell_norm,'g')
    errorbar(t,m_all_i568_tricell_norm,e_all_i568_tricell_norm,'r')
    
    %        normalized 488 vs 568
    figure; hold on;
    [m_all_i488_tricell_norm_norm,~,e_all_i488_tricell_norm_norm] = rmsnan(all_i488_tricell_norm_1);
    [m_all_i568_tricell_norm_norm,~,e_all_i568_tricell_norm_norm] = rmsnan(all_i568_tricell_norm_1);
    h1 = errorbar(t,100*m_all_i488_tricell_norm_norm,100*e_all_i488_tricell_norm_norm,'LineWidth',2);
    plot(t,100*m_all_i488_tricell_norm_norm,'k','LineWidth',4)
    set(h1,'Color',[0.6 0.6 0.6])
    set(gca, 'XTick', [0:5:25], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlim([-1 15])
    xlabel('time (min)')
    ylabel('intensity (a.u.)')

    figure; hold on;
    h2 = errorbar(t,100*m_all_i568_tricell_norm_norm,100*e_all_i568_tricell_norm_norm,'LineWidth',2);
    set(h2,'Color',[0.6 0.6 0.6])
    plot(t,100*m_all_i568_tricell_norm_norm,'k','LineWidth',4)
    set(gca, 'XTick', [0:5:25], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlim([-1 15])
    xlabel('time (min)')
    ylabel('intensity (a.u.)')
     
    ratios_488 = 100*(mean(all_i488_tricell_norm_1(:,end-1:end),2)-mean(all_i488_tricell_norm_1(:,1:2),2))./mean(all_i488_tricell_norm_1(:,1:2),2);
    ratios_1 = 100*(mean(all_i568_tricell_norm_1(:,end-1:end),2)-mean(all_i568_tricell_norm_1(:,1:2),2))./mean(all_i568_tricell_norm_1(:,1:2),2);

    
    ratios_im_i568 = mean(all_im_i568_tricell_norm_1(:,end-1:end),2)./mean(all_im_i568_tricell_norm_1(:,1:2),2);
    ratios_im_i488 = mean(all_im_i488_tricell_norm_1(:,end-1:end),2)./mean(all_im_i488_tricell_norm_1(:,1:2),2);
    figure; boxplotMWstd(0.5:1.5,ratios_1,2,'ko',1);
    xlim([0 2])
    box off;
    set(gca, 'LineWidth', 2);
    set(gca, 'XTick', [1.0], 'XTickLabel', 'E-cadherin', 'YTick', [-1000:100:1000], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('% intensity change', 'FontSize', 32, 'FontWeight', 'bold');
    plot([0 4],[0 0],'k--','LineWidth',2)
    ylim([-150 400])
    set(gcf,'Position',[235 254 383 420])
    disp('568 percent sig')
    fprintf('mean = %.2f, ste = %.2f, P = %.3f\n',mean(ratios_1),std(ratios_1)/sqrt(length(ratios_1)), signrank(ratios_1))
    
    figure; boxplotMWstd(0.5:1.5,ratios_488,2,'ko',1);
    xlim([0 2])
    box off;
    set(gca, 'LineWidth', 2);
    set(gca, 'XTick', [1.0], 'XTickLabel', 'Rap1', 'YTick', [-1000:100:1000], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('% intensity change', 'FontSize', 32, 'FontWeight', 'bold');
    plot([0 4],[0 0],'k--','LineWidth',2)
    ylim([-150 400])
    set(gcf,'Position',[235 254 383 420])
    disp('488 percent sig')
    fprintf('mean = %.2f, ste = %.2f, P = %.3f\n',mean(ratios_488),std(ratios_488)/sqrt(length(ratios_488)), signrank(ratios_488))

end
