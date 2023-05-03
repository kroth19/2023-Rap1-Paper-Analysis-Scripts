%% Read in output from SIESTA "calculate polyline intensity between fiducials" and plot

%% define data

con = {
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190524/Embryo2/analysis_linearjxn_Ecad.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo1/analysis_linearjxn_Ecad.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo2/analysis_linearjxn_Ecad.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo3/analysis_linearjxn_Ecad.mat';
    '/Users/katy/Documents/Data/2019 ECad-tdTo Rap1-GFP Wound Short/20190606/Embryo4/analysis_linearjxn_Ecad.mat';
    };

    
%% load data

dirs1 = con;
dirs2 = [];
legend_str = {'Ecad'};

t_res = 15;
xyres = 16/(60*1.5);
ntp = 63;
the_cut_t = 2;
t = (-2*t_res:t_res:(ntp-the_cut_t-1)*t_res)/60;

all_mean_intensities_1 = [];
im_mean_intensities_1 = [];
for i = 1:length(dirs1)
    load(dirs1{i});
    C = mean(int_after_bleach_corr(:,200:800,:),2);
    C = reshape(C,size(C,1),size(C,3));
    all_mean_intensities_1 = [all_mean_intensities_1 C(1:ntp,:)];
    im_mean_intensities_1 = [im_mean_intensities_1; mean(C(1:ntp,:)')];
    clear int_after_bleach_corr
end

all_mean_intensities_2 =[];
im_mean_intensities_2 = [];
for i = 1:length(dirs2)
    load(dirs2{i});
    C = mean(int_after_bleach_corr(:,200:800,:),2);
    C = reshape(C,size(C,1),size(C,3));
    all_mean_intensities_2 = [all_mean_intensities_2 C(1:ntp,:)];
    im_mean_intensities_2 = [im_mean_intensities_2; mean(C(1:ntp,:)')];
    clear int_after_bleach_corr
end

%% plot data

cm1 = jet(length(dirs1));

im_mean_intensities_1 = im_mean_intensities_1';

% plot all normalized data - group 1
figure; hold on;
im_mean_intensities_norm_1 = zeros(size(im_mean_intensities_1));
for i = 1:size(im_mean_intensities_1,2)
    im_mean_intensities_norm_1(:,i) = im_mean_intensities_1(:,i)./im_mean_intensities_1(1,i);
    plot(t,im_mean_intensities_1(:,i),'Color',cm1(i,:),'LineWidth',2)
end

ratios_1 = 100*(mean(all_mean_intensities_1(end-1:end,:))-mean(all_mean_intensities_1(1:2,:)))./mean(all_mean_intensities_1(1:2,:));
figure; boxplotMWstd(0.5:1.5,ratios_1,2,'ko',1);
hold on;
xlim([0 2])
box off;
set(gca, 'LineWidth', 2);
set(gca, 'XTick', [1.0], 'XTickLabel', legend_str, 'YTick', [-200:100:1000], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('% intensity change', 'FontSize', 32, 'FontWeight', 'bold');
plot([0 4],[0 0],'k--','LineWidth',2)
set(gcf,'Position',[235 254 383 420])
ylim([-150 400])

disp('linear percent sig')
fprintf('mean = %.2f, ste = %.2f, P = %.3f\n',mean(ratios_1),std(ratios_1)/sqrt(length(ratios_1)), signrank(ratios_1))

linear_488_end = mean(all_mean_intensities_1(end-1:end,:));
linear_488_start = mean(all_mean_intensities_1(1:2,:));

figure; hold on;
boxplotMWstd_mod(0.5:1.5,linear_488_start,2,'ko',1);
boxplotMWstd_mod(2:3,linear_488_end,2,'ko',1);