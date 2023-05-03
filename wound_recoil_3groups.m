%% To extract retraction velocities before/after cuts.
% Input:
%   - fiducials before and after ablation
% Output:
%   - retraction velocity

clear all;
close all;

%% declare directories


dirs_con = {
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210421/Embryo3/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210421/Embryo4/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210421/Embryo5/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210421/Embryo7/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210426/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210426/Embryo4/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210426/Embryo5/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210503/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210503/Embryo3/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210503/Embryo4/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210517b/Embryo1/Cut/',...
    };

dirs_DN = {
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210423/Embryo3/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210423/Embryo6/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210423/Embryo7/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210428/Embryo10/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210506b/Embryo1/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210506b/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210508/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210508/Embryo5/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210508/Embryo6/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210508/Embryo7/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210522/Embryo7/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210522/Embryo8/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210519/Embryo2/Cut/',...
    };

dirs_CA = {
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210422/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210428/Embryo3/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210428/Embryo4/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210501/Embryo1/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210501/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210505/Embryo1/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210505/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210505/Embryo4/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210510c/Embryo2/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210522/Embryo1/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210526/Embryo6/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210526/Embryo12/Cut/',...
    '/Users/katy/Documents/Data/2021 daGal4 Rap1CA EcadtdTo sqhGFP Wound Cut/20210519/Embryo3/Cut/',...
    };

dirs1 = dirs_con;
dirs2 = dirs_DN;
dirs3 = dirs_CA;

legend_str = {'control','Rap1DN','Rap1CA'};

%% declare variables
 
directories = {dirs1; dirs2; dirs3};

t_res = 4;          % 4 seconds elapse between frames
time_factor = t_res / 5.; % This script was written with a 5sec resolution in mind.
magn = 60;
um_per_pix = 16/(1.5*magn);
Da = 0.0;

ud = [];
dr_um_control = [];
dth_um_control = [];
cp_um_control = [];
cp_ang_control = [];
initial_dir = pwd;

tpo = 2;
tpf = 3;
cut_t = 2;
the_cut_t = cut_t;



marker_str = 'so^';
color_str = 'brg';
line_str = {'-'; '-'; '-'};
title_label = '';
save_flag = 0;
noplots = 0;

% calculate retraction time [s]
exp_time = 0.150;
num_slices = 11;
num_spots = 1;
num_pulses = 10; % per spot
pulse_duration = 0.067; % per pulse
num_lasers = 2;
t_retraction = ((exp_time*num_slices)*num_lasers)+(num_spots*num_pulses*pulse_duration); 

%% analysis - calculate length between fiducials in each annotated frame
h = waitbar(0, cat(2, 'Analyzing movie 1/', cat(2, num2str(numel(directories)), ' ...')));

for jj = 1:numel(directories)
    for ii = 1:numel(directories{jj})
        d = directories{jj};
        
        % Load annotations file.
        fname = cat(2, d{ii}, 'fiducials.mat');
        fname_myosin = cat(2,d{ii},'annotations_myosin.mat');
        if exist(fname, 'file')
            vars = load(fname);
            ud = vars.ud;
            cd(d{ii});

            ind = find(ud.rnfiducials > 0);
            if ~isempty(ind)
                p1 = ud.rfiducials(1, 1:2, tpo:tpf);
                tmp = nonzeros(p1); % took out "+!"
                p1 = reshape(tmp, [numel(tmp)/(size(p1, 2) * size(p1, 3)) size(p1, 2) size(p1, 3)]); % took out "-!"
                p1 = reshape(p1, [size(p1, 3), 2, 1]);
                p1 = reshape(p1, [2, size(p1, 1), 1])';

                p2 = ud.rfiducials(2, 1:2, tpo:tpf);
                tmp = nonzeros(p2); % took out "+!"
                p2 = reshape(tmp, [numel(tmp)/(size(p2, 2) * size(p2, 3)) size(p2, 2) size(p2, 3)]); % took out "-!"
                p2 = reshape(p2, [size(p2, 3), 2, 1]);
                p2 = reshape(p2, [2, size(p2, 1), 1])';

                if isfield(ud, 'rangles')
                    angles = ud.rangles;
                    [l a p angle_values] = cut_analysis(p1, p2, ud.rpolygons(:, :, tpo:tpf), angles(:, :, tpo:tpf), cut_t, t_res, 0, 0);
                else
                    angles = [];
                    [l a p angle_values] = cut_analysis(p1, p2, [], [], cut_t, t_res, 0, 0);
                end

                recoil = (l(2) - l(1)) * um_per_pix ./ t_retraction;
                l = l * um_per_pix;

                save('vars.mat', 'l', 'recoil');
            end

        else
            disp(cat(2, fname, ' does not exist.'));
        end
        
        if exist(fname_myosin,'file')
            temp = load(fname_myosin);
            im = tiffread(cat(2,'_kr-C488-D.tif'));
            slice = squeeze(im(:, :, 1));
            ud = temp.ud;
            cd(d{ii});
            trace = ud.rpolygons(1);
            msk_trace = trajectories2mask(trace, [512 512], 3, 0);
            themode = mode(reshape(double(slice), [1 prod(size(slice))]));
            avg_myosin = mean(slice(msk_trace)-themode);
            
            save('vars.mat', 'avg_myosin','-append');
        else
            disp(cat(2,fname_myosin,' does not exist.'));
        end

        if ii < numel(d)
            waitbar((ii)/numel(d), h, cat(2, 'Analyzing movie ', cat(2, cat(2, num2str(ii+1), '/'), cat(2, num2str(numel(directories)), ' ...'))));
        else
            waitbar(1, h, 'Done!');
        end
    end
end

cd(initial_dir);

close(h);


%% compile data
clear vars1;
clear l1;
clear avg_myosin;

vars1 = zeros(1,numel(dirs1));
l1 = zeros(1,numel(dirs1));
myo1 = zeros(1,numel(dirs1));
for ii = 1:numel(dirs1)
    fname = cat(2, dirs1{ii}, 'vars.mat');
    
    if exist(fname, 'file')
        load(fname);
        
        vars1(ii) = recoil;
        l1(ii) = l(1);
        myo1(ii) = avg_myosin;
    end
end
[m(1) s(1) e(1)] = rmsnan(vars1');

clear vars2;
clear l2;
clear avg_myosin;
l2 = zeros(1,numel(dirs2));
vars2 = zeros(1,numel(dirs2));
myo2 = zeros(1,numel(dirs2));
for ii = 1:numel(dirs2)
    fname = cat(2, dirs2{ii}, 'vars.mat');
    
    if exist(fname, 'file')
        load(fname);
        
        vars2(ii) = recoil;
        l2(ii) = l(1);
        myo2(ii) = avg_myosin;
    end
end
[m(2) s(2) e(2)] = rmsnan(vars2');

clear vars3;
clear l3;
clear avg_myosin;
l3 = [];
vars3 = [];
myo3 = zeros(1,numel(dirs3));
for ii = 1:numel(dirs3)
    fname = cat(2, dirs3{ii}, 'vars.mat');
    
    if exist(fname, 'file')
        load(fname);
        
        vars3(ii) = recoil;
        l3(ii) = l(1);
        myo3(ii) = avg_myosin;
    end
end
[m(3) s(3) e(3)] = rmsnan(vars3');

%% make plots - boxplot, recoil velocity
figure;
boxplotMWstd(0.5:1.5,vars1,3,'bo',1);
hold on;
boxplotMWstd(2.0:3.0,vars2,3,'ro',1);
hold on;
boxplotMWstd(3.5:4.5,vars3,3,'go',1);
hold on;

ylim([0 2]);
xlim([0 5]);
box off;
set(gca, 'LineWidth', 2);
set(gca, 'XTick', [1.0 2.5 4.0], 'XTickLabel', legend_str, 'YTickLabel', [0:0.5:10], 'YTick', [0:0.5:10], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('retraction velocity (µm/s)', 'FontSize', 32, 'FontWeight', 'bold');

rod_dunn(vars1, vars2, vars3)


%% make plots - velocity vs initial length
% first dataset - scatter plot
figure;
p = polyfit(l1, vars1,1);
f = polyval(p,l1);
hold on;
plot(l1,f,'-','Color',[0.7 0.7 0.7], 'LineWidth',3);

scatter(l1, vars1, 100, 'o', 'b', 'filled','LineWidth',1);
set(gca, 'LineWidth', 2);
ylim([0 1]);
xlim([3 17]);
box off;
set(gca, 'XTick', [1:4:21], 'XTickLabel', [1:4:21], 'YTickLabel', [0:0.5:2], 'YTick', [0:0.5:2], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
xlabel('initial length (µm)', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('retraction velocity (µm/s)', 'FontSize', 32, 'FontWeight', 'bold');

[r, pr] = corrcoef(l1, vars1);
coefficient(1) = r(1,2);

% second dataset - scatter plot
% figure;
p = polyfit(l2, vars2,1);
f = polyval(p,l2);
hold on;
plot(l2,f,'-','Color',[0.7 0.7 0.7], 'LineWidth',3);

scatter(l2, vars2, 100, 'o', 'r', 'filled','LineWidth',1);
set(gca, 'LineWidth', 2);
ylim([0 1]);
xlim([2 20]);
box off;
set(gca, 'XTick', [0:4:21], 'XTickLabel', [1:4:21], 'YTickLabel', [0:0.5:2], 'YTick', [0:0.5:2], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
xlabel('initial length (µm)', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('retraction velocity (µm/s)', 'FontSize', 32, 'FontWeight', 'bold');

[r, pr] = corrcoef(l2, vars2);
coefficient(2) = r(1,2);

% third dataset - scatter plot
figure;
p = polyfit(l3, vars3,1);
f = polyval(p,l3);
hold on;
plot(l3,f,'-','Color',[0.7 0.7 0.7], 'LineWidth',3);

scatter(l3, vars3, 100, 'o', 'r', 'filled','LineWidth',1);
set(gca, 'LineWidth', 2);
ylim([0 1]);
xlim([2 20]);
box off;
set(gca, 'XTick', [0:4:21], 'XTickLabel', [1:4:21], 'YTickLabel', [0:0.5:2], 'YTick', [0:0.5:2], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
xlabel('initial length (µm)', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('retraction velocity (µm/s)', 'FontSize', 32, 'FontWeight', 'bold');

[r, pr] = corrcoef(l3, vars3);
coefficient(3) = r(1,2);

%% make plots - boxplot, recoil velocity
figure;
boxplotMWstd(0.5:1.5,myo1,3,'bo',1);
hold on;
boxplotMWstd(2.0:3.0,myo2,3,'ro',1);
hold on;
boxplotMWstd(3.5:4.5,myo3,3,'go',1);

% ylim([0 2]);
xlim([0 5]);
box off;
set(gca, 'LineWidth', 2);
xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('myosin intensity (a.u.)', 'FontSize', 32, 'FontWeight', 'bold');

rod_dunn(myo1, myo2,myo3)

%% make plots - velocity vs myosin
% first dataset - scatter plot
figure;
p = polyfit(myo1, vars1,1);
f = polyval(p,myo1);
hold on;
plot(myo1,f,'b-', 'LineWidth',3);
scatter(myo1, vars1, 100, 'o', 'b', 'filled','LineWidth',1);

p = polyfit(myo2, vars2,1);
f = polyval(p,myo2);
hold on;
plot(myo2,f,'r-', 'LineWidth',3);
scatter(myo2, vars2, 100, 'o', 'r', 'filled','LineWidth',1);

p = polyfit(myo3,vars3,1);
f = polyval(p,myo3);
plot(myo3,f,'g-','LineWidth',3);
scatter(myo3, vars3, 100,'o','g','filled','LineWidth',1);

set(gca, 'LineWidth', 2);
box off;
xlabel('myosin intensity (a.u.)', 'FontSize', 32, 'FontWeight', 'bold');
ylabel('retraction velocity (µm/s)', 'FontSize', 32, 'FontWeight', 'bold');
