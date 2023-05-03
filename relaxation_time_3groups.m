%% calculate relaxation time by modeling Kelvin-Voigt

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


% Variables that need adjustment
dirs1 = dirs_con; % First group of folders to analyze.
dirs2 = dirs_DN; % Second.
dirs3 = dirs_CA; % Third.

noplots = 0;

legend_str = {'control'; 'Rap1DN'; 'Rap1CA'}; % Legend of the figure.
color_str = 'brg'; % Line colors for the different groups (initial letters).
title_label = ''; % A title for all the figures.
save_flag = 0; % Set to 1 if you want to save all the plots.

t_resolution = 4; % Time interval in seconds.
t_first = [(670+11*150*2)/1000. (670+11*150*2)/1000. (670+11*150*2)/1000.]; % Time in between the last image before ablation and the firt image after ablation. It is typically .67 seconds for the spot ablation and whatever the exposure time * the number of slices collected. This time is important to calculate the initial velocity of retraction. For your original experiments with dirs_LAP, dirs_IAP, dirs_DV, this value should be 4 or 5 sec, because you were taking a single movie and shooting the laser shortly after acquisition of the last movie before ablation.

ntp = 21; % Total number of time points.
the_cut_t = 2; % Index of the last time point before ablation.
um_per_pix = 16/90;

% End of variables that need adjustment.




marker_str = 'so^';
line_str = {'-'; '-'; '-'};
bins = (-25:25:200);
R_min = 0.96; % Minimum goodness of fit to consider the data from an embryo as properly fit.
T_max = 30; % Maximum relaxation time. Longer times indicate that the curve does not relax before the end of the measured time.


directories = dirs1;

rl_i = zeros(numel(directories), ntp).*nan;
rl_i_abs = zeros(numel(directories), ntp).*nan;
ra_i = zeros(numel(directories), ntp).*nan;
rp_i = zeros(numel(directories), ntp).*nan;
rsf_i = zeros(numel(directories), ntp).*nan;

binl_i = zeros(numel(directories), numel(bins) - 1) .* nan;
binv_i = zeros(numel(directories), numel(bins) - 1) .* nan;
bina_i = zeros(numel(directories), numel(bins) - 1) .* nan;
binp_i = zeros(numel(directories), numel(bins) - 1) .* nan;

maxv_i = zeros(numel(directories), 2); % [vmax d0] matrix
maxvind_i = zeros(numel(directories), 2); %[vmax ind] matrix

tau_i = zeros(numel(directories), 2); % Model parameters.
taudamp_i = zeros(numel(directories), 3); % Model parameters.
R_i = [];
Rdamp_i = [];
angv_i = zeros(numel(directories), 2); % Average angle and max velocity.
angla_i = zeros(numel(directories), 2); % Average angle and max abs retraction.
angln_i = zeros(numel(directories), 2); % Average angle and max normalized retraction.
singleangv_i = zeros(0, 2); % Angle and max velocity.

nint_i = zeros(numel(directories, 1)); % Number of interfaces in the cable.
cl_i = zeros(numel(directories, 1)); % Cable length.

for i = 1:numel(directories)
    if exist(cat(2, directories{i}, 'analysis_local.mat'), 'file')
        load(cat(2, directories{i}, 'analysis_local'));
    elseif exist(cat(2, directories{i}, 'analysis_pairdist_onemin.mat'), 'file')
        load(cat(2, directories{i}, 'analysis_pairdist_onemin'));
    else
        disp(cat(2, 'Error when loading analysis results for ', directories{i}));
    end
    
    tau_i(i, :) = tau'
    if exist('taudamp', 'var')
        taudamp_i(i, :) = taudamp';
    end
    %startr = max(cut_t - (the_cut_t - 1), 1);
    startr = cut_t - (the_cut_t - 1);
    if startr > 0
        startw = 1;
    else
        startw = abs(startr-2);
        startr = 1;
    end
    
    finishr = cut_t + ntp - the_cut_t;
    finishw = startw + finishr - startr;
    
    if finishr > numel(l)
        finishr = cut_t + numel(l) - the_cut_t;
        finishw = startw + finishr - startr;
    end
    
    rl_i(i, startw:finishw) = l(startr:finishr);
    rl_i_abs(i, :) = (rl_i(i, :) - l(cut_t)) .* um_per_pix;
    R_i(end+1) = R;
    if exist('taudamp', 'var')
        Rdamp_i(end+1) = Rdamp;
    end
    rl_i(i, :) = 100.*(rl_i(i, :) - l(cut_t))./l(cut_t);
    
    
    if ~ isempty(a)
        ra_i(i, startw:finishw) = a(startr:finishr);
        %ra_i(i, :) = (ra_i(i, :) - a(cut_t))./(4.9^2.);
        ra_i(i, :) = 100.*(ra_i(i, :) - a(cut_t))./a(cut_t);
    end
    
    if ~ isempty(p)
        rp_i(i, startw:finishw) = p(startr:finishr);
        %rp_i(i, :) = (rp_i(i, :) - p(cut_t))./4.9;
        rp_i(i, :) = 100.*(rp_i(i, :) - p(cut_t))./p(cut_t);
    end
    
    if (~ isempty(a)) && (~ isempty(p))
        sf = 4 .* pi .* a ./ (p .* p);
        rsf_i(i, startw:finishw) = sf(startr:finishr);
    end
    
    for jj = 1:(numel(bins)-1)
        ind = find((rl_i(i, :) >= bins(jj)) & (rl_i(i, :) < bins(jj+1)));
        binl_i(i, jj) = mean(rl_i(i, ind));
        bina_i(i, jj) = mean(ra_i(i, ind));
        binp_i(i, jj) = mean(rp_i(i, ind));
    end
    
    sh_l = circshift(l, [0 1]);  % We shift l (a row vector) to the right to compute velocity by substracting
    % current length from former length.
    sh_l(1) = l(1);
    v_vector = (l - sh_l) * um_per_pix ./ [repmat(t_resolution, [1 cut_t]) t_first(1) repmat(t_resolution, [1 numel(l)-(cut_t+1)])];
    maxv_i(i, :) = [max(v_vector((cut_t+1):finishr)) (l(cut_t) * um_per_pix)];
    [maxvind_i(i, 1) maxvind_i(i, 2)] = max(v_vector((cut_t+1):finishr));
    maxvind_i(i, 2) = maxvind_i(i, 2);
    
    angv_i(i, :) = [rmsnan(angle_values(cut_t, :), 2) maxv_i(i, 1)];
    angla_i(i, :) = [rmsnan(angle_values(cut_t, :), 2) max(rl_i_abs(i, :))];
    angln_i(i, :) = [rmsnan(angle_values(cut_t, :), 2) max(rl_i(i, :))];
    
    for jj = 1:size(angle_values, 2)
        singleangv_i(end+1, :) = [angle_values(cut_t, jj) maxv_i(i, 1)];
    end
    
    % If these are isolated interfaces, this condition will always be false.
    if exist(cat(2, directories{i}, 'analysis_cable.mat'), 'file')
        vars = load(cat(2, directories{i}, 'analysis_cable.mat'));
        
        nint_i(i) = vars.ninterfaces;
        cl_i(i) = vars.cable_length; % In microns already.
    else
        nint_i(i) = 1;
        cl_i(i) = l(cut_t) * um_per_pix;
    end
    
end

[ml_i, sl_i, el_i] = rmsnan(binl_i);
[ma_i, sa_i, ea_i] = rmsnan(bina_i);
[mp_i, sp_i, ep_i] = rmsnan(binp_i);




directories = dirs2;

rl_l = zeros(numel(directories), ntp).*nan;
rl_l_abs = zeros(numel(directories), ntp).*nan;
ra_l = zeros(numel(directories), ntp).*nan;
rp_l = zeros(numel(directories), ntp).*nan;
rsf_l = zeros(numel(directories), ntp).*nan;

binl_l = zeros(numel(directories), numel(bins) - 1) .* nan;
bina_l = zeros(numel(directories), numel(bins) - 1) .* nan;
binp_l = zeros(numel(directories), numel(bins) - 1) .* nan;

maxv_l = zeros(numel(directories), 2); % [vmax d0] matrix
maxvind_l = zeros(numel(directories), 2); %[vmax ind] matrix

tau_l = zeros(numel(directories), 2); % Model parameters.
taudamp_l = zeros(numel(directories), 3); % Model parameters.
R_l = [];
Rdamp_l = [];

angv_l = zeros(numel(directories), 2); % Average angle and max velocity.
angla_l = zeros(numel(directories), 2); % Average angle and max abs retraction.
angln_l = zeros(numel(directories), 2); % Average angle and max normalized retraction.
singleangv_l = zeros(0, 2); % Angle and max velocity.

nint_l = zeros(numel(directories, 1)); % Number of interfaces in the cable.
cl_l = zeros(numel(directories, 1)); % Cable length.


for i = 1:numel(directories)
    if exist(cat(2, directories{i}, 'analysis_local.mat'), 'file')
        load(cat(2, directories{i}, 'analysis_local'));
    elseif exist(cat(2, directories{i}, 'analysis_pairdist_onemin.mat'), 'file')
        load(cat(2, directories{i}, 'analysis_pairdist_onemin'));
    else
        disp(cat(2, 'Error when loading analysis results for ', directories{i}));
    end
    
    tau_l(i, :) = tau';
    if exist('taudamp', 'var')
        taudamp_l(i, :) = taudamp';
    end
    
    %startr = max(cut_t - (the_cut_t - 1), 1);
    startr = cut_t - (the_cut_t - 1);
    if startr > 0
        startw = 1;
    else
        startw = abs(startr-2);
        startr = 1;
    end
    
    finishr = cut_t + ntp - the_cut_t;
    finishw = startw + finishr - startr;
    
    
    
    if finishr > numel(l)
        finishr = cut_t + numel(l) - the_cut_t;
        finishw = startw + finishr - startr;
    end
    
    rl_l(i, startw:finishw) = l(startr:finishr);
    rl_l_abs(i, :) = (rl_l(i, :) - l(cut_t)) .* um_per_pix;
    R_l(end+1) = R;
    if exist('taudamp', 'var')
        Rdamp_l(end+1) = Rdamp;
    end
    rl_l(i, :) = 100.*(rl_l(i, :) - l(cut_t))./l(cut_t);
    
    if ~ isempty(a)
        ra_l(i, startw:finishw) = a(startr:finishr);
        %ra_l(i, :) = (ra_l(i, :) - a(cut_t))./(4.9^2.);
        ra_l(i, :) = 100.*(ra_l(i, :) - a(cut_t))./a(cut_t);
    end
    
    if ~ isempty(p)
        rp_l(i, startw:finishw) = p(startr:finishr);
        %rp_l(i, :) = (rp_l(i, :) - p(cut_t))./4.9;
        rp_l(i, :) = 100.*(rp_l(i, :) - p(cut_t))./p(cut_t);
    end
    
    if (~ isempty(a)) && (~ isempty(p))
        sf = 4 .* pi .* a ./ (p .* p);
        rsf_l(i, startw:finishw) = sf(startr:finishr);
    end
    
    for jj = 1:(numel(bins)-1)
        ind = find((rl_l(i, :) >= bins(jj)) & (rl_l(i, :) < bins(jj+1)));
        binl_l(i, jj) = mean(rl_l(i, ind));
        bina_l(i, jj) = mean(ra_l(i, ind));
        binp_l(i, jj) = mean(rp_l(i, ind));
    end
    
    sh_l = circshift(l, [0 1]);  % We shift l (a row vector) to the right to compute velocity by substracting
    % current length from former length.
    sh_l(1) = l(1);
    v_vector = (l - sh_l) * um_per_pix ./ [repmat(t_resolution, [1 cut_t]) t_first(2) repmat(t_resolution, [1 numel(l)-(cut_t+1)])];
    maxv_l(i, :) = [max(v_vector((cut_t+1):finishr)) (l(cut_t) * um_per_pix)];
    [maxvind_l(i, 1) maxvind_l(i, 2)] = max(v_vector((cut_t+1):finishr));
    maxvind_l(i, 2) = maxvind_l(i, 2);
    
    angv_l(i, :) = [rmsnan(angle_values(cut_t, :), 2) maxv_l(i, 1)];
    angla_l(i, :) = [rmsnan(angle_values(cut_t, :), 2) max(rl_l_abs(i, :))];
    angln_l(i, :) = [rmsnan(angle_values(cut_t, :), 2) max(rl_l(i, :))];
    
    for jj = 1:size(angle_values, 2)
        singleangv_l(end+1, :) = [angle_values(cut_t, jj) maxv_l(i, 1)];
    end
    
    if exist(cat(2, directories{i}, 'analysis_cable.mat'), 'file')
        vars = load(cat(2, directories{i}, 'analysis_cable.mat'));
        
        nint_l(i) = vars.ninterfaces;
        cl_l(i) = vars.cable_length; % In microns already.
    else
        nint_l(i) = 1;
        cl_l(i) = l(cut_t) * um_per_pix;
    end
end

[ml_l, sl_l, el_l] = rmsnan(binl_l);
[ma_l, sa_l, ea_l] = rmsnan(bina_l);
[mp_l, sp_l, ep_l] = rmsnan(binp_l);


if ~ isempty(dirs3)
    directories = dirs3;
    
    rl_h = zeros(numel(directories), ntp).*nan;
    rl_h_abs = zeros(numel(directories), ntp).*nan;
    ra_h = zeros(numel(directories), ntp).*nan;
    rp_h = zeros(numel(directories), ntp).*nan;
    rsf_h = zeros(numel(directories), ntp).*nan;
    
    binl_h = zeros(numel(directories), numel(bins) - 1) .* nan;
    bina_h = zeros(numel(directories), numel(bins) - 1) .* nan;
    binp_h = zeros(numel(directories), numel(bins) - 1) .* nan;
    
    maxv_h = zeros(numel(directories), 2); % [vmax d0] matrix
    maxvind_h = zeros(numel(directories), 2); %[vmax ind] matrix
    
    tau_h = zeros(numel(directories), 2); % Model parameters.
    taudamp_h = zeros(numel(directories), 3); % Model parameters.
    R_h = [];
    Rdamp_h = [];
    
    angv_h = zeros(numel(directories), 2); % Average angle and max velocity.
    angla_h = zeros(numel(directories), 2); % Average angle and max abs retraction.
    angln_h = zeros(numel(directories), 2); % Average angle and max normalized retraction.
    
    singleangv_h = zeros(0, 2); % Angle and max velocity.
    
    nint_h = zeros(numel(directories, 1)); % Number of interfaces in the cable.
    cl_h = zeros(numel(directories, 1)); % Cable length.
    
    for i = 1:numel(directories)
        if exist(cat(2, directories{i}, 'analysis_local.mat'), 'file')
            load(cat(2, directories{i}, 'analysis_local'));
        elseif exist(cat(2, directories{i}, 'analysis_pairdist_onemin.mat'), 'file')
            load(cat(2, directories{i}, 'analysis_pairdist_onemin'));
        else
            disp(cat(2, 'Error when loading analysis results for ', directories{i}));
        end
        
        tau_h(i, :) = tau';
        
        if exist('taudamp', 'var')
            taudamp_h(i, :) = taudamp';
        end
        
        %startr = max(cut_t - (the_cut_t - 1), 1);
        startr = cut_t - (the_cut_t - 1);
        if startr > 0
            startw = 1;
        else
            startw = abs(startr-2);
            startr = 1;
        end
        
        finishr = cut_t + ntp - the_cut_t;
        finishw = startw + finishr - startr;
        
        if finishr > numel(l)
            finishr = cut_t + numel(l) - the_cut_t;
            finishw = startw + finishr - startr;
        end
        
        rl_h(i, startw:finishw) = l(startr:finishr);
        rl_h_abs(i, :) = (rl_h(i, :) - l(cut_t)) .* um_per_pix;
        R_h(end+1) = R;
        if exist('taudamp', 'var')
            Rdamp_h(end+1) = Rdamp;
        end
        rl_h(i, :) = 100.*(rl_h(i, :) - l(cut_t))./l(cut_t); % relative edge length (% increase)
        
        if ~ isempty(a)
            ra_h(i, startw:finishw) = a(startr:finishr);
            %ra_h(i, :) = (ra_h(i, :) - a(cut_t))./(4.9^2.);
            ra_h(i, :) = 100.*(ra_h(i, :) - a(cut_t))./a(cut_t);
        end
        
        if ~ isempty(p)
            rp_h(i, startw:finishw) = p(startr:finishr);
            %rp_h(i, :) = (rp_h(i, :) - p(cut_t))./4.9;
            rp_h(i, :) = 100.*(rp_h(i, :) - p(cut_t))./p(cut_t);
        end
        
        if (~ isempty(a)) && (~ isempty(p))
            sf = 4 .* pi .* a ./ (p .* p);
            rsf_h(i, startw:finishw) = sf(startr:finishr);
        end
        
        for jj = 1:(numel(bins)-1)
            ind = find((rl_h(i, :) >= bins(jj)) & (rl_h(i, :) < bins(jj+1)));
            binl_h(i, jj) = mean(rl_h(i, ind));
            bina_h(i, jj) = mean(ra_h(i, ind));
            binp_h(i, jj) = mean(rp_h(i, ind));
        end
        
        sh_l = circshift(l, [0 1]);  % We shift l (a row vector) to the right to compute velocity by substracting
        % current length from former length.
        sh_l(1) = l(1);
        v_vector = (l - sh_l) * um_per_pix ./ [repmat(t_resolution, [1 cut_t]) t_first(3) repmat(t_resolution, [1 numel(l)-(cut_t+1)])];
        maxv_h(i, :) = [max(v_vector((cut_t+1):finishr)) (l(cut_t) * um_per_pix)];
        [maxvind_h(i, 1) maxvind_h(i, 2)] = max(v_vector((cut_t+1):finishr));
        maxvind_h(i, 2) = maxvind_h(i, 2);
        
        angv_h(i, :) = [rmsnan(angle_values(cut_t, :), 2) maxv_h(i, 1)];
        angla_h(i, :) = [rmsnan(angle_values(cut_t, :), 2) max(rl_h_abs(i, :))];
        angln_h(i, :) = [rmsnan(angle_values(cut_t, :), 2) max(rl_h(i, :))];
        
        for jj = 1:size(angle_values, 2)
            singleangv_h(end+1, :) = [angle_values(cut_t, jj) maxv_h(i, 1)];
        end
        
        % If these are isolated interfaces, this condition will always be false.
        if exist(cat(2, directories{i}, 'analysis_cable.mat'), 'file')
            vars = load(cat(2, directories{i}, 'analysis_cable.mat'));
            
            nint_h(i) = vars.ninterfaces;
            cl_h(i) = vars.cable_length; % In microns already.
        else
            nint_h(i) = 1;
            cl_h(i) = l(cut_t) * um_per_pix;
        end
    end
    
    [ml_h, sl_h, el_h] = rmsnan(binl_h);
    [ma_h, sa_h, ea_h] = rmsnan(bina_h);
    [mp_h, sp_h, ep_h] = rmsnan(binp_h);
end

if ~ noplots

    % Plot retracted distance (half the absolute distance).
    cut_t = the_cut_t;
    t = ((1:ntp)-cut_t).*t_resolution;
    Da = 0.0; % Length of the interface destroyed immediately by the laser (half the size of a single laser spot, 1-1.3um in Dresden, so this value should be 0.5-0.65). This is zero because you are measuring the retraction at the nodes, not at the tips created by the ablation.
    
    figure;
    
    [m, s, sem] = rmsnan(rl_i_abs);
    m(cut_t) = Da; % Length of the interface destroyed immediately by the laser.
    e = sem;
    
    h1 = errorbar(t((cut_t):(cut_t+9)), m((cut_t):(cut_t+9)), e((cut_t):(cut_t+9)), 'LineStyle', 'none', 'Marker', marker_str(1), 'MarkerFaceColor', color_str(1));
%     h1 = errorbar(t, m, e, 'LineStyle', 'none', 'Marker', marker_str(1), 'MarkerFaceColor', color_str(1));    
    set(h1, 'Color', color_str(1), 'LineWidth', 3);
    hold on;
    
    % Find tau.
    tau_isolated = nlinfit((0:t_resolution:(9*t_resolution)), m((cut_t):(cut_t+9)), @viscoelastic_cable_model, [1 3]);
%     tau_isolated = nlinfit((0:t_resolution:((length(m)-1)*t_resolution)), m, @viscoelastic_cable_model, [1 3]);
    plot((0:t_resolution:((length(m)-1)*t_resolution)), viscoelastic_cable_model(tau_isolated, (0:t_resolution:((length(m)-1)*t_resolution))), 'Color', color_str(1), 'LineWidth', 2);
    [r_isolated p_isolated] = corrcoef(m((cut_t):(cut_t+9)), viscoelastic_cable_model(tau_isolated, (0:t_resolution:(9*t_resolution))));
    
    [m, s, sem] = rmsnan(rl_l_abs);
    m(cut_t) = Da; % Length of the interface destroyed immediately by the laser.
    e = sem;
    
    h2 = errorbar(t((cut_t):(cut_t+9)), m((cut_t):(cut_t+9)), e((cut_t):(cut_t+9)), 'LineStyle', 'none', 'Marker', marker_str(2), 'MarkerFaceColor', color_str(2));
%     h2 = errorbar(t, m, e, 'LineStyle', 'none', 'Marker', marker_str(2), 'MarkerFaceColor', color_str(2));
    set(h2, 'Color', color_str(2), 'LineWidth', 3);
    
    % Find tau.
    tau_linked = nlinfit((0:t_resolution:(9*t_resolution)), m((cut_t):(cut_t+9)), @viscoelastic_cable_model, [1 4]);
%     tau_linked = nlinfit((0:t_resolution:((length(m)-1)*t_resolution)), m, @viscoelastic_cable_model, [1 4]);
    plot((0:t_resolution:((length(m)-1)*t_resolution)), viscoelastic_cable_model(tau_linked, (0:t_resolution:((length(m)-1)*t_resolution))), 'Color', color_str(2), 'LineWidth', 2);
    [r_linked p_linked] = corrcoef(m((cut_t):(cut_t+9)), viscoelastic_cable_model(tau_linked, (0:t_resolution:(9*t_resolution))));
    
    h3 = [];
    
    if ~isempty(dirs3)
        [m, s, sem] = rmsnan(rl_h_abs);
        m(cut_t) = Da; % Length of the interface destroyed immediately by the laser.
        e = sem;
        
        h3 = errorbar(t((cut_t):(cut_t+10)), m((cut_t):(cut_t+10)), e((cut_t):(cut_t+10)), 'LineStyle', 'none', 'Marker', marker_str(3), 'MarkerFaceColor', color_str(3));
        set(h3, 'Color', color_str(3), 'LineWidth', 3);
        
        % Find tau.
        tau_horizontal = nlinfit((0:t_resolution:(10*t_resolution)), m((cut_t):(cut_t+10)), @viscoelastic_cable_model, [1 1]);
        plot((0:t_resolution:(10*t_resolution)), viscoelastic_cable_model(tau_horizontal, (0:t_resolution:(10*t_resolution))), 'Color', color_str(3), 'LineWidth', 2);
        [r_horizontal p_horizontal]= corrcoef(m((cut_t):(cut_t+10)), viscoelastic_cable_model(tau_horizontal, (0:t_resolution:(10*t_resolution))));
    end
    
    
    %h = legend([h1 h2 h3], legend_str, 'Location', 'Best');
    set(gca, 'LineWidth', 2);
    set(gca, 'XTick', [0:15:60], 'YTick', [0:1:20], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold');
    xlabel('Time after the cut (seconds)', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('Distance retracted (\mum)', 'FontSize', 32, 'FontWeight', 'bold');
    title(title_label);
    
    if (save_flag)
        savefig(cat(2, 'distanceretractedmicrons_time_mean', rstrtrim(title_label)), gcf, 'jpeg');
        saveas(gcf, cat(2, 'distanceretractedmicrons_time_mean', rstrtrim(title_label)), 'fig');
    end
    
    
    
    % Take only those values that represent good fits of the data
    % (Rmovie >= R_min).
    bad_fits_i = find(R_i < R_min);
    
    if ~ isempty(find(tau_i(:, 1) > T_max))
        bad_fits_i = unique(cat(2, bad_fits_i, find(tau_i(:, 1) > T_max)'));
    end
    tau_i2 = tau_i;
    taudamp_i2 = taudamp_i;
    tau_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(tau_i, 2));
    taudamp_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(taudamp_i, 2));
    [m_i s_i e_i] = rmsnan(tau_i2);
    [mdamp_i sdamp_i edamp_i] = rmsnan(taudamp_i2);
    
    bad_fits_l = find(R_l < R_min);
    if ~ isempty(find(tau_l(:, 1) > T_max))
        bad_fits_l = unique(cat(2, bad_fits_l, find(tau_l(:, 1) > T_max)'));
    end
    tau_l2 = tau_l;
    taudamp_l2 = taudamp_l;
    tau_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(tau_l, 2));
    taudamp_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(taudamp_l, 2));
    [m_l s_l e_l] = rmsnan(tau_l2);
    [mdamp_l sdamp_l edamp_l] = rmsnan(taudamp_l2);
    
    if ~ isempty(dirs3)
        bad_fits_h = find(R_h < R_min);
        if ~ isempty(find(tau_h(:, 1) > T_max))
            bad_fits_h = unique(cat(2, bad_fits_h, find(tau_h(:, 1) > T_max)'));
        end
        tau_h2 = tau_h;
        taudamp_h2 = taudamp_h;
        tau_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(tau_h, 2));
        taudamp_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(taudamp_h, 2));
        [m_h s_h e_h] = rmsnan(tau_h2);
        [mdamp_h sdamp_h edamp_h] = rmsnan(taudamp_h2);
        
        x = 1:3;
        y = [m_h(1) m_i(1) m_l(1)];
        e = [e_h(1) e_i(1) e_l(1)];
        
        ydamp = [mdamp_h(3) mdamp_i(3) mdamp_l(3)];
        edamp = [edamp_h(3) edamp_i(3) edamp_l(3)];
        
    else
        x = 1:2;
        y = [m_i(1) m_l(1)];
        e = [e_i(1) e_l(1)];
        
        ydamp = [mdamp_i(1) mdamp_l(1)];
        edamp = [edamp_i(1) edamp_l(1)];
    end
    
    % Take only those values that represent good fits of the data
    % (Rmovie >= R_min).
    bad_fits_i = find(R_i < R_min);
    
    if ~ isempty(find(tau_i(:, 1) > T_max))
        bad_fits_i = unique(cat(2, bad_fits_i, find(tau_i(:, 1) > T_max)'));
    end
    tau_i2 = tau_i;
    tau_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(tau_i, 2));
    [m_i s_i e_i] = rmsnan(tau_i2);
    
    bad_fits_l = find(R_l < R_min);
    if ~ isempty(find(tau_l(:, 1) > T_max))
        bad_fits_l = unique(cat(2, bad_fits_l, find(tau_l(:, 1) > T_max)'));
    end
    tau_l2 = tau_l;
    tau_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(tau_l, 2));
    [m_l s_l e_l] = rmsnan(tau_l2);
    
    if ~ isempty(dirs3)
        bad_fits_h = find(R_h < R_min);
        if ~ isempty(find(tau_h(:, 1) > T_max))
            bad_fits_h = unique(cat(2, bad_fits_h, find(tau_h(:, 1) > T_max)'));
        end
        tau_h2 = tau_h;
        tau_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(tau_h, 2));
        [m_h s_h e_h] = rmsnan(tau_h2);
        
        x = 1:3;
        y = [m_h(2) m_i(2) m_l(2)];
        e = [e_h(2) e_i(2) e_l(2)];
        
    else
        x = 1:2;
        y = [m_i(2) m_l(2)];
        e = [e_i(2) e_l(2)];
    end

    
    % Plot both model parameters in the same graph.
    % Plot model parameters (not those computed from the average displacements,
    % but those computed from each movie (mean and std).
    % Take only those values that represent good fits of the data
    % (Rmovie >= R_min).
    
    bad_fits_i = find(R_i < R_min);
    
    if ~ isempty(find(tau_i(:, 1) > T_max))
        bad_fits_i = unique(cat(2, bad_fits_i, find(tau_i(:, 1) > T_max)'));
    end
    tau_i2 = tau_i;
    taudamp_i2 = taudamp_i;
    tau_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(tau_i, 2));
    taudamp_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(taudamp_i, 2));
    [m_i s_i e_i] = rmsnan(tau_i2);
    [mdamp_i sdamp_i edamp_i] = rmsnan(taudamp_i2);
    
    bad_fits_l = find(R_l < R_min);
    if ~ isempty(find(tau_l(:, 1) > T_max))
        bad_fits_l = unique(cat(2, bad_fits_l, find(tau_l(:, 1) > T_max)'));
    end
    tau_l2 = tau_l;
    taudamp_l2 = taudamp_l;
    tau_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(tau_l, 2));
    taudamp_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(taudamp_l, 2));
    [m_l s_l e_l] = rmsnan(tau_l2);
    [mdamp_l sdamp_l edamp_l] = rmsnan(taudamp_l2);
    
    if ~ isempty(dirs3)
        bad_fits_h = find(R_h < R_min);
        if ~ isempty(find(tau_h(:, 1) > T_max))
            bad_fits_h = unique(cat(2, bad_fits_h, find(tau_h(:, 1) > T_max)'));
        end
        tau_h2 = tau_h;
        taudamp_h2 = taudamp_h;
        tau_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(tau_h, 2));
        taudamp_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(taudamp_h, 2));
        [m_h s_h e_h] = rmsnan(tau_h2);
        [mdamp_h sdamp_h edamp_h] = rmsnan(taudamp_h2);
        
        x = 1:3;
        y = [m_h(1) m_i(1) m_l(1)];
        e = [e_h(1) e_i(1) e_l(1)];
        
        ydamp = [mdamp_h(1) mdamp_i(1) mdamp_l(1)];
        edamp = [edamp_h(1) edamp_i(1) edamp_l(1)];
        
    else
        x = 1:2;
        y = [m_i(1) m_l(1)];
        e = [e_i(1) e_l(1)];
        
        ydamp = [mdamp_i(1) mdamp_l(1)];
        edamp = [edamp_i(1) edamp_l(1)];
    end
    
   
    % Take only those values that represent good fits of the data
    % (Rmovie >= R_min).
    bad_fits_i = find(R_i < R_min);
    
    if ~ isempty(find(tau_i(:, 1) > T_max))
        bad_fits_i = unique(cat(2, bad_fits_i, find(tau_i(:, 1) > T_max)'));
    end
    tau_i2 = tau_i;
    taudamp_i2 = taudamp_i;
    tau_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(tau_i, 2));
    taudamp_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(taudamp_i, 2));
    [m_i s_i e_i] = rmsnan(tau_i2);
    [mdamp_i sdamp_i edamp_i] = rmsnan(taudamp_i2);
    
    bad_fits_l = find(R_l < R_min);
    if ~ isempty(find(tau_l(:, 1) > T_max))
        bad_fits_l = unique(cat(2, bad_fits_l, find(tau_l(:, 1) > T_max)'));
    end
    tau_l2 = tau_l;
    taudamp_l2 = taudamp_l;
    tau_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(tau_l, 2));
    taudamp_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(taudamp_l, 2));
    [m_l s_l e_l] = rmsnan(tau_l2);
    [mdamp_l sdamp_l edamp_l] = rmsnan(taudamp_l2);
    
    if ~ isempty(dirs3)
        bad_fits_h = find(R_h < R_min);
        if ~ isempty(find(tau_h(:, 1) > T_max))
            bad_fits_h = unique(cat(2, bad_fits_h, find(tau_h(:, 1) > T_max)'));
        end
        tau_h2 = tau_h;
        taudamp_h2 = taudamp_h;
        tau_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(tau_h, 2));
        taudamp_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(taudamp_h, 2));
        [m_h s_h e_h] = rmsnan(tau_h2);
        [mdamp_h sdamp_h edamp_h] = rmsnan(taudamp_h2);
        
        x1 = [1 6 11];
        x2 = [3 8 13];
        y1 = [m_h(1) m_i(1) mdamp_l(1)];
        e1 = [e_h(1) e_i(1) edamp_l(1)];
        y2 = [m_h(2) m_i(2) m_l(2)];
        e2 = [e_h(2) e_i(2) e_l(2)];
        
    else
        x1 = [1 6];
        x2 = [3 8];
        y1 = [m_i(1) m_l(1)];
        e1 = [e_i(1) e_l(1)];
        y2 = [m_i(2) m_l(2)];
        e2 = [e_i(2) e_l(2)];
    end
    
    % Compute and plot the peak velocities based on the two different models.
    bad_fits_i = find(R_i < R_min);
    
    if ~ isempty(find(tau_i(:, 1) > T_max))
        bad_fits_i = unique(cat(2, bad_fits_i, find(tau_i(:, 1) > T_max)'));
    end
    tau_i2 = tau_i;
    taudamp_i2 = taudamp_i;
    tau_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(tau_i, 2));
    taudamp_i2(bad_fits_i, :) = nan .* zeros(numel(bad_fits_i), size(taudamp_i, 2));
    [m_i s_i e_i] = rmsnan(tau_i2);
    [mdamp_i sdamp_i edamp_i] = rmsnan(taudamp_i2);
    
    bad_fits_l = find(R_l < R_min);
    if ~ isempty(find(tau_l(:, 1) > T_max))
        bad_fits_l = unique(cat(2, bad_fits_l, find(tau_l(:, 1) > T_max)'));
    end
    tau_l2 = tau_l;
    taudamp_l2 = taudamp_l;
    tau_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(tau_l, 2));
    taudamp_l2(bad_fits_l, :) = nan .* zeros(numel(bad_fits_l), size(taudamp_l, 2));
    [m_l s_l e_l] = rmsnan(tau_l2);
    [mdamp_l sdamp_l edamp_l] = rmsnan(taudamp_l2);
    
    if ~ isempty(dirs3)
        bad_fits_h = find(R_h < R_min);
        if ~ isempty(find(tau_h(:, 1) > T_max))
            bad_fits_h = unique(cat(2, bad_fits_h, find(tau_h(:, 1) > T_max)'));
        end
        tau_h2 = tau_h;
        taudamp_h2 = taudamp_h;
        tau_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(tau_h, 2));
        taudamp_h2(bad_fits_h, :) = nan .* zeros(numel(bad_fits_h), size(taudamp_h, 2));
        [m_h s_h e_h] = rmsnan(tau_h2);
        [mdamp_h sdamp_h edamp_h] = rmsnan(taudamp_h2);
        
        
    end
    
    maxv_model1_i = nan .* (zeros(size(tau_i2, 1), 1));
    maxv_model1_l = nan .* (zeros(size(tau_l2, 1), 1));
    maxv_model2_i = nan .* (zeros(size(tau_i2, 1), 1));
    maxv_model2_l = nan .* (zeros(size(tau_l2, 1), 1));
    maxl_model1_i = nan .* (zeros(size(tau_i2, 1), 1));
    maxl_model1_l = nan .* (zeros(size(tau_l2, 1), 1));
    maxl_model2_i = nan .* (zeros(size(tau_i2, 1), 1));
    maxl_model2_l = nan .* (zeros(size(tau_l2, 1), 1));
    tau_model1_i = nan .* (zeros(size(tau_i2, 1), 1));
    tau_model1_l = nan .* (zeros(size(tau_l2, 1), 1));
    tau_model2_i = nan .* (zeros(size(tau_i2, 1), 1));
    tau_model2_l = nan .* (zeros(size(tau_l2, 1), 1));
    tau_i = nan .* (zeros(size(tau_i2, 1), 1));
    tau_l = nan .* (zeros(size(tau_l2, 1), 1));
    
    for ii = 1:size(tau_i2, 1)
        if ~isnan(tau_i2(ii, 1))
            maxv_model1_i(ii) = viscoelastic_cable_model(tau_i2(ii, :), t_resolution) ./ t_resolution;
            maxv_model2_i(ii) = viscoelastic_cable_model_withextradamper(taudamp_i2(ii, :), t_resolution) ./ t_resolution;
            
            maxl_model1_i(ii) = viscoelastic_cable_model(tau_i2(ii, :), 60);
            maxl_model2_i(ii) = viscoelastic_cable_model_withextradamper(taudamp_i2(ii, :), 60);
            
            dtau = (1-exp(-1)) .* rl_i_abs(ii, the_cut_t+9);
            ind = find(rl_i_abs(ii, the_cut_t:end)-dtau > 0);
            ind = [ind(1)-1 ind(1)];
            tau_i(ii) = (ind(1) + (dtau-rl_i_abs(ii, the_cut_t + ind(1) - 1))/(rl_i_abs(ii, the_cut_t + ind(2) - 1)-rl_i_abs(ii, the_cut_t + ind(1) - 1))) * t_resolution;
            
            %Analytical solution to model 1.
            tau_model1_i(ii) = -1 * tau_i2(ii, 1) * log(1 - dtau/tau_i2(ii, 2));
            % Numerical solution to model 1.
            %                  vect = viscoelastic_cable_model(tau_i2(ii, :), (0:0.01:60));
            %                  [dummy ind] = min(abs(vect-dtau));
            %                  tau_model1_i(ii) = (ind - 1) * 0.01;
            
            % Numerical solution to model 2.
            vect = viscoelastic_cable_model_withextradamper(taudamp_i2(ii, :), (0:0.01:60));
            [dummy ind] = min(abs(vect-dtau));
            tau_model2_i(ii) = (ind - 1) * 0.01;
        end
    end
    
    for ii = 1:size(tau_l2, 1)
        if ~isnan(tau_l2(ii, 1))
            maxv_model1_l(ii) = viscoelastic_cable_model(tau_l2(ii, :), t_resolution) ./ t_resolution;
            maxv_model2_l(ii) = viscoelastic_cable_model_withextradamper(taudamp_l2(ii, :), t_resolution) ./ t_resolution;
            
            maxl_model1_l(ii) = viscoelastic_cable_model(tau_l2(ii, :), 60);
            maxl_model2_l(ii) = viscoelastic_cable_model_withextradamper(taudamp_l2(ii, :), 60);
            
            dtau = (1-exp(-1)) .* rl_l_abs(ii, the_cut_t+9);
            ind = find(rl_l_abs(ii, the_cut_t:end)-dtau > 0);
            ind = [ind(1)-1 ind(1)];
            tau_l(ii) = (ind(1) + (dtau-rl_l_abs(ii, the_cut_t + ind(1) - 1))/(rl_l_abs(ii, the_cut_t + ind(2) - 1)-rl_l_abs(ii, the_cut_t + ind(1) - 1))) * t_resolution;
            
            %Analytical solution to model 1.
            tau_model1_l(ii) = -1 * tau_l2(ii, 1) * log(1 - dtau/tau_l2(ii, 2));
            % Numerical solution to model 1.
            %                  vect = viscoelastic_cable_model(tau_l2(ii, :), (0:0.01:60));
            %                  [dummy ind] = min(abs(vect-dtau));
            %                  tau_model1_l(ii) = (ind - 1) * 0.01;
            
            % Numerical solution to model 2.
            vect = viscoelastic_cable_model_withextradamper(taudamp_l2(ii, :), (0:0.01:60));
            [dummy ind] = min(abs(vect-dtau));
            tau_model2_l(ii) = (ind - 1) * 0.01;
            
        end
    end
    
    
    if ~ isempty(dirs3)
        maxv_model1_h = nan .* (zeros(size(tau_h2, 1), 1));
        maxv_model2_h = nan .* (zeros(size(tau_h2, 1), 1));
        maxl_model1_h = nan .* (zeros(size(tau_h2, 1), 1));
        maxl_model2_h = nan .* (zeros(size(tau_h2, 1), 1));
        tau_model1_h = nan .* (zeros(size(tau_i2, 1), 1));
        tau_model2_h = nan .* (zeros(size(tau_i2, 1), 1));
        tau_h = nan .* (zeros(size(tau_l2, 1), 1));
        for ii = 1:size(tau_h2, 1)
            if ~isnan(tau_h2(ii, 1))
                maxv_model1_h(ii) = viscoelastic_cable_model(tau_h2(ii, :), t_resolution) ./ t_resolution;
                maxv_model2_h(ii) = viscoelastic_cable_model_withextradamper(taudamp_h2(ii, :), t_resolution) ./ t_resolution;
                
                maxl_model1_h(ii) = viscoelastic_cable_model(tau_h2(ii, :), 60);
                maxl_model2_h(ii) = viscoelastic_cable_model_withextradamper(taudamp_h2(ii, :), 60);
                
                dtau = (1-exp(-1)) .* rl_h_abs(ii, the_cut_t+10);
                ind = find(rl_h_abs(ii, the_cut_t:end)-dtau > 0);
                ind = [ind(1)-1 ind(1)];
                tau_h(ii) = (ind(1) + (dtau-rl_h_abs(ii, the_cut_t + ind(1) - 1))/(rl_h_abs(ii, the_cut_t + ind(2) - 1)-rl_h_abs(ii, the_cut_t + ind(1) - 1))) * t_resolution;
                
                %Analytical solution to model 1.
                tau_model1_h(ii) = -1 * tau_h2(ii, 1) * log(1 - dtau/tau_h2(ii, 2));
                % Numerical solution to model 1.
                %                  vect = viscoelastic_cable_model(tau_h2(ii, :), (0:0.01:60));
                %                  [dummy ind] = min(abs(vect-dtau));
                %                  tau_model1_h(ii) = (ind - 1) * 0.01;
                
                % Numerical solution to model 2.
                vect = viscoelastic_cable_model_withextradamper(taudamp_h2(ii, :), (0:0.01:60));
                [dummy ind] = min(abs(vect-dtau));
                tau_model2_h(ii) = (ind - 1) * 0.01;
            end
        end
        x = 1:3;
        
        % Relaxation time.
        [y1(1) dummy e1(1)] = rmsnan(maxv_model1_h);
        [y1(2) dummy e1(2)] = rmsnan(maxv_model1_i);
        [y1(3) dummy e1(3)] = rmsnan(maxv_model1_l);
        
        [y2(1) dummy e2(1)] = rmsnan(maxv_model2_h);
        [y2(2) dummy e2(2)] = rmsnan(maxv_model2_i);
        [y2(3) dummy e2(3)] = rmsnan(maxv_model2_l);
        
        [y3(1) dummy e3(1)] = rmsnan(maxv_h(:, 1));
        [y3(2) dummy e3(2)] = rmsnan(maxv_i(:, 1));
        [y3(3) dummy e3(3)] = rmsnan(maxv_l(:, 1));
        
        
    else
        x = 1:2;
        
        % Model 1 estimates of peak velocity.
        [y1(1) dummy e1(1)] = rmsnan(maxv_model1_i);
        [y1(2) dummy e1(2)] = rmsnan(maxv_model1_l);
        
        % Model 2 estimates.
        [y2(1) dummy e2(1)] = rmsnan(maxv_model2_i);
        [y2(2) dummy e2(2)] = rmsnan(maxv_model2_l);
        
        % Real data.
        [y3(1) dummy e3(1)] = rmsnan(maxv_i(:, 1));
        [y3(2) dummy e3(2)] = rmsnan(maxv_l(:, 1));
        
    end
    
 
    %% viscosity-to-elasticity ratio
    figure;
    var1 = tau_i2(find(~isnan(tau_i2(:,1))),1)
    var2 = tau_l2(find(~isnan(tau_l2(:,1))),1)
    var3 = tau_h2(find(~isnan(tau_h2(:,1))),1)
    
    boxplotMWstd(0.5:1.5,var1,3,'bo',1);
    hold on;
    boxplotMWstd(2.0:3.0,var2,3,'ro',1);
    hold on;
    boxplotMWstd(3.5:4.5,var3,3,'go',1);
    
    ylim([0 40]);
    xlim([0 5]);
    box off;
    set(gca, 'LineWidth', 2);
    set(gca, 'XTick', [1.0 2.5 4.0], 'XTickLabel', legend_str, 'YTickLabel', [0:10:40], 'YTick', [0:10:40], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
    xlabel('', 'FontSize', 32, 'FontWeight', 'bold');
    ylabel('relaxation time, {\it\tau_r} (s)', 'FontSize', 32, 'FontWeight', 'bold');
    
    
    
end 




