%% read in results from polygon analysis and calculate/plot closure rate constant
clear all

%% Write the path for each directory to analyze. To obtain the "area_curve.mat" matrix you need to run "polygon_analysis_script_wh_teresa.m" setting: save_area=1; and extension='mat'; 

con = {  
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo7/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201002/Embryo1/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201002/Embryo3/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo1/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo2/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo3/area_curve.mat',...
    };

Rap1DN = {
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo2/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo3/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20200930/Embryo4/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo1/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo2/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201001/Embryo3/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo6/area_curve.mat',...
    '/Users/katy/Documents/Data/2020 daGal4 Rap1DN sqhGFP EcadtdTo Wound/20201007/Embryo5/area_curve.mat',...
    };

dirs1 = con;
dirs2 = Rap1DN;
legend_str = {'control','Rap1DN'};

%% Parameters
At=30/60;% temporal resolution in minutes
t=0:At:100;

%% Initialize variables
expon_1 = zeros(1,length(dirs1));
expon_2 = zeros(1,length(dirs2));
max_areas_1 = zeros(1,length(dirs1));
max_areas_2 = zeros(1,length(dirs2));

%% Fit exponential dirs1
for ii = 1:numel(dirs1)
    
    area_normalized=[];
    Area=[];
    Area0=[];
    
    load(dirs1{ii})
    Area0=area_absolute;
    t0=[];
    [m0 t0]=max(Area0); % calculate tp of max. expansion
    
    max_areas_1(ii) = m0; % save max area
    
    Area=Area0(:,t0:end); % part of the area that will be fitted

    params_exp=[];
    params_exp=nlinfit(t(1:length(Area)),Area, @exponential_function,[Area(1) 0.1]); % fitting
    
    % Plot fitted curves
    figure()
    plot(t(1:length(Area)), Area,'b','LineWidth',2); 
    hold on
    time_exp=t(1:length(Area));
    yy=params_exp(1)*exp(-(params_exp(2)*time_exp));
    plot(time_exp,yy,'r','LineWidth',2)
    
    expon_1(ii)=params_exp(2);
    legend(num2str(expon_1(ii)))
    R=[];
    P=[];
    [R P]=corrcoef(Area(~isnan(Area)),yy(~isnan(Area)));
    correlation_1(ii)=R(1,2);
    pvalue_1(ii)=P(1,2);
    exp_fit = params_exp(2);
    save(dirs1{ii},'exp_fit','m0','-append')
    
end
%% Fit exponential dirs2
for ii = 1:numel(dirs2)
    
    area_normalized=[];
    Area=[];
    Area0=[];
    
    load(dirs2{ii})
    Area0=area_absolute;
    t0=[];
    [m0 t0]=max(Area0); % calculate tp of max. expansion
    
    max_areas_2(ii) = m0;
            
    Area=Area0(:,t0:end); % part of the area that will be fitted

    params_exp=[];
    params_exp=nlinfit(t(1:length(Area)),Area, @exponential_function,[Area(1) 0.1]); % fitting
    
    % Plot fitted curves
    figure()
    plot(t(1:length(Area)), Area,'b','LineWidth',2); 
    hold on
    time_exp=t(1:length(Area));
    yy=params_exp(1)*exp(-(params_exp(2)*time_exp));
    plot(time_exp,yy,'r','LineWidth',2)
    
    expon_2(ii)=params_exp(2);
    legend(num2str(expon_2(ii)))
    R=[];
    P=[];
    [R P]=corrcoef(Area(~isnan(Area)),yy(~isnan(Area)));
    correlation_2(ii)=R(1,2);
    pvalue_2(ii)=P(1,2);
    exp_fit = params_exp(2);
    save(dirs2{ii},'exp_fit','m0','-append')
end

%% Calculate closure rate constant means
[m_1, s_1, e_1] = rmsnan((expon_1'));
[m_2, s_2, e_2] = rmsnan((expon_2'));
%% Plot closure rate constant 
xrange=[0.1 0.6];
xt=mean([0.1 0.6]);
figure
hold on
boxplotMWstd([0.5 1.5],expon_1,3,'bo',3)
boxplotMWstd([2 3],expon_2,3,'ro',3)
set(gca, 'XTick', [1 2.5], 'XGrid', 'off', 'YGrid', 'off', 'FontSize', 32, 'FontWeight', 'bold', 'LineWidth', 2);
ylabel( 'closure rate constant (1/min)', 'FontSize', 32)
xlim([0 3.5]);
set(gca,'XTickLabel',legend_str)
%ylim([0 0.15]);
%% Run statistics
[zT, T] = rod_mannwhitney(expon_1, expon_2)