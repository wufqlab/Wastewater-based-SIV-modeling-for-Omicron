%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot/fit to Omicron and Delta data from the ADOO site
% This uses a 2-strain model; however, the focus is on the emergence of
% Omicron, so only the second strain is modeled in this file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global F
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex','defaultAxesFontSize',16) 
format long
rng('default')
bl = '#0072BD';
br = '#D95319';
%% load
load("data.mat")

% ABOO data
V_AD_O = ADflows.*ADrna_Omicron;
V_AD_D = ADflows.*ADrna_Delta;

% Cut-off for fitting/forecasting
split = 25; %other split for test: 28 27 26 25
F = 5; %100/F = the susceptible fraction of total population (F: 5 10 20 100)

allData = [ADdays V_AD_D V_AD_O];
allData = allData(split:end,:);

ADdays = allData(:,1);
V_AD_D_new = allData(:,2);
V_AD_O_new = allData(:,3);

%% Curve-fitting
%Note, we only model 2nd strain.
options = optimset('TolX',1e-8,'TolFun',1e-8,'TolCon',1E-8, ...
    'MaxIter',50000,'MaxFunEvals',100000,'display','off'); 

%     beta1 beta2  d1   d2   alpha  
lb = [1E-12 1E-12  0    0    1e9   ];
ub = [1E-4  1E-4   0    0    1e11  ]; 
p0 = [1E-8  1E-8   0    0    1e10  ]; 

% require I20 < I10 and beta2 > beta 1
A = [1 -1 0 0 0];
b = [0];

%% try Global Search/MultiStart
problem = createOptimProblem('fmincon','x0',p0,...
    'objective',@(param)obj_fun(param,ADdays,V_AD_D_new,V_AD_O_new),'lb',lb,'ub',ub,'Aineq',A,'bineq',b);

% multistart
ms = MultiStart("Display","iter");
[best_params,SSE] = run(ms,problem,100);

parameter = ["beta1";"beta2";"d1";"d2";"alpha";"SSE"];
disp(SSE)
estimated_val = [best_params';SSE];
t = table(parameter,estimated_val);

beta_report  = best_params(2);
alpha_report = best_params(5);

N0 = 1.3*10^6/F;
R01 = N0*best_params(1)*8;
R02 = N0*best_params(2)*8;

%% simulate
%Note, we only model 2nd strain.
I10 = 0;
I20 = V_AD_O_new(1)/best_params(5);
R0 = 0;
S0 = N0 - (I10 + I20 + R0);               
ICs  = [S0 I10 I20 R0 V_AD_D_new(1) V_AD_O_new(1) I10+I20];

tspan = ADdays(1):1:ADdays(end);
[~,Y] = ode45(@SIIR,tspan,ICs,[],best_params);

% run backwards
tspan = ADdays(1):-1:1;
[T,Yr] = ode45(@SIIR,tspan,ICs,[],best_params);
ind = find(flip(Yr(:,3)) < 1, 1, 'last');

%% Bootstrap (B=200) on residuals to get parameter CIs and prediction band
B = 200;

% Baseline model residuals (log10 space) at observed sample dates
tspan_fwd_bs = ADdays(1):1:ADdays(end);
[~,Y0_bs] = ode45(@SIIR,tspan_fwd_bs,ICs,[],best_params);

ind_bs = ADdays - ADdays(1);
ind_bs = ind_bs(2:end);

Omicron0_bs = Y0_bs(:,6);
daily0_bs = diff(Omicron0_bs);
daily0_bs = daily0_bs(ind_bs);

mu_log_bs  = log10(abs(daily0_bs));
obs_log_bs = log10(V_AD_O_new(2:end));
resid0_bs  = obs_log_bs - mu_log_bs;
n_bs = numel(resid0_bs);

% Prepare containers
param_boot = nan(B, numel(best_params));
SSE_boot   = nan(B, 1);

% Determine full plotted daily curve length using current simulation 'Y' & 'Yr'
daily_full = diff([flip(Yr(2:end,6))' Y(:,6)']);
L_full = numel(daily_full);
curve_boot_full = nan(B, L_full);

for bb = 1:B
    bb %this count the bootstrap
    % residual resampling
    res_b = resid0_bs( randi(n_bs, n_bs, 1) );
    ystar_log = mu_log_bs + res_b; % synthetic log10 daily Omicron
    
    data2_star = V_AD_O_new;
    data2_star(2:end) = max(eps, 10.^ystar_log);
    
    % Refit on bootstrap sample
    problem_b = createOptimProblem('fmincon','x0',p0,...
        'objective',@(param)obj_fun(param,ADdays,V_AD_D_new,data2_star),'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
    ms_b = MultiStart("Display","off");
    [param_b, SSE_b] = run(ms_b, problem_b, 100);
    param_boot(bb,:) = param_b;
    SSE_boot(bb) = SSE_b;
    
    % Simulate full curve for plotting band (backward + forward, matching main plot)
    I10_b = 0;
    I20_b = data2_star(1)/param_b(5);
    R0_b  = 0;
    S0_b  = N0 - (I10_b + I20_b + R0_b);
    ICs_b = [S0_b I10_b I20_b R0_b V_AD_D_new(1) data2_star(1) I10_b+I20_b];
    
    tspan_f = ADdays(1):1:ADdays(end);
    [~,Yb_f] = ode45(@SIIR,tspan_f,ICs_b,[],param_b);
    tspan_r = ADdays(1):-1:1;
    [~,Yb_r] = ode45(@SIIR,tspan_r,ICs_b,[],param_b);
    
    daily_b_full = diff([flip(Yb_r(2:end,6))' Yb_f(:,6)']);
    curve_boot_full(bb, :) = daily_b_full(:)';
end

% Summarize parameter distribution
param_names = {'beta1','beta2','d1','d2','alpha'};
param_mean   = mean(param_boot, 1, 'omitnan');
param_median = median(param_boot, 1, 'omitnan');
param_ci_lo  = prctile(param_boot, 2.5, 1);
param_ci_hi  = prctile(param_boot, 97.5, 1);

T_boot = table(param_names', param_mean', param_median', param_ci_lo', param_ci_hi', ...
    'VariableNames', {'parameter','mean','median','ci_lo','ci_hi'});
%% plot
time = ABflowdates(split):ABflowdates(end);

time2 = ABflowdates(1):ABflowdates(end);
figure(); box on; hold on;

    dailyOmicron = diff([flip(Yr(2:end,6))' Y(:,6)']);
    plot(time2(end-split-1:end),dailyOmicron(end-split-1:end),'LineWidth',3,'color',br);
    plot(ABflowdates(split-1:end),V_AD_O(split-1:end),'.','LineWidth',3,'MarkerSize',20,'color',br);

    
    % Bootstrap ribbons & mean/median curves
    x_plot = time2(end-split-1:end);
    seg_len = numel(x_plot);
    L_full = size(curve_boot_full, 2);
    idx = (L_full - seg_len + 1): L_full;
    seg = curve_boot_full(:,idx);
    band_lo = prctile(seg,  2.5, 1);
    band_hi = prctile(seg, 97.5, 1);
    curve_mean_bs   = mean(seg, 1, 'omitnan');
    curve_median_bs = median(seg, 1, 'omitnan');
    fill([x_plot fliplr(x_plot)], [band_lo fliplr(band_hi)], [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(x_plot, curve_mean_bs, '--', 'LineWidth', 2, 'color', br);
    plot(x_plot, curve_median_bs, ':', 'LineWidth', 2, 'color', br);
    ylabel('Cases')
    xline(ABflowdates(split),'k--','linewidth',3)
    ylabel('viral RNA copies')
    
    axis tight
    dim = [.15 .1 .3 .3];

    f1 = gcf;
    frac_S = 100/F;
    if split == 25
        data_used = 7;
    elseif split == 26
        data_used = 6;
    end

    namess = strcat('ADOO_',num2str(frac_S),'_',num2str(data_used));
    exportgraphics(f1,strcat(namess,'.png'),'Resolution',300)
    write(T_boot,'Bootstrap_ADOO.xlsx','Sheet',1)

%% functions
function err = obj_fun(param,sampleDates,data1,data2)
%Note, we only model 2nd strain.
global F
    N0 = 1.3*10^6/F;

    I10 = 0;
    I20 = data2(1)/param(5);
    R0 = 0;
    S0 = N0 - (I10 + I20 + R0);               
    ICs  = [S0 I10 I20 R0 data1(1) data2(1) I10+I20];

    tspan = sampleDates(1):1:sampleDates(end);
    [~,Y] = ode45(@SIIR,tspan,ICs,[],param);
    
    % shift back to get indices for desol
    ind = sampleDates - (sampleDates(1));
    ind = ind(2:end);

    % get daily virus for Delta
    Delta = Y(:,5);
    dailyDelta = diff(Delta);
    dailyDelta = dailyDelta(ind);

    % get daily virus for Omicron
    Omicron = Y(:,6);
    dailyOmicron = diff(Omicron);
    dailyOmicron = dailyOmicron(ind);

    temp1 = 0;
    adiff1 = rmmissing(temp1);
    
    temp2 = log10(data2(2:end)) - log10(abs(dailyOmicron));
    adiff2 = rmmissing(temp2);


    err = sum((adiff1).^2) + sum((adiff2).^2);

end

function dy = SIIR(t,y,param)
%Note, we only model 2nd strain.
    beta_1 = param(1);
    beta_2 = param(2);
    d1 = param(3);
    d2 = param(4);
    alpha = param(5);
   
    gamma1 = 1/8;
    gamma2 = 1/8;
    a = 0;%0.000027*2300000;   % MA birth rate per capita x population
    d = 0;%.00003477; % 1/78.79/365

    % per capita vaccinations (see vax_data.m)
    v = 0;%7.294594992636230e-04;

    dy = zeros(7,1);
    
    % not immune
    S = y(1);
    I1 = y(2);
    I2 = y(3);
    Rl = y(4);
    
    dy(1) = a - beta_1*S*I1 - beta_2*S*I2 - v*S - d*S;
    dy(2) = beta_1*S*I1 - gamma1*I1 - d1*I1;
    dy(3) = beta_2*S*I2 + beta_2*I2*Rl - gamma2*I2 - d2*I2;
    dy(4) = gamma1*I1- beta_2*I2*Rl - v*Rl - d*Rl;

    % track cumulative virus per strain and cases
    dy(5) = alpha*(I1);
    dy(6) = alpha*(I2);
    dy(7) = beta_1*S*I1 + beta_2*S*I2 + beta_2*Rl*I2;

end
