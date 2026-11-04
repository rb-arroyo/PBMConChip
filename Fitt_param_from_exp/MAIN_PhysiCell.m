clear all; close all; clc;

%% Reading and analysis of the experimental data
tic;
cond = 4; % Experimental condition: 2 or 4 mg/mL
D = 0; % Day 0 or day 1
B = 1; % Bootstrap number (how many resampling-with-replacement samples are generated)
fname = sprintf('./EXP_RESULTS/trajectories_collagen%01d_day%01d/',cond,D);
Filedir = fname; % experimental data directory
Filename = "all_trajectories.mat"; % file name
load(fullfile(Filedir,Filename)) % load data

%% Processing of the experimental data

for bst = 1:B
experimental_data_processor(M,cond,D,bst); % Processes the experimental data

%%
savefile = sprintf('experimental_data_B%03d.mat',bst);
load(savefile); % Load the processed data

% Time-dependent weights
lambda_mx_t = Ncells_t ./ Sxx;
lambda_my_t = Ncells_t ./ Syy;

lambda_sx_t = sqrt((Ncells_t-1) ./ (2*Sxx.^2));
lambda_sy_t = sqrt((Ncells_t-1) ./ (2*Syy.^2));

pmx = sqrt(lambda_mx_t);
pmy = sqrt(lambda_my_t);

psx = lambda_sx_t;
psy = lambda_sy_t;

QoI_real = [pmx.*Xmean;pmy.*Ymean;psx.*Sxx;psy.*Syy]; % Quantity of interest from the experiments

EDATA.QoI = QoI_real;
EDATA.LSxx = length(Sxx);
EDATA.pmx = pmx;
EDATA.pmy = pmy;
EDATA.psx = psx;
EDATA.psy = psy;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
Nreplicas = 5; % number of synthetic data replicates
nsampling = 1; % how often the information is sampled

% Initial value of the variables to optimize
b_ini = 0.15; 
s_ini = 4.0;
T_ini = 3.0;

% Random initial value within a range
% b_min = 0.0;
% b_max = 0.5;
% c_min = 1.0;
% c_max = 5.0;
% Tpers_min = 2.0;
% Tpers_max = 5.0;
% 
% b0 = b_min + (b_max - b_min)*rand();
% c0 = c_min + (c_max - c_min)*rand();
% Tpers0 = Tpers_min + (Tpers_max - Tpers_min)*rand();

bmax = 1.0;
smax = 7.0;
Tmax = 10.0;
lb = [0.0;1.0;2.0]; % Optimize 3 parameters
ub = [bmax;smax;Tmax]; % Optimize 3 parameters

% 3-parameter case (initialization and loss function)    
%x_ini = [b0,c0,Tpers0];
x_ini = [b_ini;s_ini;T_ini];
myFun = @(x) loss_function_PhysiCell(x,EDATA,nsampling,Nreplicas);

% Optimization options
step_TOL = 1e-4;
fun_TOL = 1e-4;
opt_TOL = 1e-4;
step_size = 0.01;
options = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt', ...
     'StepTolerance',step_TOL, ...
     'FunctionTolerance',fun_TOL, ...
     'FiniteDifferenceStepSize',step_size, ...
     'OptimalityTolerance',opt_TOL);

    
%% OPTIMIZER

[x_opt,resnorm,residual,exitflag,output] = lsqnonlin(myFun,x_ini,lb,ub,options);

b_opt = x_opt(1);
s_opt = x_opt(2);

resultfile = sprintf('results_c%01d_d%01d_B%03d.mat',cond,D,bst);
save(resultfile)

end

toc;
