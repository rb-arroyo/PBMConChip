function experimental_data_processor(M,cond,D,B)

Mexp = M';

% Initial and final positions of the cells
nagents = length(Mexp); % total number of cells identified in the videos

tini = zeros(nagents,1);       % initial times of each cell
Xini = zeros(nagents,1);       % x position at t=0
Yini = zeros(nagents,1);       % y position at t=0
Xfinal = zeros(nagents,1);     % final x position
Yfinal = zeros(nagents,1);     % final y position

%% Analyze the initial and final positions of the raw experimental data
%figure(1); hold on
for i=1:nagents
    tini(i) = Mexp{i}(1,1);
    Xini(i) = Mexp{i}(1,2);
    Yini(i) = Mexp{i}(1,3);
    Xfinal(i) = Mexp{i}(end,2);
    Yfinal(i) = Mexp{i}(end,3);
end 
%plot(Xini,Yini,'bo'); hold on
%plot(Xfinal,Yfinal,'ro')
%legend('Initial','Final')

%% Generate bootstrap sample (resampling with replacement)

if (B > 0)
    % Select random indices with replacement (from 1 to N)
    idx_resample = randi([1, nagents], 1, nagents);

    % 2. Preallocate the new struct for memory efficiency
    agentes_boot = cell(1,nagents);

    % 3. Build the new dataset
    for i = 1:nagents
        % Extract the sampled index
        id_original = idx_resample(i);

        % Copy the exact [t, x, y] matrix into the new dataset
        agentes_boot{1,i} = Mexp{1,id_original};
    end
Mexp = agentes_boot;

%% Analyze the initial and final positions of the bootstrap data
%figure(3); hold on
 for i=1:nagents
     tini(i) = agentes_boot{i}(1,1);
     Xini(i) = agentes_boot{i}(1,2);
     Yini(i) = agentes_boot{i}(1,3);
     Xfinal(i) = agentes_boot{i}(end,2);
     Yfinal(i) = agentes_boot{i}(end,3);
 end 
%plot(Xini,Yini,'bo'); hold on
%plot(Xfinal,Yfinal,'ro')
%legend('Initial','Final')

end
%% Experimental time data
all_times = []; % Times for which data are available in the experiment

for i = 1:nagents
    all_times = [all_times; Mexp{i}(:,1)];
end

% Sorted unique times (in case some dt have no cells)
t_unico = unique(all_times);

%% Data processing to compute statistics

% Parameter initialization
nsteps = length(t_unico); % total number of time steps
Ncells_t = zeros(nsteps,1); % number of cells at each t

Xmean = NaN(nsteps,1);
Ymean = NaN(nsteps,1);
Sxx = NaN(nsteps,1);
Syy = NaN(nsteps,1);
Sxy = NaN(nsteps,1);

% Metrics computation (mean and variance)

tol = 1e-6;  % time tolerance

for k = 1:nsteps
    t = t_unico(k);
    X = []; Y = [];
    for i = 1:nagents
        traj = Mexp{i};
        idx = find(abs(traj(:,1) - t) < tol, 1);
        if ~isempty(idx)
            X(end+1) = traj(idx,2);
            Y(end+1) = traj(idx,3);
            x{i} = X;
            y{i} = Y;
        end
    end
    Ncells_t(k) = length(X);

    if Ncells_t(k) >= 1
        Xmean(k) = mean(X);
        Ymean(k) = mean(Y);
        
        if Ncells_t(k) >= 2
            Sj = cov(X,Y);
            Sxx(k) = Sj(1,1);
            Syy(k) = Sj(2,2);
            Sxy(k) = Sj(1,2);
        end
    end
end

Xpos = Xini;
Ypos = Yini;

%% Save the data
% Matrix with the initial times and positions and the cell IDs
% CELLS = [tini,Xpos,Ypos,zeros(N,1),zeros(N,1)];
CELLS = [tini,Xpos,Ypos,zeros(nagents,1)];
C_ordenada = sortrows(CELLS,1);

% Check whether there are cells at time 0
CSVXini = C_ordenada(C_ordenada(:,1) == 0,2);
CSVYini = C_ordenada(C_ordenada(:,1) == 0,3);

[nagents_ini,~] = size(CSVXini);

if (nagents_ini > 0)
    CSVCELLS = [CSVXini,CSVYini,zeros(nagents_ini,1),zeros(nagents_ini,1)];
else
    CSVCELLS = [-100.0, -200.0, 0.0, 0];
end

vals = cellfun(@(M) M(1,1), Mexp);   % extracts the first row, first column
[~, idx] = sort(vals);               % sorted indices
Mexp_ordenada = Mexp(idx);  % reorders the cell array

%% Save the initial positions to pass them to PhysiCell
filename1 = sprintf('../PhysiCell-master/config/initial_cells_opt_c%01d_d%01d_B%03d.csv',cond,D,B);
filename2 = sprintf('../PhysiCell-master/config/cells_opt_c%01d_d%01d_B%03d.csv',cond,D,B);
filename3 = sprintf('../PhysiCell-master/config/initial_cells_opt_c%01d_d%01d.csv',cond,D);
filename4 = sprintf('../PhysiCell-master/config/cells_opt_c%01d_d%01d.csv',cond,D);
savefile = sprintf('experimental_data_B%03d.mat',B);
%writematrix(CSVCELLS,filename1);
%writematrix(C_ordenada,filename2);
writematrix(CSVCELLS,filename3);
writematrix(C_ordenada,filename4);

save (savefile);

end
