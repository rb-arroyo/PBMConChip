function [data_mean,data_cov] = data_mean_processor_PhysiCell(DATA,nsampling)

M = DATA.M;

% Data processing
N = length(M); % Number of agents

%% Times for which data are available
all_times = []; % Times for which data are available

for i= 1:N
    all_times = [all_times; M{i}(:,1)];
end

% Sorted unique times (in case some dt have no cells)
t_unico = unique(all_times);
nsteps = length(t_unico); % total number of time steps
%% Statistics computation
tol = 1e-6;  % time tolerance

for i = 1:nsteps
    t = t_unico(i);
    Pos_x = [];
    Pos_y = [];
    if rem(i,nsampling) == 0
        for j = 1:N 
            idx = find(abs(M{j}(:,1) - t) < tol, 1);
            
            if ~isempty(idx)
                pos_x = M{j}(idx,2);
                pos_y = M{j}(idx,3);
                Pos_x = [Pos_x;pos_x];
                Pos_y = [Pos_y;pos_y];
            end
        end 
       % I = I + 1;
        Ncells_t(i) = length(Pos_x);

        if Ncells_t(i) >= 1
        MEANX(i) = mean(Pos_x);
        MEANY(i) = mean(Pos_y);
        
        if Ncells_t(i) >= 2
            COV = cov(Pos_x,Pos_y);
            Sxx(i) = COV(1,1);
            Syy(i) = COV(2,2);
            Sxy(i) = COV(1,2);
        end
    end

        data_mean = [MEANX,MEANY];
        data_cov = [Sxx,Syy];

    end
  
end

end
