function Error_mask = loss_function_PhysiCell(x,DATA,nsampling,nrealizations)

% Experimental data
QoI_real = DATA.QoI;
LSxx = DATA.LSxx;
pmx = DATA.pmx;
pmy = DATA.pmy;
psx = DATA.psx;
psy = DATA.psy;

% Model parameters
param.b = x(1);
param.s = x(2);
param.T = x(3);

% Generate the data with the corresponding number of replicates
data_generator_PhysiCell(param,nrealizations);

% Analyze the data from the replicates
% Reading the simulated data
dir = '../PhysiCell-master/Results/OPT';
QoI = [];
for i=1:nrealizations
    MP = load_data_PhysiCell([dir int2str(i)],LSxx);
    OPTDATA.M = MP;
    OPTDATA.nsteps(1) = length(MP{1});
    OPTDATA.dt = MP{1}(2,1)-MP{1}(1,1);
    QoI = [QoI,quantity_of_interest_PhysiCell(OPTDATA,nsampling)];  %% QoI for cells in PhysiCell
end
QoI_sim = mean(QoI,2);
QoI_sim = [pmx.*QoI_sim(1:LSxx);pmy.*QoI_sim(LSxx+1:2*LSxx); psx.*QoI_sim(2*LSxx+1:3*LSxx);psy.*QoI_sim(3*LSxx+1:end)];
OPTDATA.QoI_sim = QoI_sim;

save PHY_DATA OPTDATA

Error = QoI_sim - QoI_real;
mask = ~isnan(Error);
Error_mask = Error(mask);

save('OPT.mat', 'Error', 'Error_mask')

end
