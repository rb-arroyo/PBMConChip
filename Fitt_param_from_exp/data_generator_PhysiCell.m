function data_generator_PhysiCell (param,Nreal)

% Generating synthetic data
% Run PhysiCell simulations and copy the results
command1 = '../PhysiCell-master/BashScripting/OPTphysicell.sh "'; % Runs the bash script
command2 =  command1 + string(param.b) + '" ' + '"' + string(param.s) + '" ' + '"' + string(param.T) + '" ' + '"' + string(Nreal) + '"';
status = unix(command2);

pause(1);

disp('Synthetic PhysiCell data generated');

end
