function [A] = load_data_PhysiCell(dir,nsteps)
% Reading the simulated data
dt = 2; % Data saved every dt min
%! Automate this. Perhaps read it from the .xml?
t = nsteps-1; % Number of time steps

for i = 1:t+1
    filename = sprintf('output00000%03d_cells.mat',i-1);
    T{i} = load ([dir '/output/' filename]);
    ti(i) = (i-1)*dt;
    M{i} = T{i}.cells(1:4,:);
    [~,C] = size(M{i});
    for j = 1:C
        ind = M{i}(1,j);
        A{ind+1}(i,1) = ti(i);
        A{ind+1}(i,2) = M{i}(2,j);
        A{ind+1}(i,3) = M{i}(3,j);
        A{ind+1}(i,4) = M{i}(4,j);
    end
end
I = 1;
for i = 1:length(A) 
    if ~isempty(A{i})
        Amod{I} = A{i};
        I = I+1;
    end
end
A = Amod;
for i = 1:length(A)
    A{i}(all(A{i} == 0,2),:) = [];
end

nameFile = 'PHYSICELL_DATA';
save(nameFile,'A','t','dt');
end
