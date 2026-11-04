function QoI = quantity_of_interest_PhysiCell(DATA,nsampling)

[mean,cova] = data_mean_processor_PhysiCell(DATA,nsampling);
QoI = [mean(:);cova(:)];

end
