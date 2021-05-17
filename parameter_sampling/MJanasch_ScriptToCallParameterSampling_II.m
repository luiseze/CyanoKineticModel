% This script calls the parameter sampling and coefficient calculating 
% functions for the kinetic models
% Markus Janasch, Ph.D. Student, KTH; edited by Luise Zeckey, KTH
% Created: 2017-03-22, last modified: 2021-05-16

Iterations = 1000;
InputDataStructure = 'Expanded_Model.mat';
MetConcSamplingData = 'Metabolite_Pool_Concentrations.txt';
MetConcData_RAW = importdata(MetConcSamplingData);

tic
for i = 1:length(MetConcData_RAW.data)
    NrOfMetDataSet = i;
    [DataOut, N] = MJanasch_Parameter_Sampling_II(Iterations,InputDataStructure,MetConcData_RAW.data(NrOfMetDataSet,:),MetConcData_RAW.textdata);
    set = num2str(i);
    filename = strcat('skm_sample_', set, '_2021-05-06.mat');
    filedir = strcat('parameter_sampling_results/', filename);
    DataSetOut = filedir;
    save(DataSetOut,'DataOut');
end
toc
