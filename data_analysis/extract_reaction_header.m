% adapted from Janasch, M., Asplund-Samuelsson, J., Steuer, R., & Hudson, E. P. (2019). Kinetic modeling of the Calvin cycle identifies flux control and stable metabolomes in Synechocystis carbon fixation. Journal of Experimental Botany, 70(3), 1017–1031. https://doi.org/10.1093/jxb/ery382

% Extract the reactions header from a CBB SKM model file
model_file = 'Expanded_Model.mat';
outdir = 'parameter_sampling_results/';

load(model_file);
header = char({N.reaction.id}.');
outfile = fullfile(outdir, 'cbb_reaction_header.long.txt');
dlmwrite(outfile, header, 'delimiter', '');
