% adapted from Janasch, M., Asplund-Samuelsson, J., Steuer, R., & Hudson, E. P. (2019). Kinetic modeling of the Calvin cycle identifies flux control and stable metabolomes in Synechocystis carbon fixation. Journal of Experimental Botany, 70(3), 1017–1031. https://doi.org/10.1093/jxb/ery382

% Extract parameters from a set of CBB SKM mat files

% Create list of infiles; assuming only one set of sampled parameters
% is present in the input directory
indir = 'parameter_sampling_results/';
outdir = 'parameter_sampling_results/par';
infiles = dir(fullfile(indir,'*.mat'));
infiles = {infiles.name}.';
N = length(infiles);
tic
% Iterate over the infiles
for n = 1:N
  fprintf(2, '%3.1f%%\r', n/N*100)
  % Load data
  load(char(fullfile(indir, infiles(n))));
  % Acquire metabolite concentration set ID
  infile_split = strsplit(char(infiles(n)), '_');
  infile_n = str2num(char(infile_split(3)));
  % Extract parameters and vector indicating stability or not
  parameters = horzcat(DataOut.StabilityIndicator, DataOut.Parameters);
  % Write to tab-delimited file
  outfile = fullfile(outdir, strcat(num2str(infile_n), '.tab'));
  dlmwrite(outfile, parameters,'delimiter', '\t');
end
toc
% Save a header
outfile = fullfile(outdir, 'par_header.long.txt');
dlmwrite(outfile, char(DataOut.ParID), 'delimiter', '')

%fprintf(2, '%3.1f%%\n', n/N*100)
