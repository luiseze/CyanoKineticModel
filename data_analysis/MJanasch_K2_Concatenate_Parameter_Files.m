% Markus Janasch, Ph.D. Student, KTH

% Concatenate Parameter files

% Create list of infiles; assuming only one set of sampled parameters
% is present in the input directory
indir = 'parameter_sampling_results/par';
outdir = 'parameter_sampling_results/';
infiles = dir(fullfile(indir,'*.tab'));
infiles = {infiles.name}.';
N = length(infiles);
Data_Total=[];
% Iterate over the infiles
r1 = 1;
r2 = r1+999; % 1000 Iterations
tic
for n = 1:N
  fprintf(2, '%3.1f%%\r', n/N*100)
  % Load data
  Data = load(char(fullfile(indir, infiles(n))));
  % Acquire metabolite concentration set ID
  infile_split = strsplit(char(infiles(n)), '.');
  
  infile_n = str2double(char(infile_split{1}));
  Data_Temp(:,2:246) = Data(:,:); % 245 parameters
  Data_Temp(:,1) = ones(1000,1)*infile_n; 

  
  Data_Total(r1:r2,:) = Data_Temp(:,:);
  r1 = r1+1000;
  r2 = r2+1000;
  
end
toc
outfile = 'concset_stability_parameters_noheaders.tab';
dlmwrite(fullfile(outdir, outfile), Data_Total, 'delimiter', '\t');