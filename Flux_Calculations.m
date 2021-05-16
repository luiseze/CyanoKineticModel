%% Flux Simulations using Synechocystis Genome-Scale Metabolic Models
% Markus Janasch, Ph.D. Student KTH; edited by Luise Zeckey, KTH
% created:21/02-23, last modified: 21/05-16

%% Tidy up
clear all
close all
clc

%% Load COBRA
initCobraToolbox()

%% Load the Model
load('Models/SYN_GEM_iSyn731_CycleSyn.mat'); % Model from Sarkar, D., Mueller, T. J., Liu, D., Pakrasi, H. B., & Maranas, C. D. (2019). A diurnal flux balance model of Synechocystis sp. PCC 6803 metabolism. PLoS Computational Biology, 15(1), e1006692. https://doi.org/10.1371/journal.pcbi.1006692

%% Set Solver
changeCobraSolver('gurobi','all');


%% Set Carbon Input
iSyn731 = changeRxnBounds(iSyn731,'EX_H2CO3',-10,'l');     % Maximum bicarbonate uptake rate 10 mmol/h/gDW
iSyn731 = changeRxnBounds(iSyn731,'EX_CO2',0,'l');          % No CO2 uptake, only bicarbonate
iSyn731 = changeRxnBounds(iSyn731,'rxn01870',0,'b');        % Glucose uptake shut off


%% Smaller bounds
iSyn731.ub(iSyn731.ub==10000)=100;
iSyn731.lb(iSyn731.lb==-10000)=-100;
iSyn731.ub(iSyn731.ub==200)=100;
iSyn731.lb(iSyn731.lb==-200)=-100;

%% change RuBisCO reaction
iSyn731 = addReaction(iSyn731, 'RuBisCO', 'reactionFormula', 'cpd00871[ca] + 0.97 cpd00011[ca] + 0.97 cpd00001[ca] + 0.03 cpd00007[ca] --> 1.97 cpd00169[ca] + 0.03 cpd00727[ca] + 2 cpd00067[ca]');

% block the other two RuBisCO reactions
iSyn731 = changeRxnBounds(iSyn731,'rxn00018',0,'b'); % RuBisC
iSyn731 = changeRxnBounds(iSyn731,'rxn02251',0,'b'); % RuBisO

%% make (1) NAD[c] + (1) NADPH[c] --> (1) NADH[c] + (1) NADP[c] (rxn00083) reversible
iSyn731 = changeRxnBounds(iSyn731,'rxn00083',-100,'l');
iSyn731.rev(561) = 1;

%% change boundaries based on Gopalakrishnan, S., Pakrasi, H. B., & Maranas, C. D. (2018). Elucidation of photoautotrophic carbon flux topology in Synechocystis PCC 6803 using genome-scale carbon mapping models. Metabolic Engineering, 47, 190–199. https://doi.org/10.1016/j.ymben.2018.03.008
iSyn731 = changeRxnBounds(iSyn731,'rxn00973',0.156812142,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00973',0.219052132,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00786',-0.268701252,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00786',0.401999547,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00256',0.156812252,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00256',0.224162989,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01476',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01476',0.0256,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00459',1.090850802,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00459',1.276653938,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01334',-7.431407342,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01334',-6.936703027,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00549',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00549',0.402,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00782',17.35164049,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00782',56.08435918,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00333',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00333',1.249210127,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01115',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01115',0.0256,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00558',-0.2988858875,'u');
iSyn731 = changeRxnBounds(iSyn731,'rxn00558',-0.351186186,'l');

iSyn731 = changeRxnBounds(iSyn731,'rxn00198',1e-07,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00198',0.219052032,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00248',-0.732862671,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00248',9.691973171,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00154',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00154',0.339419926,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01100',17.35164045,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01100',18.28065603,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01106',-1.276653014,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01106',-0.226538021,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00251',1.13,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00251',1.315803116,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01111',9.785684711,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01111',10.34309409,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00173',-1.030143109,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00173',-0.751438435,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00148',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00148',0.185803216,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01116',-6.653579061,'u');
iSyn731 = changeRxnBounds(iSyn731,'rxn01116',-7.025185293,'l');

iSyn731 = changeRxnBounds(iSyn731,'rxn00777',3.132321862,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00777',3.318124988,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn01345',6.93670313,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn01345',29.39840994,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00747',7.019412262,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00747',7.391016865,'u');

iSyn731 = changeRxnBounds(iSyn731,'rxn00604',0.0000001,'l');
iSyn731 = changeRxnBounds(iSyn731,'rxn00604',0.0256,'u');

%% Set Objective 
iSyn731 = changeObjective(iSyn731,'Ec_biomass_SynAuto');

%% Simple Solve FBA
FBA_Solution = optimizeCbModel(iSyn731,'max');

%% Relevant Outputs:

CBB = {'rxn00018';  
    'rxn01100';     
    'rxn00782';     
    'rxn00747';     
    'rxn00786';     
    'rxn00785';     
    'rxn01200';     
    'rxn00549';     
    'rxn01334';     
    'rxn01345';     
    'rxn00777';     
    'rxn01116';     
    'rxn01111';     
    'rxn01106';     
    'rxn00459';     
    'rxn00148';     
    'rxn00154';     
    'rxnPhosphoketolase';   
    'rxnPhosphoketolase2';  
    'rxn00173';     
    'ATPm';         
    'H2ASE_syn';    
    'rxn00558';     
    'rxn00604';     
    'rxn01476';     
    'rxn01115';     
    'rxn00256';     
    'rxn00973';     
    'rxn00198';     
    'rxn00248';     
    'rxn00251';     
    'rxn00159';     
    'rxn00980';     
    'rxn00333';     
    'rxn00799';     
    'rxn01333';     
    'rxn02251';    
    'RuBisCO'
    };

CBB_Names = {'RuBisC';
    'PGK';
    'GAPD';
    'TPI';
    'ALD';
    'TKT1';
    'TKT2';
    'FBPase';
    'FBA';
    'SBPase';
    'RPI';
    'RPE';
    'PRK';
    'PGM';
    'ENO';
    'PYK';
    'PDH';
    'XFPK1';
    'XFPK2';
    'PTA';
    'ATPSyn';
    'NADPase';
    'PGI';
    'ZWF';
    'DEVB';
    'GND';
    'CS';
    'ACO';
    'ICD';
    'MDH';
    'PPK';
    'ME';
    'PGP';
    'GLCO';
    'FUMC';
    'TAL';
    'RuBisO';
    'RuBisCO'
    };    

CBB_Fluxes = FBA_Solution.v(KX_FindIndex(iSyn731.rxns,CBB)); 

KX_FindIndex(iSyn731.rxns,CBB)
%% Get nicer output format

for i = 1:length(CBB)
    Output.CBB(i).Reaction_ID = CBB_Names(i);
    Output.CBB(i).Flux = CBB_Fluxes(i);
end

%%

iSyn731.metCharges(1005) = NaN;
iSyn731.metCharges(1006) = NaN;
iSyn731.metCharges(1007) = NaN;
iSyn731.metCharges(1008) = NaN;
iSyn731.metCharges(1009) = NaN;

% set lower bound to 95% of the value of the objective function in FBA
iSyn731 = changeRxnBounds(iSyn731,'Ec_biomass_SynAuto',0.95*FBA_Solution.f,'l');

tic
solutions_Sampling = randomSampling(iSyn731, 1000, true, false, false, [], true);
toc

solutions_Sampling_T = transpose(solutions_Sampling);

solutions_Sampling_T_full = full(solutions_Sampling_T);

% Make a structure that we can export to Excel for easier evaluation
rsOut_Cyano(:, 1) = iSyn731.rxns; % The first column we want reaction IDs
rsOut_Cyano(:, 2) = num2cell(mean(solutions_Sampling_T_full)); % Second column we want the means of the fluxes. As this structure is no longer an all-numerical matrix, we need to convert from numbers to cells (num2cell).
rsOut_Cyano(:, 3) = num2cell(std(solutions_Sampling_T_full)); % Standard deviation
rsOut_Cyano(:, 4) = num2cell(min(solutions_Sampling_T_full)); % Minimum flux
rsOut_Cyano(:, 5) = num2cell(max(solutions_Sampling_T_full)); % Maximum flux

for i = 1:length(CBB)
    rsOut_CBB(i, 1) = iSyn731.rxns(KX_FindIndex(iSyn731.rxns,CBB{i}));
    rsOut_CBB(i, 2) = CBB_Names(i);
    rsOut_CBB(i, 3) = rsOut_Cyano(KX_FindIndex(rsOut_Cyano(:, 1), CBB{i}), 2);
    rsOut_CBB(i, 4) = rsOut_Cyano(KX_FindIndex(rsOut_Cyano(:, 1), CBB{i}), 3);
    rsOut_CBB(i, 5) = rsOut_Cyano(KX_FindIndex(rsOut_Cyano(:, 1), CBB{i}), 4);
    rsOut_CBB(i, 6) = rsOut_Cyano(KX_FindIndex(rsOut_Cyano(:, 1), CBB{i}), 5);
end
