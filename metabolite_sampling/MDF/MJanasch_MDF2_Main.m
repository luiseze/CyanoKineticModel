%% Max-Min Driving Force 2.0
% Based on "Pathway Thermodynamics Highlights Kinetic Obstacles in Central Metabolism" by Noor et al., 2014
% Markus Janasch, Ph.D. Student
% Microbial Metabolic Engineering Group
% KTH School of Biotechnology, Science For Life Laboratory, Stockholm
% edited by Luise Zeckey, KTH
% Created: 2019/09/06, Last modified: 2021/05/16


%% Add Support functions
addpath('Support')

%% Tidy up
clear all
close all
clc

%% Define Parameters
T = 303.15;                     % Temperature [K]
R = 0.008314;                   % Universal gas constant [kJ/(mol*K)]
RT=R*T;


%% Load in Data
MDF2ModelRaw = MJanasch_MDF2_LoadData('MDF_Data.xlsx'); % Loads dG0, concentration boundaries and creates S

model = MDF2ModelRaw;

% Get number of reactions m and number of metabolites n
[n,m] = size(model.rawS);

% Get standard delta G's (transpose to get from column vector to row vector
dG0trans = transpose(model.reactions.dG0(:,1));

%% Set Solver
changeCobraSolver('gurobi','LP'); % Employ COBRA's function to easily change LP solvers (I can choose between matlab, gurobi and mosek

%% Employ MDF Algorithm

solutionLP = MJanasch_MDF2_SolveLP(model,RT,n,m);%,beq);

A = full(solutionLP.A);
b = solutionLP.b;
c = solutionLP.c;

dG = model.reactions.dG0+RT*transpose(model.rawS)*solutionLP.z(1:n);
dGtrans = transpose(dG);

conc = exp(solutionLP.z(1:n));
MDF = solutionLP.z(end,1)*RT;


%% Plot
% %fprintf('\n MDF [kJ/mol]: %d \n',round(MDF,4)) % Output
% 
disp(['Max-min Driving Force [kJ/mol] : ',num2str(round(MDF,4))])
% 
LeveldG0    = 0;
%LeveldG     = zeros(r-1,1);
LeveldG    = 0;
%for l = 1:r-1;
    for i = 1:length(dG0trans)
        LeveldG0(i+1) = LeveldG0(i) + dG0trans(i);
        %LeveldG(l,i+1) = LeveldG(l,i) + dGtrans(i,l);
        LeveldG(i+1) = LeveldG(i) + dGtrans(i);
    end
%end

min = min(LeveldG);
max = max(LeveldG0);

figure (2)
plot(0:m,LeveldG0,'-o',0:m,LeveldG,'-x')
hold on;
%plot(0:m,LeveldG(1,:),'-x','linewidth',3)
title('Thermodynamic landscape - Central Carbon Metabolism', 'FontSize',15,'FontName', 'Arial')
xlabel('Reactions', 'FontSize',15,'FontName', 'Arial')
ylabel('\Delta_rG''-Level', 'FontSize',15,'FontName', 'Arial')
%legend('\Delta_rG''^{0}','optimized \Delta_rG''', 'FontSize',15,'FontName', 'Arial')
xlim([0 length(LeveldG0)-1])
%ylim([-163 max+1])
grid on;