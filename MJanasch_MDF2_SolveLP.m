% edited by Luise Zeckey, KTH

function [solution]=MJanasch_MDF2_SolveLP(model,RT,n,m)%,beq)

% Define Objective
c = [
    zeros(n,1);
    1
    ];

% Define S_Ratio matrices

S_ratio_max = zeros(3,n);

% NADPH/NADP
S_ratio_max(1,10) = 1;
S_ratio_max(1,22) = -1;
% ATP/ADP
S_ratio_max(2,21) = 1;
S_ratio_max(2,4) = -1;
% NADH/NAD
S_ratio_max(3,27) = 1;
S_ratio_max(3,28) = -1;

S_ratio_min = -1*S_ratio_max;

% Define max and min vector for ratios

R_max = [
    27.7303;
    15.463;
    0.349414     % 0.349415
    ];

R_min = [
    0.0378572;
    0.0594318;
    0.00463872
    ];


% Define A
A = [
     transpose(model.rawS), ones(m,1);
     eye(n),                zeros(n,1);                     % eye is matlab function for identity matrix, diagonal 1
    -eye(n),                zeros(n,1);
    S_ratio_max,            zeros(3,1);
    S_ratio_min,            zeros(3,1)
    ];

% Define b
b = [
    -model.reactions.dG0/RT;
     log(model.metabolites.uconc);
    -log(model.metabolites.lconc);
    log(R_max);
    -log(R_min)
    ];

solution.c=c;
solution.A=A;
solution.b=b;


% Primal Problem: Finding MDF
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(transpose(-c),A,b);

solution.z=X;

end