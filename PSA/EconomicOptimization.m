% Script file to run Process Optimisation
clc; clear; close all;
format long
% addpath('GA_files\');
% addpath('CycleSteps\');

load('Params')
N = 10 ;
type = 'EconomicEvaluation' ;

% IsothermParams = [16 x 13] matrix (16 adsorbents) 
% (1)  Co-MOF-74     (2)  Cu-BTTRi    (3)  Cu-TDPAT    (4)  Mg-MOF-74
% (5)  MOF-177       (6)  Ni-MOF-74   (7)  NTU-105     (8)  Sc2BDC3 
% (9)  SIFSIX-2-Cu-i (10) SIFSIX-3-Ni (11) Ti-MIL-91   (12) UTSA-16
% (13) UiO-66(OH)2   (14) ZIF-8       (15) Zn-MOF-74   (16) Zeolite-13X
for i = 16 % Zeolite-13X
    
    % load parameters
    IsothermParams     = IsothermPar(i, :) ;
    material_propertry = SimParam(i, :)    ;
    
    material    = {}                 ;
    material{1} = material_propertry ;
    material{2} = IsothermParams     ;
    
    Function = @(x) PSACycleSimulation( x, material, type, N ) ; % Function to simulate the PSA cycle
    
    % initial variables
    [~, vars] = sortt(loadpopfile('Zeolite-13X_Process.txt'));
    % vars = [vars, ones(length(vars), 1), 1e4*ones(length(vars), 1), ones(length(vars), 1)];
    
    options            = nsgaopt();                          % create default options structure
    options.popsize    = 60;                                 % population size
    options.outputfile = 'Zeolite-13X_Economic.txt';
    options.maxGen     = 120;                                 % max generation
    
    options.vartype    = [1, 1, 1, 1, 1, 1] ;
    
    options.initfun={@Pop_Override, vars}   ;                % Supply variables from previous results
    
    options.numObj  = 2 ;                                    % number of objectives
    options.numVar  = 6 ;                                    % number of design variables
    options.numCons = 3 ;                                    % number of constraints
    
    %*********************** parameters in x *********************************
    % P_0         = x(1)       ;   % Adsorption pressure [Pa]
    % t_ads       = x(2)       ;   % Time of adsorption step [s]
    % alpha       = x(3)       ;   % Light product reflux ratio [-]
    % v_feed      = x(4)       ;   % Feed inlet velocity [m/s]
    % beta        = x(5)       ;   % Heavy product reflux ratio [-]
    % P_l         = x(6)       ;   % Purge Pressure [Pa]
    options.lb = [1e5,  10, 0.01, 0.1, 0, 1e4]   ;           % lower bound of x
    options.ub = [10e5, 1000, 0.99, 2, 1, 5e4]   ;           % upper bound of x
    %*************************************************************************
    options.nameObj = {'-productivity','energy'} ;           % the objective names are showed in GUI window.
    options.objfun  = Function                   ;           % objective function handle
    
    options.useParallel = 'yes' ;                            % parallel computation is non-essential here
    options.poolsize     = 32   ;                            % number of worker processes
    
    result = nsga2(options)     ;                            % begin the optimization!
    
    
    % re-optimize using 30 finite volumes
    N = 30 ;
    Function = @(x) PSACycleSimulation( x, material, type, N ) ; % Function to simulate the PSA cycle
    
    options.objfun     = Function           ;                % objective function handle
    options.initfun    = {@initpop, result} ;                % Supply variables from previous results
    options.maxGen     = 80                ;                 % populaion size
    options.outputfile = 'Zeolite-13X_Economic_2.txt'         ;
    result2            = nsga2(options)     ;                % begin the optimization!
    
end

