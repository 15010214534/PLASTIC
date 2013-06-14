% SETUP for PLASTIC: This file contains
% major simulation variables that specify 
% the network.
% Guy Billings 17/02/2007

% Load a previous data file?

%load('');

% If you loaded a file you need to set this switch to '1'
% This ensures that weight matricies are not reinitialised
% Parameters will STILL BE OVERWRITTEN

loadfromold=0;

%% NETWORK SETUP:
%----------------

% Network size: number of integrate and fire neurons and inputs

N=1; 		
inputs=800; 

% Set the fraction of feedforward connections be used. If less than
% 100% then connections are randomly assigned

input_frac=1;   % Fraction of feed forward synapses to connect (1 for all connected)

% Set up the recurrent connections (if more than one IF neuron)

in=0; 		% all to all recurrent inhibition 
ex=0;		% recurrent excitation with some radius
wdrec=0;	% inhibitory recurrent connection strength (nS) 
radius=0; 	% number of neighbours for antagonist connections
an_str=0;	% excitatory antagonists strength (nS)
mult=1;     % excitatory recurrent multiplier

% Set up weight seeding

seedme=0;	% Seed feedforward weights for map formation
width=0;    % width of seeding	
rfac=0;     % seeding strength factor (multiplies original values)

% Switch determining the weight dependence of the Spike timing
% dependent plasticity learning rule

wdep=1;		% Set to 0 for weight independent rule

%% INPUT SETUP:
%----------------

% Base rate and variance in that rate for protocol 0:

zero_base=10;
zero_var=0.3;

% Peak and base rates for spatially periodic input stimulus
% Only important for protocols 1,2,3 and 4

input_peak=80;  
input_base=10;  

% Foreground and background rates and variances for protocol 5

fg_rate=90;	
bg_rate=10;     
fg_var=0.3;
bg_var=0.3;

% Specify the stimulus protocol

protocol=0	; %(0->5)

% 0: correlated or uncorrelated poisson inputs
% 1: Rate modulated periodic stimulus (randomly moving)
% 2: Receptive field test mode 
% 3: Program stimulus (see later)
% 4: Rate modulated periodic stimulus (stationary)
% 5: Mixed frequency Poisson inputs

% Setup for protocol 4 and 3
stim_loc=400;		% Location of stationary stimulus
stim_dur=0.02;		% Correlation timescale

% Programable protocol, using the spatially periodic input stimulus:
% -ve value for input provides an interval with no input
% +ve value is time in seconds
% (Do not set interval times less than 1ms)
% Progam format: 
% cycle=[stimulus location,stimulus duration(s);repeat...];
% repeats=number of times to repeat the specified cycle

cycle=[100,0.01;-1,0.02;200,0.01]'; 
repeats=100;

% Set input correlation: corrf is an array that specifies two groups of synaptic inputs, one
% uncorrelated and one correlated (rate correlations). Set the element of corrf to '1' in order
% assign that input to the correlated group. corr_com sets the fraction of the base rate 
% that correlated inputs have in common.

corrf(1:inputs)=zeros(1,inputs);  
corr_com=0.3;

%% SIMULATION SETUP
%------------------

% Test mode to measure tuning curve of loaded recurrent network

TEST=0;  

% Track: Set to '1' to record a snapshot of all weights after every
% 'snap_int' outerloop intervals (i.e. after every correlation time interval
% of the inputs: mean_int)

TRACK=0; 
snap_int=100; 

% Loop setup

mean_int=0.02; 	% Correlation time interval (s)
sim_steps=100;  	% Total repetitions of correlation time interval

% Simulation segmentation: For long simualtions you may want PLASTIC to 
% write to a .mat file every so often. This safegaurds data in the event of
% a crash or a computer failure

sim_repeats=1;  % number of repeats 

% Set up the simulation output:

save_me=1; 				% Set to 1 to save results to a file
save_ini=1;				% Save the initial state of the network
filen='plastic';		% Specify a string for the filename
extension='.mat';		% Specify file extension
