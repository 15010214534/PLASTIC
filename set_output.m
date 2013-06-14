% OUTPUTS for PLASTIC: This file specifies
% what data PLASTIC should output. The data
% appears in the MATLAB workspace.
% Guy Billings 17/02/2007

% Output flags set to one for output
% of that data in to workspace

output(1)=0; % Pre-synaptic spike reconstruction data
output(2)=0; % Voltage output: 0 for none, or the neuron index 
output(3)=0; % Postsynaptic spike reconstruction data
output(4)=0; % Stimulus location and duration
output(5)=0; % overlap of network with initial state
output(6)=0; % average weights for the neurons
output(7)=1; % frequency in sampling period
output(8)=0; % average excitatory reccurents
output(9)=0; % store weight variances along with averages
output(11)=0; % Euclidean distance
output(12)=0; % Pearson's R

output(13)=0; 			% Weight tracking (number of wieghts; 0 if none)
output(14)=0; 			% ==1 for reduced tracking at each sample point rather than every timestep
output(15)=0; 			% current input recording flag (*net* recurrent current and input)
output(16)=snap_int; 	% Wieght snapshot interval
nflags=ones(1,N); 		% Flags for neurons to be included in the average

% Setup how often the variables should be sampled:

stime=0.5; %(s)	
sample_period=round(stime/const_master);
output(10)=sample_period; % sampling period for auxillary variables
                          % (timesteps)

