function Data_ERP_N3()
% This function creates the model data depicted in Figure 7 of
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

% To ensure availability of the simulation routine and the utility functions 
% from fieldtrip it should be called from Create_Data.m

% Simulation parameters
load('Data/Parameter_N3');   
            
% Simulation duration
T = 3600;         % duration of the simulation

% stimulation parameters
% first number is the mode of stimulation
% 0 == none
% 1 == semi-periodic
% 2 == phase dependend
    
var_stim    = [ 2;          % mode of stimulation
                70;         % strength of the stimulus              in Hz (spikes per second)
                80;       	% duration of the stimulus              in ms
                5;          % time between stimulation events       in s  (ISI)
                0;          % range of ISI                          in s  [ISI-range,ISI+range]  
                2;          % Number of stimuli per event
                1050;        % time between stimuli within a event   in ms         
                450];       % time until stimuli after minimum      in ms
 

% Model Output is given as: 
% 1. Pyramidal membrane voltage         in mV
% 2. Thalamic relay membrane voltage    in mV
% 3. Thalamic calcium concentration     in mu mol
% 4. Thalamic I_h activation            unitless
% 5. Stimulation markers                in sampling rate
[Vp, Vt, Ca, ah, Marker_Stim] = TC_mex(T, Param_Cortex, Param_Thalamus, Connectivity, var_stim);

% Power of spindle bands (BP filtered)
Fs          = length(Vp)/T;
Vp_FSP      = ft_preproc_hilbert(ft_preproc_bandpassfilter(Vp, Fs, [12, 15], 513, 'fir'), 'abs').^2;
Vp_SSP      = ft_preproc_hilbert(ft_preproc_bandpassfilter(Vp, Fs, [09, 12], 513, 'fir'), 'abs').^2;

% Check whether stimuli are too early
Range_ERP   = [-1, 3];
    
Marker_Stim = Marker_Stim(Marker_Stim>-Range_ERP(1)*Fs);
Marker_Stim = Marker_Stim(Marker_Stim< (T-Range_ERP(2))*Fs);

% Define the matrices
N_Stim      = length(Marker_Stim);
time_event  = linspace(Range_ERP(1), Range_ERP(2), (Range_ERP(2)-Range_ERP(1))*Fs+1);
Events      = zeros(length(time_event), N_Stim);
Events_T    = zeros(length(time_event), N_Stim);
Events_FSP  = zeros(length(time_event), N_Stim);
Events_SSP  = zeros(length(time_event), N_Stim);

for i=1:N_Stim
    Events(:,i)     =  Vp((Marker_Stim(i)+Range_ERP(1)*Fs)+1:(Marker_Stim(i)+Range_ERP(2)*Fs+1));
    Events_T(:,i)   =  Vt((Marker_Stim(i)+Range_ERP(1)*Fs)+1:(Marker_Stim(i)+Range_ERP(2)*Fs+1));
    Events_FSP(:,i) =  Vp_FSP((Marker_Stim(i)+Range_ERP(1)*Fs)+1:(Marker_Stim(i)+Range_ERP(2)*Fs+1));
    Events_SSP(:,i) =  Vp_SSP((Marker_Stim(i)+Range_ERP(1)*Fs)+1:(Marker_Stim(i)+Range_ERP(2)*Fs+1));
end

save('Data/ERP_Stim_Model', 'Events', 'Events_T', 'Events_FSP', 'Events_SSP', 'Marker_Stim');
save('Data/Ca_Depletion', 'Vp', 'Vt', 'Ca', 'ah', 'Marker_Stim');
end