function Data_SO_Average(type)
% This function creates the model data depicted in Figure 6 of
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

% To ensure availability of the simulation routine and the utility functions 
% from fieldtrip it should be called from Create_Data.m
if nargin==0    
    type = 2;
end

switch type
    case 1;
        load_model = 'Data/Time_Series_N2.mat';
        load_data  = 'Data/KC_Average_data.mat';
        savename   = 'Data/KC_Average.mat';
    case 2;
        load_model = 'Data/Time_Series_N3.mat';
        load_data  = 'Data/SO_Average_data.mat';
        savename   = 'Data/SO_Average.mat';
end

% Load data
load(load_model, 'Vp');
load(load_data);

% Rename data
mean_ERP_data = mean_ERP_sham;
mean_FSP_data = mean_FSP_sham;
sem_ERP_data  = sem_ERP_sham;
sem_FSP_data  = sem_FSP_sham;

% Process model Data
T       = 3600;
Fs      = length(Vp)/T;

% Power of fast spindle band (BP filtered) and slow wave activity for peak
% detection
Vp_low      = ft_preproc_bandpassfilter(Vp, Fs, [0.25,4], 513,  'fir') + mean(Vp);
Vp_FSP      = ft_preproc_hilbert(ft_preproc_bandpassfilter(Vp, Fs, [12, 15], 513, 'fir'), 'abs').^2;

% Search for peaks
[~, x_SO]   = findpeaks(-Vp_low, 'MINPEAKHEIGHT', 68, 'MINPEAKDISTANCE', 0.2*Fs);

% Remove those events, that are too close to begin/end
x_SO        = x_SO(x_SO<(T-2)*Fs); 
x_SO        = x_SO(x_SO>    2*Fs); 

% Set the variables 
N_Stim      = length(x_SO);
Range_ERP   = [-1.25, 1.25];
Events      = zeros(length(time_events), N_Stim);
Events_FSP  = zeros(length(time_events), N_Stim);

% Segmentation
for i=1:N_Stim
    Events(:,i)     =  Vp    ((x_SO(i)+Range_ERP(1)*Fs)+1:(x_SO(i)+Range_ERP(2)*Fs+1));
    Events_FSP(:,i) =  Vp_FSP((x_SO(i)+Range_ERP(1)*Fs)+1:(x_SO(i)+Range_ERP(2)*Fs+1));
end

% Averaging
mean_ERP_model= mean(Events,    2); %#ok<*NASGU>
mean_FSP_model= mean(Events_FSP,2);
sd_ERP_model  = std (Events,    0, 2);
sd_FSP_model  = std (Events_FSP,0, 2);

% Save the data
save(savename, 'N_Stim', 'time_events', 'mean_ERP_data', 'mean_FSP_data', 'sem_ERP_data', 'sem_FSP_data', 'mean_ERP_model', 'mean_FSP_model', 'sd_ERP_model', 'sd_FSP_model');
end