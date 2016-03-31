function Data_Time_Series(type)
% This function creates the model data depicted in Figure 4-8 of
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

% To ensure availability of the simulation routine and the utility functions 
% from fieldtrip it should be called from Create_Data.m

% Load the parameter settings
if type == 1    
    load('Data/Parameter_N2');           
    Protocol_Name = 'N2';                              
else    
    load('Data/Parameter_N3');       
    Protocol_Name = 'N3';                   
end   
 
% there is no stimulation so set var_stim to zero
var_stim            = zeros(8,1);

% Duration of the simulation
T                   = 3600;         

% Run the simulation (Can take some time)
[Vp, Vt, Ca, ah, ~] = TC_mex(T, Param_Cortex, Param_Thalamus, Connectivity, var_stim);

% Save the full data
save(['Data/Time_Series_',Protocol_Name], 'Vp', 'Vt', 'Ca', 'ah');

% Save a smaller snipplet for example time series plot
Vp = Vp(1:3000);
Vt = Vt(1:3000);
Ca = Ca(1:3000);
ah = ah(1:3000);
save(['Data/Time_Series_Short_',Protocol_Name], 'Vp', 'Vt', 'Ca', 'ah');
end