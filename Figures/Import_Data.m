function Import_Data(Type)
% This function imports the experimental data from Ngo et al 2013 utilized
% in
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

if nargin ==0
    Type = 1;
end

load('Data/Orig/Experimental_Data');

switch (Type)
    case 1; 
        Data = data.KC_Average; 
        sn = 'Data/KC_Average_data';    
    case 2; 
        Data = data.SO_Average; 
        sn = 'Data/SO_Average_data';   
    case 3; 
        Data = data.ERP_Average; 
        sn = 'Data/ERP_Average_data';   
end

time_events     = Data(:,1); %#ok<*NASGU>

mean_ERP        = Data(:,2);
mean_ERP_sham   = Data(:,3);

sem_ERP         = Data(:,4);
sem_ERP_sham    = Data(:,5);

mean_FSP        = Data(:,6);
mean_FSP_sham   = Data(:,7);

sem_FSP         = Data(:,8);
sem_FSP_sham    = Data(:,9);

save(sn, 'time_events', 'mean_ERP', 'mean_ERP_sham', 'sem_ERP', 'sem_ERP_sham', 'mean_FSP', 'mean_FSP_sham', 'sem_FSP', 'sem_FSP_sham');
end
