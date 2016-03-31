function Plot_Compare_SO_Average()
% This function creates Figure 6 of
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

% Path to utility functions
if(isempty(strfind(path, [pwd, '/Tools'])))
    addpath([pwd, '/Tools']);
end

% Load the data
load('Data/KC_Average');
load('Data/SO_Average');


% Image settings
Data_Range_ERP      = [-120, 60];
xRange              = [-1.25,1.25];
xTicks              = -1:0.5:1;

Model_Range_ERP_N2  = [-75, -41];
Model_Ticks_ERP_N2  = linspace(-75,-45,4);
Model_Range_FSP_N2  = [-0.25, 2];
Data_Range_FSP_N2   = [-1, 5];

Model_Range_ERP_N3  = [-75, -45];
Model_Ticks_ERP_N3  = linspace(-75,-45,4);
Model_Range_FSP_N3  = [-0.5, 1.];
Data_Range_FSP_N3   = [2, 5];

mean_ERP_data_N2    = mean_ERP_data;
mean_ERP_model_N2   = mean_ERP_model;
mean_FSP_data_N2    = mean_FSP_data;
mean_FSP_model_N2   = mean_FSP_model;
sem_ERP_data_N2     = sem_ERP_data;
sem_FSP_data_N2     = sem_FSP_data;
sd_ERP_model_N2     = sd_ERP_model;
sd_FSP_model_N2     = sd_FSP_model;

mean_ERP_data_N3    = mean_ERP_data;
mean_ERP_model_N3   = mean_ERP_model;
mean_FSP_data_N3    = mean_FSP_data;
mean_FSP_model_N3   = mean_FSP_model;
sem_ERP_data_N3     = sem_ERP_data;
sem_FSP_data_N3     = sem_FSP_data;
sd_ERP_model_N3     = sd_ERP_model;
sd_FSP_model_N3     = sd_FSP_model;

% Define handle for plotting
BL_model =@(y,x) boundedline(y,x(:,1), x(:,2), 'alpha', 'transparency', 0.1, 'r');
BL_data  =@(y,x) boundedline(y,x(:,1), x(:,2), 'alpha', 'transparency', 0.1, 'black');

% Option array for set
Option_Name = { 'YLim';
                'YTick';
                'YColor';
                'XTick';
                'XLim'}';
Option_Data_ERP  = {Data_Range_ERP;
                    linspace(Data_Range_ERP(1), Data_Range_ERP(2), 4);
                    'black';
                    xTicks;
                    xRange}';

Option_Model_ERP_N2 = { Model_Range_ERP_N2;
                        Model_Ticks_ERP_N2;
                        'black';
                        xTicks;
                        xRange}'; %#ok<*NBRAK>

Option_Model_FSP_N2 = { Model_Range_FSP_N2;
                        linspace(Model_Range_FSP_N2(1), Model_Range_FSP_N2(2), 4);
                        'black';
                        xTicks;
                        xRange}';

Option_Data_FSP_N2  = { Data_Range_FSP_N2;
                        linspace(Data_Range_FSP_N2(1), Data_Range_FSP_N2(2), 4);
                        'black';
                        xTicks;
                        xRange}';

Option_Model_ERP_N3 = { Model_Range_ERP_N3;
                        Model_Ticks_ERP_N3;
                        'black';
                        xTicks;
                        xRange}'; %#ok<*NBRAK>

Option_Model_FSP_N3 = { Model_Range_FSP_N3;
                        linspace(Model_Range_FSP_N3(1), Model_Range_FSP_N3(2), 4);
                        'black';
                        xTicks;
                        xRange}';

Option_Data_FSP_N3  = { Data_Range_FSP_N3;
                        linspace(Data_Range_FSP_N3(1), Data_Range_FSP_N3(2), 4);
                        'black';
                        xTicks;
                        xRange}';

% Create figure
figure(1);
clf, shg

% Create panel
p = panel('no-manage-font');
p.pack(2,2);

% set margins
p.de.marginbottom = 10;
p.de.margintop    = 10;
p.de.marginleft   = 25;
p.de.marginright  = 25;

% ERP
p(1,1).select();
[AX1_N2, ~, ~] = plotyy(time_events,[mean_ERP_data_N2, sem_ERP_data_N2],time_events,[mean_ERP_model_N2, sd_ERP_model_N2], BL_data, BL_model);
set(AX1_N2(1), Option_Name, Option_Data_ERP, 'box', 'off', 'XTickLabel', []);
set(AX1_N2(2), Option_Name, Option_Model_ERP_N2, 'XTickLabel', []);
ylabel(AX1_N2(1),'EEG [\muV]');
ylabel(AX1_N2(2),'V_{p} [mV]');
p(1,1).title('(A) N2');

% Spindle power
p(2,1).select();
[AX2_N2, ~, ~] = plotyy(time_events,[mean_FSP_data_N2, sem_FSP_data_N2],time_events,[mean_FSP_model_N2, sd_FSP_model_N2], BL_data, BL_model);
set(AX2_N2(1), Option_Name, Option_Data_FSP_N2);
set(AX2_N2(2), Option_Name, Option_Model_FSP_N2);
xlabel('Time [s]');
ylabel(AX2_N2(1),'Spindle Power [\muV^{2}]');
ylabel(AX2_N2(2),'Spindle Power [mV^{2}]');

% ERP
p(1,2).select();
[AX1_N3, ~, ~] = plotyy(time_events,[mean_ERP_data_N3, sem_ERP_data_N3],time_events,[mean_ERP_model_N3, sd_ERP_model_N3], BL_data, BL_model);
set(AX1_N3(1), Option_Name, Option_Data_ERP, 'XTickLabel', []);
set(AX1_N3(2), Option_Name, Option_Model_ERP_N3, 'XTickLabel', []);
ylabel(AX1_N3(1),'EEG [\muV]');
ylabel(AX1_N3(2),'V_{p} [mV]');
p(1,2).title('(B) N3');

% Spindle power
p(2,2).select();
[AX2_N3, ~, ~] = plotyy(time_events,[mean_FSP_data_N3, sem_FSP_data_N3],time_events,[mean_FSP_model_N3, sd_FSP_model_N3], BL_data, BL_model);
set(AX2_N3(1), Option_Name, Option_Data_FSP_N3);
set(AX2_N3(2), Option_Name, Option_Model_FSP_N3);
xlabel('Time [s]');
ylabel(AX2_N3(1),'Spindle Power [\muV^{2}]');
ylabel(AX2_N3(2),'Spindle Power [mV^{2}]');
end