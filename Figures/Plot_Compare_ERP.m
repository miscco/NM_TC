function Plot_Compare_ERP()
% This function creates Figure 7 of
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
load('Data/ERP_Average_data');
load('Data/ERP_double_model', 'Events', 'Events_FSP');

% Calculate average/sem
mean_ERP_model= mean(Events,    2);
mean_FSP_model= mean(Events_FSP,2);
sd_ERP_model  = std (Events,    0, 2);
sd_FSP_model  = std (Events_FSP,0, 2);

mean_ERP_data = mean_ERP;
mean_FSP_data = mean_FSP;
sem_ERP_data  = sem_ERP;
sem_FSP_data  = sem_FSP;

% Option array for set
Model_Range_ERP = [-75, -45];
Data_Range_ERP  = [-80, 50];
Model_Range_FSP = [-0.25, 1.25];
Data_Range_FSP  = [2, 8];
xRange          = -1:0.5:3;

Option_Name     = { 'YLim';
                    'YTick';
                    'YTickLabel';
                    'YColor';
                    'XTick';
                    'XLim'}';
            
Option_Model_ERP = {Model_Range_ERP;
                    -75:10:-40;
                    -75:10:-40;
                    'black';
                    xRange;
                    [xRange(1),xRange(end)]}'; %#ok<*NBRAK>

Option_Data_ERP  = {Data_Range_ERP;
                    -80:40:40;
                    -80:40:40;
                    'black';
                    xRange;
                    [xRange(1),xRange(end)]}';

Option_Model_FSP = {Model_Range_FSP;
                    linspace(Model_Range_FSP(1), Model_Range_FSP(2), 4);
                    linspace(Model_Range_FSP(1), Model_Range_FSP(2), 4);
                    'black';
                    xRange;
                    [xRange(1),xRange(end)]}';

Option_Data_FSP  = {Data_Range_FSP;
                    linspace(Data_Range_FSP(1), Data_Range_FSP(2), 4);
                    linspace(Data_Range_FSP(1), Data_Range_FSP(2), 4);
                    'black';
                    xRange;
                    [xRange(1),xRange(end)]}';

                
Lines = [0, 1.05; 0, 1.05];
                
% Define handle for plotting
BL_model =@(y,x) boundedline(y,x(:,1), x(:,2), 'alpha', 'transparency', 0.1, 'r');
BL_data  =@(y,x) boundedline(y,x(:,1), x(:,2), 'alpha', 'transparency', 0.1, 'black');

% Create figure
figure(1)
clf

% Create panel
p = panel('no-manage-font');
p.pack(2,1);

% set margins
p.de.margintop = 10;

% ERP
p(1,1).select();
[AX1, ~, ~] = plotyy(time_events,[mean_ERP_data, sem_ERP_data],time_events,[mean_ERP_model, sd_ERP_model], BL_data, BL_model);
set(AX1(1),Option_Name, Option_Data_ERP, 'box', 'off');
set(AX1(2),Option_Name, Option_Model_ERP);
line(Lines, get(AX1(1), 'YLim'), 'Color', 'black', 'LineStyle', ':');
line(Lines, get(AX1(2), 'YLim'), 'Color', 'black', 'LineStyle', ':');
ylabel(AX1(1),'EEG [\muV]');
ylabel(AX1(2),'V_{p} [mV]');

% Spindle power
p(2,1).select();
[AX2, ~, ~] = plotyy(time_events,[mean_FSP_data, sem_FSP_data],time_events,[mean_FSP_model, sd_FSP_model], BL_data, BL_model);
set(AX2(1),Option_Name, Option_Data_FSP);
set(AX2(2),Option_Name, Option_Model_FSP);
line(Lines, get(AX1(1), 'YLim'), 'Color', 'black', 'LineStyle', ':');
line(Lines, get(AX1(2), 'YLim'), 'Color', 'black', 'LineStyle', ':');
line(Lines, get(gca, 'YLim'), 'Color', 'black', 'LineStyle', ':');
ylabel(AX2(1),'Spindle Power [\muV^{2}]');
ylabel(AX2(2),'Spindle Power [mV^{2}cd ../NM_Th]');
xlabel('Time [s]');
end