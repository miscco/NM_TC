% mex command is given by: 

function Test_Stimulation(type)
if nargin == 0
    type = 4;
end


mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3 -lgopm" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

% Path to fieltrip preprocessing function
if(isempty(strfind(path, '/nfshome/schellen/Documents/MATLAB/Tools/fieldtrip/preproc')))
    addpath('~/Documents/MATLAB/Tools/fieldtrip/preproc');
end

% Path to helper function
if(isempty(strfind(path, '/nfshome/schellen/Documents/MATLAB/Tools/boundedline')))
    addpath('~/Documents/MATLAB/Tools/boundedline');
end
  
Param_Cortex        = [6;          % sigma_e
                       2.05;         % g_KNa
                       120E-3];     % dphi
                       
Param_Thalamus      = [0.052;       % g_h
                       0.02];       % g_LK

Connectivity        = [ 2.6;        % N_et
                        2.6;        % N_er
                        5;          % N_te
                        10];        % N_ti   
                       
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
                1075;        % time between stimuli within a event   in ms         
                450];       % time until stimuli after minimum      in ms

T       	= 3600;           % duration of the simulation

load('/nfshome/schellen/Documents/MATLAB/TC_model/Data/ERP_Average_data');

Model_Range_ERP = [-75, -45]; 
Data_Range_ERP  = [-80, 50]; 
Model_Range_FSP = [-0.25, 1.25]; 
Data_Range_FSP  = [2, 8]; 
xRange          = -1:0.5:3;

% Option array for set
Option_Name     = { 'ylim';
                    'ytick';
                    'yticklabel';
                    'ycolor';
                    'xtick';
                    'xlim'}';
            
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

[Ve, Vt, Ca, ah, Marker_Stim] = TC(T, Param_Cortex, Param_Thalamus, Connectivity, var_stim);
Fs          = length(Ve)/T;
Ve_FSP      = ft_preproc_hilbert(ft_preproc_bandpassfilter(Ve, Fs, [12, 15], 513, 'fir'), 'abs').^2;
xRange      = [-1, 3]; 

% Search for peaks
x_SO        = Marker_Stim;

% Remove those events, that are too close to begin/end
x_SO        = x_SO(x_SO<(T-xRange(end))*Fs); 
x_SO        = x_SO(x_SO>  -xRange(1)*Fs);

% Set the variables 
N_Stim      = length(x_SO);
time_event  = linspace(xRange(1), xRange(end), (xRange(end)-xRange(1))*Fs+1);
Events      = zeros(length(time_event), N_Stim);
Events_FSP  = zeros(length(time_event), N_Stim);

% Segmentation
for i=1:N_Stim
    Events(:,i)     =  Ve    ((x_SO(i)+xRange(1)*Fs):(x_SO(i)+xRange(end)*Fs));
    Events_FSP(:,i) =  Ve_FSP((x_SO(i)+xRange(1)*Fs):(x_SO(i)+xRange(end)*Fs));
end

mean_ERP_model= mean(Events,    2); %#ok<*NASGU>
mean_FSP_model= mean(Events_FSP,2);
sd_ERP_model  = std (Events,    0, 2);
sd_FSP_model  = std (Events_FSP,0, 2);

% Define handle for plotting
BL_model =@(y,x) boundedline(y,x(:,1), x(:,2), 'alpha', 'transparency', 0.1, 'r');
BL_data  =@(y,x) boundedline(y,x(:,1), x(:,2), 'alpha', 'transparency', 0.1, 'black');
    
figure(1)
subplot(411)
plot(linspace(0,30,3000),Ve(101:3100));
title(['Ve with a mean of :',num2str(mean(Ve))]);
subplot(412)
plot(linspace(0,30,3000),Vt(101:3100));
title(['Vt with a mean of :',num2str(mean(Vt))]);
subplot(413)
plot(linspace(0,30,3000),Ca(101:3100));
title(['Ca with a mean of :',num2str(mean(Ca))]);
subplot(414)
plot(linspace(0,30,3000),ah(101:3100));
title(['ah with a mean of :',num2str(mean(ah))]);
                
% Create figure
figure(2)
clf
subplot(211)
[AX1, ~, ~] = plotyy(time_events,[mean_ERP, sem_FSP],time_events,[mean_ERP_model, sd_ERP_model], BL_data, BL_model);
set(AX1(1),Option_Name, Option_Data_ERP, 'box', 'off');
set(AX1(2),Option_Name, Option_Model_ERP);
ylabel(AX1(1),'EEG [$\mu$V]');
ylabel(AX1(2),'V$_{p}$ [mV]');

subplot(212)
[AX2, ~, ~] = plotyy(time_events,[mean_FSP, sem_FSP],time_events,[mean_FSP_model, sd_FSP_model], BL_data, BL_model);
set(AX2(1),Option_Name, Option_Data_FSP);
set(AX2(2),Option_Name, Option_Model_FSP);
ylabel(AX2(1),'Spindle Power [$\mu$V$^{2}$]');
ylabel(AX2(2),'Spindle Power [mV$^{2}$]');
title([num2str(N_Stim), ' Events'])

% Marker for stimulation
for i=1:2
    hx1 = graph2d.constantline((i-1)*1.05+0.125*(i-1)*(i-2)/2,'ydata', get(AX1(1),'ylim'), 'parent', AX1(1), 'color', 'black', 'LineStyle', ':');
    hx2 = graph2d.constantline((i-1)*1.05+0.125*(i-1)*(i-2)/2,'ydata', get(AX2(1),'ylim'), 'parent', AX2(1), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx1,'x');
    changedependvar(hx2,'x');
end

end