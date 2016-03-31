function Plot_Time_Series(type)
% This function creates Figure 4 and 5 of
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

if nargin==0  
    type = 1;
end


% Path to utility functions
if(isempty(strfind(path, [pwd, '/Tools'])))
    addpath([pwd, '/Tools']);
end
                
if type == 1   
    Name   = 'N2';
else
    Name   = 'N3';
end

% load the data
load(['Data/Time_Series_Short_',Name]);

T        = 30;
L        = max(size(Vt));
timeaxis = linspace(0,T,L);

% Create figure
figure(1)
clf

% Create panel
p = panel('no-manage-font');
p.pack(2,1);

% set margins
p.de.margintop = 10;

p(1,1).select();
plot(timeaxis, Vp, 'Color', 'black');
ylabel(gca, 'V_{p} [mV]');
set(gca, 'XTick', [], 'YTick', [-80, -60, -40], 'YLim', [-80, -40]);

% thalamic relay membrane voltage
p(2,1).select(); 
plot(timeaxis, Vt, 'Color', 'black');
p(2,1).xlabel('Time [s]'); 
ylabel(gca, 'V_{t} [mV]');
switch (type)
    case 1
    ylim([-70, -30]);
    set(gca, 'ytick', [-70, -50, -30]);
    case 2
    ylim([-70, -40]);
    set(gca, 'ytick', [-70, -55, -40]);
end
end