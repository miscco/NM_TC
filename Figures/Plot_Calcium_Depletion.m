function Plot_Calcium_Depletion()
% This function creates Figure 8 of
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
 
% Load the endogenous data
load('Data/Time_Series_N3.mat');
Vp_endo = Vp(10900:11400);
Vt_endo = Vt(10900:11400);
ah_endo = ah(10900:11400);
timeaxis = linspace(0,5,501);

Range = [-1.25, 3.75];
Fs    = 100;
i     = 5;

% Load the stimulationd data
load('Data/Ca_Depletion.mat');
Vp_stim = Vp((Marker_Stim(i)+Range(1)*Fs)+1:(Marker_Stim(i)+Range(2)*Fs+1));
Vt_stim = Vt((Marker_Stim(i)+Range(1)*Fs)+1:(Marker_Stim(i)+Range(2)*Fs+1));
ah_stim = ah((Marker_Stim(i)+Range(1)*Fs)+1:(Marker_Stim(i)+Range(2)*Fs+1));

Lines = [1.26, 2.31; 1.26, 2.31];

% Create figure
figure(2)
clf

% Create panel
p = panel('no-manage-font');
p.pack(3,2);

% Set margins
p.de.marginbottom = 10;
p.de.margintop    = 5;
p.de.marginleft   = 5;
p.de.marginright  = 5;

% Plot the data
p(1,1).select();
plot(timeaxis, Vp_endo, 'color', 'black');
set(gca, 'XTick', 1:4, 'XTickLabel', {'','','',''}, 'YTick', [-80, -60, -40], 'YLim', [-80, -40], 'XLim', [0,5]);
p(1,1).ylabel('V_{p} [mV]');
p(1,1).title('(A) Endogenous SOs');

p(2,1).select(); 
plot(timeaxis, Vt_endo, 'color', 'black');
set(gca, 'XTick', 1:4, 'XTickLabel', {'','','',''}, 'YTick', [-70, -55, -40], 'YLim', [-70, -40], 'XLim', [0,5]);
p(2,1).ylabel('V_{t} [mV]');

p(3,1).select(); 
plot(timeaxis, ah_endo.*51, 'color', 'black');
set(gca, 'XTick', 1:4, 'XTickLabel', {'','','',''}, 'YTick', [15, 18, 21], 'YLim', [15, 21], 'XLim', [0,5]);
p(3,1).xlabel('Time [s]'); 
p(3,1).ylabel('g_{h} [\muS/cm^{2}]');

p(1,2).select();
plot(timeaxis, Vp_stim, 'color', 'black');
set(gca, 'XTick', 1:4, 'XTickLabel', {'','','',''}, 'YAxisLocation', 'Right', 'YTick', [-80, -60, -40], 'YLim', [-80, -40], 'XLim', [0,5]);
line(Lines, get(gca, 'YLim'), 'color', 'black', 'LineStyle', ':');
p(1,2).ylabel('V_{p} [mV]');
p(1,2).title('(B) Closed loop stimulation');

p(2,2).select(); 
plot(timeaxis, Vt_stim, 'color', 'black');
set(gca, 'XTick', 1:4, 'XTickLabel', {'','','',''}, 'YAxisLocation', 'Right', 'YTick', [-70, -55, -40], 'YLim', [-70, -40], 'XLim', [0,5]);
line(Lines, get(gca, 'YLim'), 'color', 'black', 'LineStyle', ':');
p(2,2).ylabel('V_{t} [mV]');

p(3,2).select(); 
plot(timeaxis, ah_stim.*51, 'color', 'black');
set(gca, 'XTick', 1:4, 'XTickLabel', {'','','',''}, 'YAxisLocation', 'Right', 'YTick', [15, 18, 21], 'YLim', [15, 21], 'XLim', [0,5]);
line(Lines, get(gca, 'YLim'), 'color', 'black', 'LineStyle', ':');
p(3,2).xlabel('Time [s]'); 
p(3,2).ylabel('g_{h} [\muS/cm^{2}]');

end