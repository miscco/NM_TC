% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(type)
if nargin == 0
    type = 2;
end

if type == 1
    
    Param_Cortex        = [4.7;         % sigma_e
                           1.35;        % g_KNa
                           120E-3];     % dphi
                       
    Param_Thalamus      = [0.037;        % g_h
                           0.024;        % g_LK_t
                           0.024];        % g_LK_r
                       
else    
    Param_Cortex        = [6.;         % sigma_e
                           2.1;          % g_KNa
                           120E-3];     % dphi
                       
    Param_Thalamus      = [0.042;       % g_h
                           0.022;       % g_LK_t
                           0.022;       % g_LK_r  
                           2.5;         % k1
                           4;           % k2
                           1;           % k3
                           1];          % k4
end
                        
Connectivity            = [ 4;          % N_et
                            4;          % N_er
                            5;          % N_te
                            8];        % N_ti   

% stimulation parameters
% first number is the mode of stimulation
% 0 == none
% 1 == semi-periodic
% 2 == phase dependend
    
var_stim    = [ 2;          % mode of stimulation
                80;         % strength of the stimulus              in Hz (spikes per second)
                150;       	% duration of the stimulus              in ms
                5;          % time between stimulation events       in s  (ISI)
                0;          % range of ISI                          in s  [ISI-range,ISI+range]  
                3;          % Number of stimuli per event
                1050;        % time between stimuli within a event   in ms         
                450];       % time until stimuli after minimum      in ms

T       	= 30;           % duration of the simulation

[Ve, Vt, ah, Marker_Stim] = TC(T, Param_Cortex, Param_Thalamus, Connectivity, var_stim);

L        = length(Vt);
timeaxis = linspace(0,T,L);

figure(1)
subplot(311), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), 
xlabel('Time in s'), 
ylabel('V_{e} in mV')
xlim([0,30]);
ylim([-80, -40]);
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'),'xdata', get(gca,'xlim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end

subplot(312), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), 
xlabel('Time in s'), 
ylabel('V_{t} in mV')
xlim([0,30]);
ylim([-70,-35]);
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'),'xdata', get(gca,'xlim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end

subplot(313), plot(timeaxis,ah)
title('Thalamic relay I_h activation'), 
xlabel('Time in s'), 
ylabel('m_h in mV')
xlim([0,30]);
ylim(get(gca,'ylim'));
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'),'xdata', get(gca,'xlim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end
% [Pxx,f]   = pwelch(Ve-mean(Ve),hamming(5*L/T), 2*L/T, [], L/T);
% n         = find(f<=30, 1, 'last' );
% 
% figure(2)
% plot(f(1:n),log(Pxx(1:n)))
% title('Powerspectrum with pwelch'), xlabel('frequency in Hz'), ylabel('Power (log)')
%save('Timeseries', 'Ve', 'Vt');
end