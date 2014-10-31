% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(type)
if nargin == 0
    type = 1;
end

if type == 1
    
    Param_Cortex        = [4.7;         % sigma_e
                           1.33;        % g_KNa
                           120E-3];     % dphi
                       
    Param_Thalamus      = [0.048;        % g_h
                           0.0245;        % g_LK_t
                           0.0245];        % g_LK_r
                       
else    
    Param_Cortex        = [6.05;         % sigma_e
                           2.1;          % g_KNa
                           120E-3];     % dphi
                       
    Param_Thalamus      = [0.052;       % g_h
                           0.02;       % g_LK_t
                           0.02];       % g_LK_r                    
end
                        
Connectivity            = [ 3;          % N_et
                            3;          % N_er
                            5;          % N_te
                            10];        % N_ti   

% stimulation parameters
% first number is the mode of stimulation
% 0 == none
% 1 == semi-periodic
% 2 == phase dependend
    
var_stim    = [ 0;          % mode of stimulation
                50;         % strength of the stimulus              in Hz (spikes per second)
                120;       	% duration of the stimulus              in ms
                5;          % time between stimulation events       in s  (ISI)
                0;          % range of ISI                          in s  [ISI-range,ISI+range]  
                1;          % Number of stimuli per event
                950;        % time between stimuli within a event   in ms         
                5];       % time until stimuli after minimum      in ms

T       	= 30;           % duration of the simulation

[Ve, Vt, Marker_Stim] = TC(T, Param_Cortex, Param_Thalamus, Connectivity, var_stim);

L        = length(Vt);
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), 
xlabel('Time in s'), 
ylabel('V_{e} in mV')
ylim([-80, -40])
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end

subplot(212), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), 
xlabel('Time in s'), 
ylabel('V_{t} in mV')
ylim([-70,-35]);
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end
[Pxx,f]   = pwelch(Ve-mean(Ve),hamming(5*L/T), 2*L/T, [], L/T);
n         = find(f<=30, 1, 'last' );

figure(2)
plot(f(1:n),log(Pxx(1:n)))
title('Powerspectrum with pwelch'), xlabel('frequency in Hz'), ylabel('Power (log)')
%save('Timeseries', 'Ve', 'Vt');
end