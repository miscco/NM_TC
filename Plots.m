% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" TC_mex.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(type)
if nargin == 0
    type = 2;
end

%mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" TC_mex.cpp Cortical_Column.cpp Thalamic_Column.cpp

if type == 1    
    Param_Cortex        = [4.7;         % sigma_p
                           1.33;        % g_KNa
                           2.];         % dphi
                       
    Param_Thalamus      = [0.034;       % g_LK
                           0.052];      % g_h                   
else    
    Param_Cortex        = [5.8;           % sigma_p
                           1.8;        % g_KNa
                           2.];          % dphi
                       
    Param_Thalamus      = [0.038;       % g_LK
                           0.052];      % g_h  
end
                        
Connectivity            = [ 2.5;        % N_tp
                            2.5;        % N_rp
                            5;          % N_pt
                            10];        % N_it    

% stimulation parameters
% first number is the mode of stimulation
% 0 == none
% 1 == semi-periodic
% 2 == phase dependend
    
var_stim    = [ 2;          % mode of stimulation
                200;         % strength of the stimulus              in Hz (spikes per second)
                100;       	% duration of the stimulus              in ms
                5;          % time between stimulation events       in s  (ISI)
                0;          % range of ISI                          in s  [ISI-range,ISI+range]  
                3;          % Number of stimuli per event
                1050;        % time between stimuli within a event   in ms         
                400];       % time until stimuli after minimum      in ms

T       	= 30;           % duration of the simulation

[Vp, Vt, Ca, ah, Marker_Stim] = TC_mex(T, Param_Cortex, Param_Thalamus, Connectivity, var_stim);

L        = length(Vt);
timeaxis = linspace(0,T,L);

figure(2)
subplot(411), plot(timeaxis,Vp)
title('Pyramidal membrane voltage'), 
xlabel('Time in s'), 
ylabel('V_{p} in mV')
xlim([0,30]);
ylim([-80, -40]);
% vertical line for markers
for i=1:var_stim(6)
    %hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'),'xdata', get(gca,'xlim'), 'color', 'black', 'LineStyle', ':');
    %changedependvar(hx,'x');
end

subplot(412), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), 
xlabel('Time in s'), 
ylabel('V_{t} in mV')
xlim([0,30]);
%ylim([-80,-35]);
% vertical line for markers


subplot(413), plot(timeaxis,ah)
title('Thalamic relay I_h activation'), 
xlabel('Time in s'), 
ylabel('m_h in mV')
xlim([0,30]);
ylim(get(gca,'ylim'));
% vertical line for markers


subplot(414), plot(timeaxis,Ca)
title('Thalamic relay ca concentration'), 
xlabel('Time in s'), 
ylabel('[Ca] in \mu mol')
xlim([0,30]);
ylim(get(gca,'ylim'));
% vertical line for markers

% [Pxx,f]   = pwelch(Ve-mean(Ve),hamming(10*L/T), 2*L/T, [], L/T);
% n         = find(f<=30, 1, 'last' );
% 
% figure(2)
% plot(f(1:n),log(Pxx(1:n)))
% title('Powerspectrum with pwelch'), xlabel('frequency in Hz'), ylabel('Power (log)')
%save('Timeseries', 'Ve', 'Vt');
end