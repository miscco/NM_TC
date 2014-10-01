% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(T)

if nargin == 0
    
    Param_Cortex_N2     = [4.7;         % sigma_e
                           1.43;         % g_KNa
                           120E-3];     % dphi
    
    Param_Cortex_N3     = [6.3;          % sigma_e
                           2.;         % g_KNa
                           120E-3];     % dphi                      
                        
                        
    Param_Thalamus_N2   = [0.025;       % g_LK_t
                           0.025;       % g_LK_r
                           0.08];       % g_h
                        
                        
    Param_Thalamus_N3	= [0.021;       % g_LK_t
                           0.021;       % g_LK_r
                           0.08];       % g_h
                
    Connectivity        = [2.4;         % N_et
                           2.6;         % N_er
                           5;           % N_te
                           10];         % N_ti   

    % stimulation parameters
    % first number is the mode of stimulation
    % 0 == none
    % 1 == semi-periodic
    % 2 == phase dependend
    
    var_stim    = [ 2;          % mode of stimulation
                    25;         % strength of the stimulus              in Hz (spikes per second)
                    120;       	% duration of the stimulus              in ms
                    5;          % time between stimulation events       in s  (ISI)
                    0;          % range of ISI                          in s  [ISI-range,ISI+range]  
                    2;          % Number of stimuli per event
                    950;        % time between stimuli within a event   in ms         
                    500];       % time until stimuli after minimum      in ms

    T       	= 30;           % duration of the simulation
end

[Ve, Vt, Marker_Stim] = TC(T, Param_Cortex_N2, Param_Thalamus_N2, Connectivity, var_stim);
%[Ve, Vt, Marker_Stim] = TC(T, Param_Cortex_N3, Param_Thalamus_N3, Connectivity, var_stim);

L        = max(size(Vt));
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), xlabel('time in s'), ylabel('Ve in mV')
ylim([-80, -40])
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end
hx2 = graph2d.constantline((Marker_Stim/1E2 -var_stim(8)/1E3), 'color', 'red');
changedependvar(hx2,'x');

subplot(212), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), xlabel('time in s'), ylabel('Vt in mV')
ylim(get(gca,'ylim'));
% vertical line for markers
for i=1:var_stim(6)
    hx = graph2d.constantline(Marker_Stim/1E2+(i-1)*var_stim(7)/1E3,'ydata', get(gca,'ylim'), 'color', 'black', 'LineStyle', ':');
    changedependvar(hx,'x');
end
% [Pxx,f]   = pwelch(Ve-mean(Ve),hamming(L/30), 4*L/T, 2048, L/T);
% n         = find(f<=30, 1, 'last' );
% 
% figure(2)
% plot(f(1:n),log(Pxx(1:n)))
% title('Powerspectrum with pwelch'), xlabel('frequency in Hz'), ylabel('Power (log)')
%save('Timeseries', 'Ve', 'Vt');
end