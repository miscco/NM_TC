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
    % 1 == periodic
    % 2 == phase dependend up state
    % 3 == phase dependend down state
    
    var_stim    = [ 0;          % mode of stimulation
                    60;          % strength of the stimulus             in Hz (spikes per second)
                    200;       	% duration of the stimulus              in ms
                    5;          % time between stimuli                  in s    
                    300];       % time until stimuli after negativ peak	in ms

    T       	= 600;           % duration of the simulation
end

%[Ve, Vt, Marker_Stim] = TC(T, Param_Cortex_N2, Param_Thalamus_N2, Connectivity, var_stim);
[Ve, Vt, Marker_Stim] = TC(T, Param_Cortex_N3, Param_Thalamus_N3, Connectivity, var_stim);

L        = max(size(Vt));
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), xlabel('time in s'), ylabel('Ve in mV')
ylim([-80, -40])
% vertical line for markers
hx1 = graph2d.constantline(Marker_Stim(2,:), 'color', 'black');
hx2 = graph2d.constantline(Marker_Stim(1,:), 'color', 'red');
changedependvar(hx1,'x');
changedependvar(hx2,'x');

subplot(212), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), xlabel('time in s'), ylabel('Vt in mV')
% vertical line for markers
hx = graph2d.constantline(Marker_Stim(2,:));
changedependvar(hx,'x');

[Pxx,f]   = pwelch(Ve-mean(Ve),hamming(L/30), 4*L/T, 2048, L/T);
n         = find(f<=30, 1, 'last' );

figure(2)
plot(f(1:n),log(Pxx(1:n)))
title('Powerspectrum with pwelch'), xlabel('frequency in Hz'), ylabel('Power (log)')
%save('Timeseries', 'Ve', 'Vt');
end