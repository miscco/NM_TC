% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(T)

if nargin == 0
    Input_N3    = [ 6.5;        % sigma_e
                    2.1;        % g_KNa
                    120E-3];	% dphi
                        
                        
    Input_N2    = [ 4.7;        % sigma_e
                    1.5;        % g_KNa
                    120E-3];	% dphi
                
    Connectivity= [2.4;         % N_et
               	   2.5;         % N_er
                   5;           % N_te
                   5];          % N_ti   

    % stimulation parameters
    % first number is the mode of stimulation
    % 0 == none
    % 1 == periodic
    % 2 == phase dependend up state
    % 3 == phase dependend down state
    
    var_stim    = [ 0;          % mode of stimulation
                    60;          % strength of the stimulus             in Hz (spikes per second)
                    100;       	% duration of the stimulus              in ms
                    5;          % time between stimuli                  in s    
                    650];       % time until stimuli after negativ peak	in ms

    T       	= 30;           % duration of the simulation
end

[Ve, Vt, Marker_Stim] = TC(T, Input_N2, Connectivity, var_stim);

L        = max(size(Vt));
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), xlabel('time in s'), ylabel('Ve in mV')
% vertical line for markers
%hx = graph2d.constantline(Marker_Stim(3,:));
%changedependvar(hx,'x');

subplot(212), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), xlabel('time in s'), ylabel('Vt in mV')
% vertical line for markers
%hx = graph2d.constantline(Marker_Stim(3,:));
%changedependvar(hx,'x');


% [Pxx,f]   = pwelch(Ve-mean(Ve),[], [], [], L/T);
% n         = find(f<=30, 1, 'last' );
% 
% figure(2)
% plot(f(1:n),log(Pxx(1:n)))
% title('Powerspectrum with pwelch'), xlabel('frequency in Hz'), ylabel('Power (log)')
% save('Timeseries', 'Ve', 'Vt');
end