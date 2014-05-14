% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=gnu++0x -fpermissive" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(T)

if nargin == 0
    Input_N3    = [ 8.7;          % sigma_e
                    2.6;        % alpha_Na
                    3;          % tau_Na
                    1.6;        % g_KNa
                    60E-3];     % dphi
                        
                        
    Input_N2    = [ 4.6;        % sigma_e
                    2;          % alpha_Na
                    1.2;        % tau_Na
                    1.33;       % g_KNa
                    60E-3];     % dphi

    Con     	= [0;           % N_et
               	   0;           % N_er
                   0;           % N_te
                   0];          % N_ti    

    % stimulation parameters
    % first number is the mode of stimulation
    % 0 == none
    % 1 == periodic
    % 2 == phase dependend up state
    % 3 == phase dependend down state
    
    var_stim    = [ 0           % mode of stimulation
                    100.0;      % strength of the stimulus      in Hz (spikes per second)
                    100;       	% duration of the stimulus      in ms
                    8;          % time between stimuli          in s    
                    550];       % time until stimuli after min 	in ms

    T       	= 30;           % duration of the simulation
end

[Ve, Vt] = TC(T, Input_N3, Con, var_stim);

L        = max(size(Vt));
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), xlabel('time in s'), ylabel('Ve in mV')
subplot(212), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), xlabel('time in s'), ylabel('Vt in mV')
%save('TC.mat','Ve','Vt')
