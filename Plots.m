% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=gnu++0x -fpermissive" TC.cpp Cortical_Column.cpp Thalamic_Column.cpp

function Plots(T, onset)

if nargin == 0
    Input_N3    = [ 2.6;        % alpha_Na
                    3;          % tau_Na
                    1.33	% g_KNa
                    -63;        % theta_e
                    8;          % sigma_e
                    30E-3];     % dphi
                        
                        
    Input_N2    = [ 2;          % alpha_Na
                    1;          % tau_Na
                    1.33;	% g_KNa
                    -58.5;      % theta_e
                    4.25;       % sigma_e
                    30E-3];     % dphi

    Con     	= [10;		% N_tr
               	   20;		% N_rt
                   40];		% N_rr    

    var_stim    = [ 0;          % strength of the stimulus 	in Hz (spikes per second)
                    0;          % time between stimuli 		in s    
                    0;          % time until first stimuli 	in s
                    0];        	% duration of the stimulus 	in ms

    T       	= 30;  		% duration of the simulation
end

[Ve, Vt] = TC(T, Input_N2, Con, var_stim);

L        = max(size(Vt));
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,Ve)
title('Pyramidal membrane voltage'), xlabel('time in s'), ylabel('Ve in mV')
subplot(212), plot(timeaxis,Vt)
title('Thalamic relay membrane voltage'), xlabel('time in s'), ylabel('Vt in mV')
%save('TC.mat','Ve','Vt')
