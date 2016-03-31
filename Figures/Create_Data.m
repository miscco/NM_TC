function Create_Data()
% This function creates the model data depicted in Figure 4-8 of
%
% A thalamocortical neural mass model of the EEG during NREM sleep and its
% response to auditory stimulation.
% M Schellenberger Costa, A Weigenand, H-VV Ngo, L Marshall, J Born, T Martinetz,
% JC Claussen.
% PLoS Computational Biology (in review).

% Move to source folder(assuming it contains the Figures folder
cd ..;

% Check if the executable exists and compile if needed
if(exist('Thalamus_mex.mesa64', 'file')==0)
    mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" TC_mex.cpp Cortical_Column.cpp Thalamic_Column.cpp;
end

% Add the path to the simulation routine
addpath(pwd);

% Go back into figures folder
cd Figures;

% NOTE: The data routines utilize various functions from the fieldtrip
% toolbox http://www.fieldtriptoolbox.org/
% You have to make sure, that the path to the fieldtrip routines is known
% Edit the below statement to fit to your system.
if(isempty(strfind(path, '/Tools/fieldtrip/preproc')))
    addpath('/Tools/fieldtrip/preproc');
end

% Import the original data from Ngo et al 2013
Import_Data(1);
Import_Data(2);
Import_Data(3);

% Time series
Data_Time_Series(1);
Data_Time_Series(2);

% Extract KC/SO averages
Data_SO_Average(1);
Data_SO_Average(2);
 
% Data of rhytmic two-click stimulation
Data_ERP_N3();
end