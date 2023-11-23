% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 22/12/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name:       Input_PRACTISE
%   Purpose:    Organisation file for the input of PRACTISE
%   Comment:    This file defines and loads the input m-file. Edit the 
%               variable with the path and file name of the input m-file as
%               needed, but leave the run-expression untouched. Comment and
%               uncommenting line 20 and 22 here facilitates switching fast
%               between different input m-files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Define path and file name of the input m-file 
%     for configuration A 
Input_PRACTISE_m='output\PRACTISE\test_run'; % uncomment/comment in case
%     for configuration B 
% Input_PRACTISE_m='..\17112011_configB\Input_17112011_B'; % uncomment/comment in case

% Load the input m-file
Input_PRACTISE_m=Unix_Slash_PRACTISE(Input_PRACTISE_m);                    %KEEP
run(Input_PRACTISE_m)                                                      %KEEP