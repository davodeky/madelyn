% Pharmacokinetics in space

%source code:
%https://ascpt.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpsp4.12399&file=psp412399-sup-0002-Code.pdf 
%from https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12399 

%Main function: Calls each function of a PBPK model
clear all; %clear all variables from previous simulations
close all; %close all windows from previous simulations
clc; %clear command window
DefineParameters; %model structure and parameters are defined
Population; %the defined population is generated
Drug; %load drug files / in vitro-to-in vivo extrapolation
SolveODE; %solve the ordinary differential equations
PostProcessing; %process the data from the ODE solution / output results