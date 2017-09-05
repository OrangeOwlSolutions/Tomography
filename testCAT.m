clear all
close all
clc

rho1 = -1;
rho2 =  1;
numProjections  = 110;
numDetectors    = 127;

dataFileName = 'dataFileName.mat';   
simulateCAT(dataFileName, rho1, rho2, numProjections, numDetectors);

%%%%%%%%%%%%%%%%%%%
% RECONSTRUCTIONS %
%%%%%%%%%%%%%%%%%%%
x0 = rho1 : .01 : rho2;
y0 = rho1 : .01 : rho2;

FBP(dataFileName, 'ReconstructionFPB', x0, y0);
ConvBP(dataFileName, 'ReconstructionConvPB', x0, y0);

