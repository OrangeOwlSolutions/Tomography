clear all
close all
clc

N = 100;                                                % --- Number of input points

xin = linspace(0, pi, N);                               % --- Input sampling points
yin = sin(xin);                                         % --- Input samples

xout = pi * rand(1, N);                                 % --- Output sampling points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEAREST NEIGHBOR INTERPOLATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
youtNearestNeighbor         = nearestNeighborInterpolation1D(xin, yin, xout);
youtNearestNeighborMatlab   = interp1(xin, yin, xout, 'nearest').';
fprintf('RMS error Nearest Neighbor = %2.15f\n', 100 * sqrt(sum(abs(youtNearestNeighbor - youtNearestNeighborMatlab).^2)) / sum(abs(youtNearestNeighborMatlab).^2));

%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR INTERPOLATION %
%%%%%%%%%%%%%%%%%%%%%%%%
youtLinear                  = linearInterpolation1D(xin, yin, xout);
youtLinearMatlab            = interp1(xin, yin, xout, 'linear').';
fprintf('RMS error Linear = %2.15f\n', 100 * sqrt(sum(abs(youtLinear - youtLinearMatlab).^2)) / sum(abs(youtLinearMatlab).^2));

