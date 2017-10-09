clear all
close all
clc

num_points = 50;

runge = @(x) 1 ./ (1 + 25 * x.^2);

x = linspace(-1, 1, num_points);
y = runge(x);

[s0, s1, s2, s3] = cubicNaturalSpline(x, y);

xout = linspace(-1, 1, 3 * num_points);

yout = constructCubicSpline(x, xout, s0, s1, s2, s3);

% ymatlab = interp1(x, y, xout, 'spline'); does not provide natural spline
% due to different boundary conditions
% ymatlab = spline(x, y, xout); does not provide natural spline
% due to different boundary conditions
pp = csape(x, y, 'second');
ymatlab = ppval(pp, xout);

figure
plot(x, y, '*');
hold on
plot(xout, yout, 'r');
