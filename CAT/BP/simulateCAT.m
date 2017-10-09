function simulateCAT(dataFileName, rho1, rho2, numProjections, numDetectors)

fontSize = 18;

% --- Each row describes an ellipse.
% --- At points where ellipses overlap, the refractive indices are cumulated.

% --- Shepp-Logan phantom
% --- From Murrell, "Computer-Aided Tomography," The Mathematical Journal (1996), p. 60

% -----------------------------------------------------------
%    center          major   minor   rotation   refractive
%    coordinate      axis    axis    angle      index
%     xtemp      ytemp                       (deg)
AI = [0      0       0.92    0.69     90        2.0
      0     -0.0184  0.874   0.6624   90       -0.9
      0.22   0       0.31    0.11     72       -0.1
     -0.22   0       0.41    0.16    108       -0.1
      0      0.35    0.25    0.21     90        0.3
      0      0.1     0.046   0.046    0         0.3
      0     -0.1     0.046   0.046    0         0.3
     -0.08  -0.605   0.046   0.023    0         0.3
      0     -0.605   0.023   0.023    0         0.3
      0.06  -0.605   0.046   0.023    90        0.3];

x1 = AI(:, 1);                  % --- Centre coordinates xtemp
y1 = AI(:, 2);                  % --- Centre coordinates ytemp
A  = AI(:, 3);                  % --- Semi major axes
B  = AI(:, 4);                  % --- Semi minor axes
a1 = AI(:, 5);                  % --- Rotation angles
ri = AI(:, 6);                  % --- Refractive indices
numEllipses = length(x1);       % --- Number of ellipses

% --- Image sampling
xmin    = -1;                   % --- Left side
xmax    =  1;                   % --- Right side
ymin    = -1;                   % --- Bottom side
ymax    =  1;                   % --- Upper side
Deltax  = 0.01;                 % --- Sampling step along xtemp
Deltay  = 0.01;                 % --- Sampling step along ytemp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE REFRACTIVE INDICES MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = xmin : Deltax : xmax;       % --- Sampling points along xtemp
y = ymin : Deltay : ymax;       % --- Sampling points along ytemp
Nx = length(x);                 % --- Number of sampling points along xtemp
Ny = length(y);    
[XX, YY] = meshgrid(x, -y);
f0 = zeros(Ny, Nx);             % --- Refractive index matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATING THE REFERENCE IMAGE: SUMMING UP THE CONTRIBUTIONS FROM ALL THE ELLIPSES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : numEllipses,
  alpha1 = a1(k) / 180 * pi;
  xtemp =  (XX - x1(k)) * cos(alpha1) + (YY - y1(k)) * sin(alpha1);
  ytemp = -(XX - x1(k)) * sin(alpha1) + (YY - y1(k)) * cos(alpha1);
  f0 = f0 + ri(k) * heaviside(1 - (xtemp / A(k)).^2 - (ytemp / B(k)).^2);
end
fmax = max(max(f0));

close all
figure
imagesc(x, -y, f0 / fmax, [0 1])
set(gca, 'YDir', 'Normal')
colormap(gray)
axis square
axis([xmin xmax ymin ymax])
set(gca, 'FontSize', fontSize, 'FontWeight', 'b')
title('Reference image')
xlabel('x', 'FontSize', fontSize, 'FontWeight', 'b')
ylabel('y', 'FontSize', fontSize, 'FontWeight', 'b')

print('-djpeg', 'phantom.jpg', '-r360');

% --- Angles at which the projections are taken - Uniform partition of the
% interval (0, pi)
angles = (1 : numProjections)' / numProjections * pi;

% --- Inter-detector distance
drho = (rho2 - rho1) / (numDetectors - 1);
% --- Detector detectorPositions
rho  = rho1 + (0 : (numDetectors - 1)) * drho;
detectorPositions = rho;

%%%%%%%%%%%%%%%%%%%%%%%
% GENERATING SINOGRAM %
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Generating projections...\n');

sinogram = zeros(numProjections, numDetectors);
for currentProjection = 1 : numProjections,         % --- Loop over the projections
    
    p1 = zeros(1, numDetectors);
    
    for k = 1 : numEllipses,                        % --- Loop over the ellipses
        theta = angles(currentProjection) - a1(k)/180*pi;
        
        a = sqrt((A(k) * cos(theta))^2 + (B(k) * sin(theta))^2);
        
        s1 = sqrt(x1(k)^2 + y1(k)^2);
        gamma1 = atan2(y1(k), x1(k));
        
        for currentDetector = 1 : numDetectors,
            t1 = rho(currentDetector) - s1 * cos(gamma1 - angles(currentProjection));
            if (abs(t1) < a),
                p1(currentDetector) = p1(currentDetector) + 2 * ri(k) * A(k) * B(k) / a^2 * sqrt(a^2 - t1^2);
            end
        end
    end
    sinogram(currentProjection, :) = p1;
end

%%%%%%%%%%
% FIGURE %
%%%%%%%%%%
figure
imagesc(rho, angles / pi * 180, sinogram / max(max(sinogram)), [0 1]);
set(gca,'YDir','Normal')
colormap(gray)
axis square
title('Sinogram', 'FontSize', fontSize)
set(gca, 'FontSize', fontSize, 'FontWeight', 'b')
xlabel('\rho','FontSize', fontSize, 'FontWeight', 'b')
ylabel('\theta (deg)','FontSize', fontSize, 'FontWeight', 'b')

print('-djpeg', 'sinogram.jpg', '-r360');

save(dataFileName, 'numProjections', 'angles', 'numDetectors', 'detectorPositions', 'sinogram');

