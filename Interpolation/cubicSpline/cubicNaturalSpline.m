function [s0, s1, s2, s3] = cubicNaturalSpline(x, f)

% --- g(x) = sk0 + sk1 * (x - x(k)) + sk2 * (x - x(k))^2 + sk3 * (x -
% x(k))^3   in [xk_1, xk]
%
% --- The coefficents sk0, sk1, sk2, and sk3 are returned in the vectors
% s0,s1,s2, and s3 for each interval [xk_1, xk]

% --- Forces vectors to be columns
x = x(:); f = f(:); 

n = length(x) - 1;

h = x((1 : n) + 1) - x((0 : (n - 1)) + 1);     % --- +1 is due to Matlab's indexing; h has length n

% --- Constructing matrix A
lowerDiagonal = h(1 : (n - 1));
mainDiagonal  = 2 * (h(1 : (n - 1)) + h(2 : n));
upperDiagonal = h(2 : n);
A = spdiags([lowerDiagonal mainDiagonal upperDiagonal], [-1 0 1], n - 1, n - 1);

% --- Constructing Hf
d = (f((1 : n) + 1) - f((0 : (n - 1)) + 1)) ./ h; % --- +1 is due to Matlab's indexing
Hf = 6 * (d(2 : n) - d(1 : (n - 1)));

% --- Solving the linear system
% m = A \ Hf;
[L, U] = LUdecompositionCrout(A);
m = LUsolutionCrout(L, U, Hf).';

% --- Applying natural boundary conditions
m = [0; m; 0];

% --- Determining the cubic spline coefficients
s0 = f;
s1 = d - h .* (2 * m(1 : n) + m(2 : (n + 1))) / 6;
s2 = m / 2;
s3 =(m(2 : (n + 1)) - m(1 : n)) ./ (6 * h);
