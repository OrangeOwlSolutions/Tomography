function Zi = nearestNeighborInterpolation2D(X, Y, Z, xi, yi, methodflag)
%QINTERP2 2-dimensional fast interpolation
% qinterp2 provides a speedup over interp2 in the same way that
% qinterp1 provides a speedup over interp1
%
% Usage:
%   yi = qinterp2(X,Y,Z,xi,yi)  - Same usage as interp2, X and Y should be
%                                 "plaid" (e.g., generated with meshgrid).
%                                 Defaults to bilinear interpolation
%   yi = qinterp2(...,flag)
%           flag = 0       - Nearest-neighbor
%           flag = 1       - Triangular-mesh linear interpolation.
%           flag = 2       - Bilinear (equivalent to MATLAB's 'linear')
%
% Usage restrictions
%   X(:,n) and Y(m,:) must be monotonically and evenly increasing
%   e.g.,  [X,Y] = meshgrid(-5:5,0:0.025:1);
%
% Examples:
%   % Set up the library data for each example
%   [X,Y] = meshgrid(-4:0.1:4,-4:0.1:4);
%   Z = exp(-X.^2-Y.^2);
%
%   % Interpolate a line
%   xi = -4:0.03:4; yi = xi;
%   Zi = qinterp2(X,Y,Z,xi,yi);
%   % Plot the interpolant over the library data
%   figure, mesh(X,Y,Z), hold on, plot3(xi,yi,Zi,'-r');
%
%   % Interpolate a region
%   [xi,yi] = meshgrid(-3:0.3:0,0:0.3:3);
%   Zi = qinterp2(X,Y,Z,xi,yi);
%   % Plot the interpolant
%   figure, mesh(X,Y,Zi);
%
% Error checking
%   WARNING: Little error checking is performed on the X or Y arrays. If these
%   are not proper monotonic, evenly increasing plaid arrays, this
%   function will produce incorrect output without generating an error.
%   This is done because error checking of the "library" arrays takes O(mn)
%   time (where the arrays are size [m,n]).  This function is
%   algorithmically independent of the size of the library arrays, and its
%   run time is determine solely by the size of xi and yi
%
% Using with non-evenly spaced arrays:
%   See qinterp1

% Search array error checking
if size(xi)~=size(yi)
    error('%s and %s must be equal size',inputname(4),inputname(5));
end

% Library array error checking (size only)
if size(X)~=size(Y)
    error('%s and %s must have the same size',inputname(1),inputname(2));
end
librarySize = size(X);

% Decide the interpolation method
% if nargin>=6
%     method = methodflag;
% else
%     method = 2; % Default to bilinear
% end

% Get X and Y library array spacing
ndx = 1/(X(1,2)-X(1,1));    ndy = 1/(Y(2,1)-Y(1,1));
% Begin mapping xi and yi vectors onto index space by subtracting library
% array minima and scaling to index spacing
xi = (xi - X(1,1))*ndx;       yi = (yi - Y(1,1))*ndy;

% Fill Zi with NaNs
Zi = NaN*ones(size(xi));

        % Find the nearest point in index space
        rxi = round(xi)+1;  ryi = round(yi)+1;
        % Find points that are in X,Y range
%         flag = rxi>0 & rxi<=librarySize(2) & ~isnan(rxi) &...
%             ryi>0 & ryi<=librarySize(1) & ~isnan(ryi);
        flag = rxi < 1 | rxi > librarySize(2) | ryi < 1 | ryi > librarySize(1);
%         flagx = rxi < 1 | rxi > librarySize(2);
%         flagy = ryi < 1 | ryi > librarySize(1);
%         flag = flagx | flagy;
        nflag = ~flag;                
%         nflagx = ~flagx;                
%         nflagy = ~flagy;                
        
%         [RXI, RYI] = meshgrid(rxi, ryi);
        
        % Map subscripts to indices
        ind = ryi + librarySize(1)*(rxi-1);
%         Zi(flag) = Z(ind(flag));
        Zi(nflag) = Z(ind(nflag));
%         Zi(nflag) = Z(RXI(nflag), RYI(nflag));
        
