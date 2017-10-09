function yi = nearestNeighborInterpolation1D(x, y, xi)

% --- NOTE: Input sample spacing must be uniform

% --- Forces vectors to be columns
x = x(:); xi = xi(:); y = y(:);

% --- Gets the x spacing. One over to perform divide only once
ndx = 1 / (x(2) - x(1)); % 

% --- Subtract minimum of x from xi
xi = xi - x(1);      

% --- Fills output with NaNs
yi = NaN * ones(size(xi));

% --- Indices of nearest-neighbors
rxi = round(xi * ndx) + 1;        

% --- Finds out of bounds indices
flag = rxi < 1 | rxi > length(x);
                                      
% --- Finds in-bound indices
nflag = ~flag;                
        
yi(nflag) = y(rxi(nflag));

