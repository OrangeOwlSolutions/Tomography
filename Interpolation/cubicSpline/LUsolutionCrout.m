function x = LUsolutionCrout(L, U, c)
        
[numRows, numCols] = size(L);

y = zeros(1, numRows);

% --- Solving L y = c
y(1) = c(1) / L(1, 1);
for k = 2 : numRows
   y(k) = (c(k) - L(k, k - 1) * y(k - 1)) / L(k, k);
end

% --- Solving U x = y
x(numRows) = y(numRows);
for k = (numRows - 1) : -1 : 1
    x(k) = y(k) - U(k, k + 1) * x(k + 1);
end

