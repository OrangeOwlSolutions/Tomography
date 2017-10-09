function [L, U] = LUdecompositionCrout(A)
        
[numRows, numCols] = size(A);

% --- Setting first element of L
L(1, 1) = A(1, 1);

% --- Setting lower diagonal of L
for k = 2 : numRows
    L(k, k - 1) = A(k, k - 1);
end

% --- Setting the main diagonal of U
for k = 1 : numRows
    U(k, k) = 1;
end

% --- Alternately setting the upper diagonal of L and the lower diagonal of
% U
for k = 1 : numRows - 1,
   U(k, k + 1) = A(k, k + 1) / L(k, k);
   L(k + 1, k + 1) = A(k + 1, k + 1) - L(k + 1, k) * U(k, k + 1);
end
