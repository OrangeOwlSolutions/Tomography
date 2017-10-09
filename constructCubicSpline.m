function y = constructCubicSpline(xin, xout, s0, s1, s2, s3)

y = zeros(1, length(xout));

n = length(xin);

[~, xout_pos] = histc(xout, xin);

xout_pos = max(xout_pos, 1);            % --- To avoid index = 0 when xin < xout(1)
xout_pos = min(xout_pos, n - 1);        % --- To avoid index = n + 1 when xin > xout(end).

for k = 1 : length(xout)
    y(k) = s0(xout_pos(k)) + s1(xout_pos(k)) * (xout(k) - xin(xout_pos(k))) + ...
           s2(xout_pos(k)) * (xout(k) - xin(xout_pos(k))).^2 + s3(xout_pos(k)) * (xout(k) - xin(xout_pos(k))).^3;
end