function [obj, x] = Dent(x, lambda)
if nargin < 2
    lambda = 0.85;
end;
[N n] = size(x);
obj = zeros(N, 2);

d = lambda * exp(-(x(:,1) - x(:,2)).^2);

obj(:,1) = 0.5 * (sqrt(1 + (x(:,1) + x(:,2)).^2) + sqrt(1 + (x(:,1) - x(:,2)).^2) + x(:,1) - x(:,2)) + d;
obj(:,2) = 0.5 * (sqrt(1 + (x(:,1) + x(:,2)).^2) + sqrt(1 + (x(:,1) - x(:,2)).^2) - x(:,1) + x(:,2)) + d;