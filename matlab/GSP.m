function [obj, x] = GSP(x, gamma)
if nargin < 2
    gamma = 1;
end;
[N n] = size(x);
obj = zeros(N, 2);
alpha = 1./(2*gamma);
obj(:,1) = (1./(n.^alpha)) .* (sum(x.^2, 2).^alpha);
obj(:,2) = (1./(n.^alpha)) .* (sum((1-x).^2, 2).^alpha);