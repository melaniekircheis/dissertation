function [w,t] = dcf_sinc_square(x,M)
%function [w,t] = dcf_sinc_square(x,M)
% Function for computing sinc^2 weights, cf. [Greengard 06]
% 
% INPUT:
% x - Vector of nonequispaced nodes in spatial domain
% 
% OUTPUT:
% w - Computed weights
% t - Computation time
% 
% Author: Melanie Kircheis


tic

% Input check - column vectors needed
if size(x,1) < size(x,2)
    x = x.';
end%if

% Computation of the weights
w = zeros(length(x),1);
for j=1:size(x,1)
    z = 1;
    for h=1:size(x,2)
        z = z .* (M*my_sinc(M*pi,(x(j,h)'-x(:,h))).^2);
    end%for
    w(j) = 1./sum(z);
end%for
t = toc;

end%function

%% Nested functions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end
