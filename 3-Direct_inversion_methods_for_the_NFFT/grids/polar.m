function [x,w] = polar(R,T)
% Computes nodes of the 2-dimensional polar grid with parameters R and T
% Additionally, also the normalized weights associated to these nodes are computed

% Index sets
IR = -R/2:R/2-1;
IT = ceil(-T/2):floor((T-1)/2);

% Computation of the nodes
x = kron(IR'/R,[cos(pi*IT/T);sin(pi*IT/T)].');
x = mod(x+0.5,1)-0.5; % Wrap points onto the torus

% Computation of the weights
ind = (IR==0);
W = zeros(1,length(IR));
W(ind) = 1/(T*(1+R^2));
W(~ind) = 4*abs(IR(~ind))/(T*(1+R^2));
w = kron(W',ones(T,1));

end%function