function [x,w] = mpolar(R,T)
% Computes nodes of the 2-dimensional modified polar grid with parameters R and T
% Additionally, also the weights associated to these nodes are computed

% Index sets
IR = ceil(-sqrt(2)*R/2):floor(sqrt(2)*R/2);
IT = ceil(-T/2):floor((T-1)/2);

% Computation of the nodes on bigger circles
x = kron(IR'/R,[cos(pi*IT/T);sin(pi*IT/T)].');

% Exclude nodes outside the torus
ind = (x(:,1)>=-0.5)&(x(:,1)<0.5)&(x(:,2)>=-0.5)&(x(:,2)<0.5);
x = x(ind,:);

% Computation of the weights
ind2 = (IR==0);
W = zeros(1,length(IR));
W(ind2) = 1/4;
W(~ind2) = abs(IR(~ind2));
W = kron(W',ones(T,1));
w = W(ind);
w = w./sum(w);

end%function