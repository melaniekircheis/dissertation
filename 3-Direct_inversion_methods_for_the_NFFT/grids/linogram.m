function [x,w] = linogram(R,T)
% Computes nodes of the 2-dimensional linogram grid with parameters R and T
% Additionally, also the weights associated to these nodes are computed

% Index sets
IR = -R/2:R/2-1;
IT = ceil(-T/4):floor((T-1)/4);

% Initialize vectors
x1 = zeros(T*R/2,2);
x2 = zeros(T*R/2,2);

% Computation of grid sets
for t=1:length(IT)
    x1((t-1)*length(IR)+1:t*length(IR),:) = [IR./R;4*IT(t)/T.*IR./R].';
    x2((t-1)*length(IR)+1:t*length(IR),:) = [-4*IT(t)/T.*IR./R;IR./R].';
end

% Get complete grid
x = [x1;x2];
x = mod(x+0.5,1)-0.5; % Wrap points onto the torus

% Computation of the weights
ind = (IR==0);
W = zeros(1,length(IR));
W(ind) = 1/(T*(1+R^2));
W(~ind) = 4*abs(IR(~ind))/(T*(1+R^2));
w = kron(ones(T,1),W');

end%function