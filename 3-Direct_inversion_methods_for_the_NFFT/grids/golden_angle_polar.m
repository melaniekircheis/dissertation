function x = golden_angle_polar(R,T)
% Computes nodes of the 2-dimensional golden angle polar grid with parameters R and T

% Index sets
IR = -R/2:R/2-1;

% Computation of the golden angles
golden_angle = pi/((1+sqrt(5))/2);
start_angle = pi/2;
a = start_angle*ones(1,T);
a = a + (0:T-1) .* golden_angle;
a = mod(a,pi)-pi/2;
a = sort(a);

% Computation of the grid points
x = kron(IR'/R,[cos(a);sin(a)].');
x = mod(x+0.5,1)-0.5; % Wrap points onto the torus

end%function