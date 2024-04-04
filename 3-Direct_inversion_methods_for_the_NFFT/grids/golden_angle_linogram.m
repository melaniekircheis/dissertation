function x = golden_angle_linogram(R,T)
% Computes nodes of the 2-dimensional golden angle linogram grid with parameters R and T

% Index set
IR = -R/2:R/2-1;

% Set additional parameter
nu = 1/(2*R);

% Computation of the golden angles
golden_angle = pi/((1+sqrt(5))/2);
start_angle = pi/2;
a = start_angle*ones(1,T);
a = a + (0:T-1) .* golden_angle;
a = mod(a,pi)-pi/2;
a = sort(a);

% Initialize vectors
x1 = zeros(T*R/2,2);
x2 = zeros(T*R/2,2);

% Computation of the grid points
for t=1:T/2
    x1((t-1)*length(IR)+1:t*length(IR),:) = [IR./R+nu;(IR./R+nu)*tan(a(t+T/2)-pi/4)].';
    x2((t-1)*length(IR)+1:t*length(IR),:) = [(IR./R+nu)*cot(a(t)-pi/4);IR./R+nu].';
end

% Get complete grid
x = [x1;x2];
x = mod(x+0.5,1)-0.5; % Wrap points onto the torus

end%function