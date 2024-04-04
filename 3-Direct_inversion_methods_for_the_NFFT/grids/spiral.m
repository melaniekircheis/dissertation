function x = spiral(N)
%function x = spiral(N)
% Computes modified interleaved Archimedean spiral with parameter N
%
% Input:
% N ... number of nodes

% Set parameters
N = ceil(1.6*N);
inter = 3; % number of interleaved spirals
K = sqrt(N/8)-1; % number of cycles
K = K/inter;

% Computation of the nodes on bigger spiral
h = sqrt((0:N/inter-1)'/(N/inter)); 
x = 1/sqrt(2)*h.*[cos(3*pi*K*h),sin(3*pi*K*h)];

% Compute interleaves
phi = 2*pi/inter;
R = [cos(phi),-sin(phi);sin(phi),cos(phi)];
y = (R*x')';
z = (R*y')';

% Combine interleaves spirals
x = [x;y;z];

% Exclude nodes outside the torus
ind = (x(:,1)>=-0.5)&(x(:,1)<0.5)&(x(:,2)>=-0.5)&(x(:,2)<0.5);
x = x(ind,:);

end%function