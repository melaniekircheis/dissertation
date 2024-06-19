% Code file for the verification of the accuracy of the tensor decomposition

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('..\nfft\matlab\nfft')  % NFFT software
error('NFFT software needs to be installed on your system first!') % path above needs to be adjusted and this line to be deleted

%% Two-dimensional setting -- Clenshaw-Curtis

% Set evaluation points
R = 6; % Number of grid points
x_1d = (-1:2/R:1-1/R)'; % 1d grid points
[X1,X2] = ndgrid(x_1d);
x_2d = [X2(:) X1(:)]; % 2d grid points
x_2d_kron = [kron(x_1d,ones(R,1)),kron(ones(R,1),x_1d)];
verification_x_2d_kron = norm(x_2d - x_2d_kron);

% Set evaluations of sinc functions
N = 24;
tic; ev_x_1d = sin(N*pi*x_1d)./(N*pi*x_1d); ev_x_1d(x_1d==0)=1; t_ev_1d = toc;
tic; ev_x_2d = sin(N*pi*x_2d)./(N*pi*x_2d); ev_x_2d(x_2d==0)=1; ev_x_2d = prod(ev_x_2d,2); t_ev_2d = toc;
tic; ev_x_2d_kron = kron(ev_x_1d,ev_x_1d); t_ev_2d_kron = t_ev_1d + toc;
verification_ev_x_2d_kron = norm(ev_x_2d - ev_x_2d_kron);

% Compute Clenshaw-Curtis weights
n = 4*N;
k = (0:n)';
eps = ones(n+1,1); eps([1,end])=1/2;
alpha = (2./(1-4*(0:n/2).^2)).';
% Fast computation of the weights using FFT
beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
y = [beta;flip(beta(2:end-1))];
w_cc_1d = ifft(y);
w_cc_1d = eps.*real(w_cc_1d(1:n+1));
w_cc_2d = kron(w_cc_1d,w_cc_1d);

% Set Chebychev points
z_cc_1d = 1/2*(cos(k*pi/n));
[Z1,Z2] = ndgrid(z_cc_1d);
z_cc_2d = [Z2(:) Z1(:)];
z_cc_2d_kron = [kron(z_cc_1d,ones(n+1,1)),kron(ones(n+1,1),z_cc_1d)];
verification_z_cc_kron_2d = norm(z_cc_2d - z_cc_2d_kron);
fprintf(['Error in the tensorization of 2d Chebyshev points = ',num2str(verification_z_cc_kron_2d),'\n']);

% Compute the approximation by means of an adjoint NFFT
tic; plan = nfft(2,[R,R],(n+1)^2); 
plan.x = -z_cc_2d*2*N/R;
plan.f = w_cc_2d; 
nfft_adjoint(plan); 
approx_x_cc_2d = plan.fhat; t_approx_cc_2d = toc;
% Verification using 1d approximation
tic; plan = nfft(1,R,n+1); 
plan.x = -z_cc_1d*2*N/R;
plan.f = w_cc_1d; 
nfft_adjoint(plan); 
approx_x_cc_1d = plan.fhat; t_approx_cc_1d = toc;
tic; approx_x_cc_2d_kron = kron(approx_x_cc_1d,approx_x_cc_1d); t_approx_cc_2d_kron = t_approx_cc_1d + toc;
verification_approx_cc_kron_2d = norm(approx_x_cc_2d - approx_x_cc_2d_kron);
fprintf(['Error in the tensorization of 2d approximation = ',num2str(verification_approx_cc_kron_2d),'\n']);

% Compute approximation error
err_max_cc_2d = max(abs(ev_x_2d-approx_x_cc_2d));
err_max_cc_1d = max(abs(ev_x_1d-approx_x_cc_1d));
err_cc_2d = zeros(R,1);
for r = 1:R
    err_cc_2d(r) = max(abs(ev_x_1d(r)*ev_x_1d-approx_x_cc_1d(r)*approx_x_cc_1d));
end%for
err_max_cc_2d_kron = max(err_cc_2d);
verification_err_max_cc_kron_2d = norm(err_max_cc_2d - err_max_cc_2d_kron);

%% Three-dimensional setting -- Clenshaw-Curtis

% Set evaluation points
[X1,X2,X3] = ndgrid(x_1d);
x_3d = [X3(:) X2(:) X1(:)];
x_3d_kron = [kron(x_1d,ones(R^2,1)),kron(ones(R,1),[kron(x_1d,ones(R,1)),kron(ones(R,1),x_1d)])];
verification_x_3d_kron = norm(x_3d - x_3d_kron);

% Set evaluations of sinc functions
ev_x_3d = sin(N*pi*x_3d)./(N*pi*x_3d); ev_x_3d(x_3d==0)=1; ev_x_3d = prod(ev_x_3d,2);
tic; ev_x_3d_kron = kron(ev_x_1d,kron(ev_x_1d,ev_x_1d)); t_ev_3d_kron = t_ev_1d + toc;
verification_ev_x_3d_kron = norm(ev_x_3d - ev_x_3d_kron);

% Compute Clenshaw-Curtis weights
w_cc_3d = kron(w_cc_1d,kron(w_cc_1d,w_cc_1d));

% Set Chebychev points
[Z1,Z2,Z3] = ndgrid(z_cc_1d);
z_cc_3d = [Z3(:) Z2(:) Z1(:)];
z_cc_3d_kron = [kron(z_cc_1d,ones((n+1)^2,1)),kron(ones(n+1,1),[kron(z_cc_1d,ones(n+1,1)),kron(ones(n+1,1),z_cc_1d)])];
verification_z_cc_kron_3d = norm(z_cc_3d - z_cc_3d_kron);
fprintf(['Error in the tensorization of 3d Chebyshev points = ',num2str(verification_z_cc_kron_3d),'\n']);

% Compute the approximation by means of an adjoint NFFT
tic; plan = nfft(3,[R,R,R],(n+1)^3); 
plan.x = -z_cc_3d*2*N/R;
plan.f = w_cc_3d; 
nfft_adjoint(plan); 
approx_x_cc_3d = plan.fhat; t_approx_cc_3d = toc;
% Verification using 1d approximation
tic; approx_x_cc_3d_kron = kron(approx_x_cc_1d,kron(approx_x_cc_1d,approx_x_cc_1d)); t_approx_cc_3d_kron = t_approx_cc_1d + toc;
verification_approx_cc_kron_3d = norm(approx_x_cc_3d - approx_x_cc_3d_kron);
fprintf(['Error in the tensorization of 3d approximation = ',num2str(verification_approx_cc_kron_3d),'\n']);

% Compute approximation error
err_max_cc_3d = max(abs(ev_x_3d-approx_x_cc_3d));
err_cc_3d = zeros(R,1);
for r = 1:R
    h = zeros(R,1);
    for p = 1:R
        h(p) = max(abs(ev_x_1d(p)*ev_x_1d(r)*ev_x_1d-approx_x_cc_1d(p)*approx_x_cc_1d(r)*approx_x_cc_1d));
    end%for
    err_cc_3d(r) = max(h);
end%for
err_max_cc_3d_kron = max(err_cc_3d);
verification_err_max_cc_kron_3d = norm(err_max_cc_3d - err_max_cc_3d_kron);

%% Two-dimensional setting -- Gauss-Legendre

% Compute Legendre weights
n = 2*N;
[z_leg_1d,w_leg_1d] = legpts(n+1,[-1/2,1/2]); % legendre points and weights, source: https://github.com/chebfun/chebfun
w_leg_1d = w_leg_1d';
w_leg_2d = kron(w_leg_1d,w_leg_1d);

% Set Legendre points
[Z1,Z2] = ndgrid(z_leg_1d);
z_leg_2d = [Z2(:) Z1(:)];
z_leg_2d_kron = [kron(z_leg_1d,ones(n+1,1)),kron(ones(n+1,1),z_leg_1d)];
verification_z_leg_kron_2d = norm(z_leg_2d - z_leg_2d_kron);
fprintf(['Error in the tensorization of 2d Legendre points = ',num2str(verification_z_leg_kron_2d),'\n']);

% Compute the approximation by means of an adjoint NFFT
tic; plan = nfft(2,[R,R],(n+1)^2); 
plan.x = -z_leg_2d*2*N/R;
plan.f = w_leg_2d; 
nfft_adjoint(plan); 
approx_x_leg_2d = plan.fhat; t_approx_leg_2d = toc;
% Verification using 1d approximation
tic; plan = nfft(1,R,n+1); 
plan.x = -z_leg_1d*2*N/R;
plan.f = w_leg_1d; 
nfft_adjoint(plan); 
approx_x_leg_1d = plan.fhat; t_approx_leg_1d = toc;
tic; approx_x_leg_2d_kron = kron(approx_x_leg_1d,approx_x_leg_1d); t_approx_leg_2d_kron = t_approx_leg_1d + toc;
verification_approx_leg_kron_2d = norm(approx_x_leg_2d - approx_x_leg_2d_kron);
fprintf(['Error in the tensorization of 2d approximation = ',num2str(verification_approx_leg_kron_2d),'\n']);

% Compute approximation error
err_max_leg_2d = max(abs(ev_x_2d-approx_x_leg_2d));
err_max_leg_1d = max(abs(ev_x_1d-approx_x_leg_1d));
err_leg_2d = zeros(R,1);
for r = 1:R
    err_leg_2d(r) = max(abs(ev_x_1d(r)*ev_x_1d-approx_x_leg_1d(r)*approx_x_leg_1d));
end%for
err_max_leg_2d_kron = max(err_leg_2d);
verification_err_max_leg_kron_2d = norm(err_max_leg_2d - err_max_leg_2d_kron);

%% Three-dimensional setting -- Gauss-Legendre

% Compute Legendre weights
w_leg_3d = kron(w_leg_1d,kron(w_leg_1d,w_leg_1d));

% Set Legendre points
[Z1,Z2,Z3] = ndgrid(z_leg_1d);
z_leg_3d = [Z3(:) Z2(:) Z1(:)];
z_leg_3d_kron = [kron(z_leg_1d,ones((n+1)^2,1)),kron(ones(n+1,1),[kron(z_leg_1d,ones(n+1,1)),kron(ones(n+1,1),z_leg_1d)])];
verification_z_leg_kron_3d = norm(z_leg_3d - z_leg_3d_kron);
fprintf(['Error in the tensorization of 3d Legendre points = ',num2str(verification_z_leg_kron_3d),'\n']);

% Compute the approximation by means of an adjoint NFFT
tic; plan = nfft(3,[R,R,R],(n+1)^3); 
plan.x = -z_leg_3d*2*N/R;
plan.f = w_leg_3d; 
nfft_adjoint(plan); 
approx_x_leg_3d = plan.fhat; t_approx_leg_3d = toc;
% Verification using 1d approximation
tic; approx_x_leg_3d_kron = kron(approx_x_leg_1d,kron(approx_x_leg_1d,approx_x_leg_1d)); t_approx_leg_3d_kron = t_approx_leg_1d + toc;
verification_approx_leg_kron_3d = norm(approx_x_leg_3d - approx_x_leg_3d_kron);
fprintf(['Error in the tensorization of 3d approximation = ',num2str(verification_approx_leg_kron_3d),'\n']);

% Compute approximation error
err_max_leg_3d = max(abs(ev_x_3d-approx_x_leg_3d));
err_leg_3d = zeros(R,1);
for r = 1:R
    h = zeros(R,1);
    for p = 1:R
        h(p) = max(abs(ev_x_1d(p)*ev_x_1d(r)*ev_x_1d-approx_x_leg_1d(p)*approx_x_leg_1d(r)*approx_x_leg_1d));
    end%for
    err_leg_3d(r) = max(h);
end%for
err_max_leg_3d_kron = max(err_leg_3d);
verification_err_max_leg_kron_3d = norm(err_max_leg_3d - err_max_leg_3d_kron);