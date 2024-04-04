% Code file for Figure 3.1
% Visualizes the amount of nonzero entries of the matrix B for d=1

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Set parameters
M = 16; % Number of Fourier coefficients
sigma = 1.0;    % Oversampling factor
m = ceil(sigma.*M);
ind = (mod(m,2)==1);  % if it's odd, make it even
m(ind) = m(ind)+1;
n = 2;  % Truncation parameter 
window = 'BSpline';
N = 2*M; % Number of nodes

% Initialization of vectors
smax = 4;
delta = zeros(smax,1); q = zeros(smax,1); t = zeros(smax,1);
q_tilde = zeros(smax,1); t_tilde = zeros(smax,1);
nnz_cols = zeros(smax,M); mean_nnz_cols = zeros(smax,2);
nnz_rows = zeros(N,smax);

%% Main
for s = 1:smax
    % Initialization of nodes
    switch s
        case 1
        % Equispaced nodes 
            x = (-0.5:1/N:0.5-1/N);
            text = 'Equispaced points'; X_equi = x;
        case 2
        % Jittered equispaced nodes 
            x = (-0.5:1/N:0.5-1/N) + 1/(4*N)*rand(1,N);
            text = 'Jittered points'; X_jittered = x;
        case 3
        % Chebyshev nodes
            x = 1/2*cos((2*(N-(0:N-1))-1)/(2*N)*pi);
            text = 'Chebyshev points'; X_cheb = x;
        case 4
        % Randomly spaced nodes
            x = sort(rand(1,N)-0.5);
            text = 'Random points'; X_rand = x;
    end%switch
    
    % Compute mesh norm and separation distance
    R = 200*(2*N);%10^6; % Number of testing nodes
    z = (-0.5:1/R:0.5-1/R)'; % Testing grid    
    dist = @(x,y) abs(mod(y-x+0.5,1)-0.5); % Distance function for periodic points
    delta(s) = 2*max(min(dist(z,x),[],2)); % Compute mesh norm
    d = [x(2:end),x(1)+1]-x; % Compute distance between neighboring points
    q(s) = min(d); t(s) = max(d); % Compute separation distance
    q_tilde(s) = 1;
    t_tilde(s) = 1;
    for j = 1:N-1
        d = dist(x(j),x(j+1:end));
        q_tilde(s) = min(q_tilde(s),min(d));
    end%for

    % Visualization of grids and respective mesh norm
    figure(1); subplot(2,smax/2,s); plot(x,zeros(size(x)),'o'); xlim([-0.5,0.5]);
    title({text,['$\delta=$',num2str(delta(s))]}); hold on

    % Generate matrix B
    B = sparse(N,prod(m));
    l = -m/2 : m/2-1;
    for j = 1:N
      y = x(j)-l./m';
      y = mod(y+0.5,1)-0.5;
      phi1 = phi(y,M(1),m(1),window,n);
      B = B + sparse(j,1:prod(m),phi1,N,prod(m));
    end%for

    % Visualization of sparsity
    figure(2); subplot(1,smax,s); spy(B); xticks([]); yticks([]); title(text); hold on

    % Count the nonzeros 
    nnz_cols(s,:) = sum(full(B>0)); % in each column
    nnz_rows(:,s) = sum(full(B>0),2); % in each row
end%for

% Finalize the plots
figure(1); sgtitle('Considered points and corresponding mesh norm');
figure(2); sgtitle({'Figure 3.1: Nonzero entries of the matrix $B \in R^{N\times |\mathcal I_{{M_{\sigma}}}|}$ for several',' choices of the nodes $x_j\in T^d$, $j=1,\dots,N,$ ',' with $d=1$, $M_\sigma=M=16$, $N=2M$ and $m=2$.'}); %