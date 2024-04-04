% Code file for Figure 3.11
% Compares the density compensation method and the optimization method
% for the reconstruction of the 2-dimensional Shepp-Logan phantom

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software
addpath('..\nfft\applications\polarFFT')  % Shepp Logan phantom
error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted
addpath('grids')  % grids
addpath('opt_B')  % optimization methods for sparse matrix B

% Disable the warning that tolerance for PCG might be too small
warning('off','MATLAB:pcg:tooSmallTolerance')

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Switch flag for selection of grid
s = 4;

% Set phantom sizes to consider
a = 1:10;

% Init of vectors
t_opt = zeros(size(a));
t_trafo_opt = zeros(size(a));
t_exact = zeros(size(a));
t_trafo_exact = zeros(size(a));
err_reconstr_opt = zeros(size(a));
err_reconstr_voronoi = zeros(size(a));
err_reconstr_exact = zeros(size(a));

%% Main
for l = a

% Set parameters
z = 2^l;
M = [z;z];
sigma = [1.0;1.0];    % Oversampling factor
m = ceil(sigma.*M);
ind = (mod(m,2)==1);  % if it's odd, make it even
m(ind) = m(ind)+1;
n = 4; % cut-off parameter

% Set grid
R = 2*z;
switch s
    case 1
        % 2d-jittered grid
        R = floor(sqrt(2)*R);
        x = (-0.5+1/(2*R):1/R:0.5-1/(2*R));
        [X1,X2] = ndgrid(x,x);
        x = [X1(:) X2(:)];
        x = x + 1/(4*R)*(2*rand(R^2,2)-1);
        clear X1 X2;

    case 2
        % Polar grid
        [x,w] = polar(R,2*R);

    case 3
        % Modified polar grid
        [x,w] = mpolar(R,2*R);

    case 4
        % Linogram grid
        [x,w] = linogram(R,2*R);
        
    case 5
        % Random grid
        x = rand((floor(sqrt(2)*R))^2,2)-0.5;
end%switch

help=[x,w];
help = unique(help,'stable','rows');
x = help(:,1:2);
w = help(:,3);
N = size(x,1); % Number of nodes

%% Fast computation of optimized matrix B

[O,t_opt(l)] = optB_overdet_2d(N,M,m,x,'Dirichlet',n);

%% Computation of exact weights

tic;
plan_2M = nfft(2,2*M,N); % create plan of class type nfft
plan_2M.x = -x; % set nodes in plan and perform precomputations
if N<prod(2*M)
    [w_exact,~] = pcg(@(x) AAstar_mult(plan_2M,x,'notransp'),ones(N,1),1e-15,1000);
else
    % Init right side
    e0 = zeros(prod(2*M),1); e0(2*M(1)*M(1)+M(1)+1) = 1;
    % Iteration using function handle
    [y,~] = pcg(@(x) AAstar_mult(plan_2M,x,'transp'),e0,1e-15,1000);
    % NFFT
    plan_2M.fhat = y;
    nfft_trafo(plan_2M);
    w_exact = plan_2M.f;
end%if
t_exact(l) = toc;

%% Reconstruction of Shepp-Logan phantom

% Initialization of the phantom
P = phantom(M(1));
fhat = P(:);

% Computation of corresponding values f     
plan = nfft(2,M,N); % create plan of class type nfft
plan.x = -x; % set nodes in plan and perform precomputations
plan.fhat = fhat(:); % set Fourier coefficients
nfft_trafo(plan); % compute nonequispaced Fourier transform
f = plan.f;

% Reconstruction using optimization method
[fcheck_opt,t_trafo_opt(l)] = nfft_adjoint_opt_2d(f,O,m,M,'Dirichlet',n);
err_reconstr_opt(l) = norm(fcheck_opt - fhat,2)/norm(fhat,2);

% Reconstruction using Voronoi weights
plan.f = w.*f;
nfft_adjoint(plan);
fcheck_voronoi = plan.fhat;
err_reconstr_voronoi(l) = norm(fcheck_voronoi - fhat,2)/norm(fhat,2);

% Reconstruction using exact weights
tic; plan.f = w_exact.*f;
nfft_adjoint(plan);
fcheck_exact = plan.fhat;
t_trafo_exact(l) = toc;
err_reconstr_exact(l) = norm(fcheck_exact - fhat,2)/norm(fhat,2);

%% Visualization

% Reshape the reconstructions
Pcheck_voronoi = reshape(fcheck_voronoi,M(1),M(2)); % Voronoi weights
Pcheck_exact = reshape(fcheck_exact,M(1),M(2)); % Exact weights
Pcheck = reshape(fcheck_opt,M(1),M(2)); % Optimization

% Visualization of the phantom and its reconstructions
fig1 = figure(1); set(fig1,'units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
    imshow(P)
    title('Original phantom')
subplot(1,3,2)
    imshow(real(Pcheck_exact),[])
    title('Algorithm 3.2')
subplot(1,3,3)
    imshow(real(Pcheck))
    title('Algorithm 3.20')
sgtitle({'Figure 3.11 (top): Reconstruction of the Shepp-Logan phantom of size $1024\times 1024$',' via the density compensation method from Algorithm 3.2 using weights',' computed by Algorithm 3.9 compared to the optimization method from',' Algorithm 3.20 using the modified matrix computed by Algorithm 3.24',' for the linogram grid (3.90) of size $R=M=1024$, $T=2R$.'});

% Display progress
fprintf([datestr(now,'HH:MM:SS'),'   z=',num2str(z),' done\n']);
end%for

% Visualization of a single line
row = floor(13/16*z);
figure(2)
subplot(1,3,1)
    plot(P(row,:))
    title('Original phantom')
    xlim([1,z])
    ylim([0,1])
subplot(1,3,2)
    plot(real(Pcheck_exact(row,:)))
    title('Algorithm 3.2')
    xlim([1,z])
    ylim([0,1])
subplot(1,3,3)
    plot(real(Pcheck(row,:)))
    title('Algorithm 3.20')
    xlim([1,z])
    ylim([0,1])
sgtitle({'Figure 3.11 (bottom): Detailed view of the 832nd row each.'});
   
%% Generate tables for tikz 

if (save_results == 1)
v = char('jittered','polar','mpolar','linogram','random');
fileID = fopen(['phantom_sizes_',strtrim(v(s,:)),'.txt'],'w');
fprintf(fileID,'Table\n\n');
fprintf(fileID,'phantom size & relative error (voronoi) & relative error (exact) & relative error (opt) & precompute exact & optimization time & reconstruction time (exact) & reconstruction time (opt) \n');
format = '%g & %1.4e & %1.4e & %1.4e & %1.4e & %1.4e & %1.4e & %1.4e \\\\ \n';
t = [t_exact',t_opt',t_trafo_exact',t_trafo_opt'];
matrix = [2.^a',err_reconstr_voronoi',err_reconstr_exact',err_reconstr_opt',t];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n\n');
fprintf(fileID,'------------------------------------------------\n');
fprintf(fileID,'Data for row plots');
p = [P(row,:).',real(Pcheck_voronoi(row,:)).',real(Pcheck_exact(row,:)).',real(Pcheck(row,:)).'];
format = '%g %1.4e \n';
for k = 1:size(p,2)
    fprintf(fileID,'\n\n'); 
	matrix = [(1:z)',p(:,k)];
	fprintf(fileID,format,transpose(matrix));
    fprintf(fileID,'\n');
end%for
fclose(fileID);
end%if

%% Nested functions

% Multiplication by A and A* using NFFT software
function x = AAstar_mult(plan,x,transp_flag)
    if strcmp(transp_flag,'notransp')
        plan.f = x;
        nfft_adjoint(plan);
        nfft_trafo(plan);
        x = plan.f;
    else
        plan.fhat = x;
        nfft_trafo(plan);
        nfft_adjoint(plan);
        x = plan.fhat;
    end%if
end