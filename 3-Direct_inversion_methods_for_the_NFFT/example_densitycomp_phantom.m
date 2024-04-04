% Code file for Figure 3.8
% Testfile for the reconstruction of the 2-dimensional Shepp-Logan phantom

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software
addpath('..\nfft\applications\polarFFT')  % Shepp Logan phantom
error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted
addpath('density_comp_factors')  % density compensation methods
addpath('grids')  % grids

% Switch flag for saving results to txt-file
save_results = 0;

% Disable the warning that tolerance for PCG might be too small
warning('off','MATLAB:pcg:tooSmallTolerance')

%% Setup

% Set parameters
z = 2^6;
M = [z;z];
sigma = [1.0;1.0];    % Oversampling factor
m = ceil(sigma.*M);
ind = (mod(m,2)==1);  % if it's odd, make it even
m(ind) = m(ind)+1;

% Set grid
R = ceil(z);
s = 4;
switch s
    case 1
        % Polar grid
        x = polar(R,2*R);

    case 2
        % Modified polar grid
        x = mpolar(R,2*R);

    case 3
        % Linogram grid
        x = linogram(R,2*R);
        
    case 4
        % Spiral grid
        x = spiral(2*R^2);
end%switch
x = unique(x,'stable','rows');
N = length(x); % Number of nodes

% Visualization of the grid
figure(1); plot(x(:,1),x(:,2),'o'); title('Considered grid points $x_j$');

%% Setup of the Shepp-Logan phantom

% Initialization of the phantom
P = phantom(M(1));
fhat = P(:);

% Computation of corresponding values f     
plan = nfft(2,M,N); % create plan of class type nfft
plan.x = -x; % set nodes in plan and perform precomputations
plan.fhat = fhat(:); % set Fourier coefficients
nfft_trafo(plan); % compute nonequispaced Fourier transform
f = plan.f;

%% Computation of weights

% Exact weights
tic;
plan_2M = nfft(2,2*M,N); % create plan of class type nfft
plan_2M.x = x; % set nodes in plan and perform precomputations
e0 = zeros(prod(2*M),1); e0(prod(2*M)/2+M(1)+1) = 1;
[w_exact,~] = lsqr(@(x,transp_flag) nfft_transp_mult(plan_2M,x,transp_flag),e0,1e-15,10000);
t_exact = toc;

% Least squares approx
tic;
[w_ls,~] = pcg(@(x) AAstar_mult(plan_2M,x,'notransp'),ones(N,1),1e-15,1000);
t_ls = toc;

% Density compensation factors using pseudoinverse, cf. [Sedarat et al.]
tic;
w_dcf_pinv = dcf_pinv(x,M,M,'nfft');
t_dcf_pinv = toc;
w_dcf_pinv(w_dcf_pinv==inf) = 0;

% Density compensation factors using normal equations of first kind, cf. [Sedarat et al.], [Rosenfeld]
tic;
w_dcf_neq_first = dcf_neq_first(x,M,M,'nfft');
t_dcf_neq_first = toc;

% Density compensation factors using normal equations of second kind, cf. [Pipe, Menon]
tic;
w_dcf_neq_second = dcf_neq_second(x,M,M,'nfft');
t_dcf_neq_second = toc;
figure(2); 

% Visualization of the weights
figure(2); subplot(2,3,2); plot(real(w_exact)); title('(3.16)');
subplot(2,3,3); plot(real(w_ls)); title('(3.17)');
subplot(2,3,4); plot(real(w_dcf_pinv)); title('(3.37)');
subplot(2,3,5); plot(real(w_dcf_neq_first)); title('(3.45)');
subplot(2,3,6); plot(real(w_dcf_neq_second)); title('(3.50)');
sgtitle('Density compensation factors'); 

%% Reconstruction of Shepp-Logan phantom

% Reconstruction using exact weights
tic; plan.f = w_exact.*f;
nfft_adjoint(plan);
fcheck_exact = plan.fhat;
t_trafo_exact = toc;
err_reconstr_exact = norm(fcheck_exact - fhat,2)/norm(fhat,2);

% Reconstruction using least squares approx weights
tic; plan.f = w_ls.*f;
nfft_adjoint(plan);
fcheck_ls = plan.fhat;
t_trafo_ls = toc;
err_reconstr_ls = norm(fcheck_ls - fhat,2)/norm(fhat,2);

% Reconstruction using dcf_pinv
tic; plan.f = w_dcf_pinv.*f;
nfft_adjoint(plan);
fcheck_dcf_pinv = plan.fhat;
t_trafo_dcf_pinv = toc;
err_reconstr_dcf_pinv = norm(fcheck_dcf_pinv - fhat,2)/norm(fhat,2);

% Reconstruction using dcf_neq_first
tic; plan.f = w_dcf_neq_first.*f;
nfft_adjoint(plan);
fcheck_dcf_neq_first = plan.fhat;
t_trafo_dcf_neq_first = toc;
err_reconstr_dcf_neq_first = norm(fcheck_dcf_neq_first - fhat,2)/norm(fhat,2);

% Reconstruction using dcf_neq_second
tic; plan.f = w_dcf_neq_second.*f;
nfft_adjoint(plan);
fcheck_dcf_neq_second = plan.fhat;
t_trafo_dcf_neq_second = toc;
err_reconstr_dcf_neq_second = norm(fcheck_dcf_neq_second - fhat,2)/norm(fhat,2);

%% Visualization

% Reshape reconstructions
Pcheck_exact = reshape(fcheck_exact,M(1),M(2)); % Exact weights
Pcheck_ls = reshape(fcheck_ls,M(1),M(2)); % Least squares approx
Pcheck_dcf_pinv = reshape(fcheck_dcf_pinv,M(1),M(2)); % dcf_pinv
Pcheck_dcf_neq_first = reshape(fcheck_dcf_neq_first,M(1),M(2)); % dcf_neq_first
Pcheck_dcf_neq_second = reshape(fcheck_dcf_neq_second,M(1),M(2)); % dcf_neq_second

% Visualization of phantom and reconstructions
fig1 = figure(3); set(fig1,'units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
    imshow(P,[])
    title('Original phantom')
subplot(2,3,2)
    imshow(real(Pcheck_exact),[])
    title('Use of (3.16)')
subplot(2,3,3)
    imshow(real(Pcheck_ls),[])
    title('Use of (3.17)')
subplot(2,3,4)
    imshow(real(Pcheck_dcf_pinv),[])
    title('Use of (3.37)')
subplot(2,3,5)
    imshow(real(Pcheck_dcf_neq_first),[])
    title('Use of (3.45)')
subplot(2,3,6)
    imshow(real(Pcheck_dcf_neq_second),[])
    title('Use of (3.50)')
sgtitle({'Figure 3.8 (top): Reconstruction of the Shepp-Logan phantom of size $M= 64$ ','via Algorithm 3.2 using density compensation factors computed',' by (3.16), (3.17), (3.45) and (3.50) for the spiral grid (3.93) of size $R=M$, $T=2R$.'});

% Visualization of single row
row = floor(13/16*z);
figure(4)
subplot(2,3,1)
    plot(P(row,:))
    title('Original phantom')
    xlim([1,z])
    ylim([0,1])
subplot(2,3,2)
    plot(real(Pcheck_exact(row,:)))
    title('Use of (3.16)')
    xlim([1,z])
    ylim([0,1])
subplot(2,3,3)
    plot(real(Pcheck_ls(row,:)))
    title('Use of (3.17)')
    xlim([1,z])
    ylim([0,1])
subplot(2,3,4)
    plot(real(Pcheck_dcf_pinv(row,:)))
    title('Use of (3.37)')
    xlim([1,z])
    ylim([0,1])
subplot(2,3,5)
    plot(real(Pcheck_dcf_neq_first(row,:)))
    title('Use of (3.45)')
    xlim([1,z])
    ylim([0,1])
subplot(2,3,6)
    plot(real(Pcheck_dcf_neq_second(row,:)))
    title('Use of (3.50)')
    xlim([1,z])
    ylim([0,1])
sgtitle({'Figure 3.8 (bottom): Detailed view of the 52nd row.'});
   
%% Generate table

if (save_results == 1)
v = char('polar','mpolar','linogram','spiral');
fileID = fopen(['example_dcf_phantom',strtrim(v(s,:)),'.txt'],'w');
fprintf(fileID,'\n\n');
p = [P(row,:).',real(Pcheck_exact(row,:)).',real(Pcheck_ls(row,:)).',real(Pcheck_dcf_neq_first(row,:)).',real(Pcheck_dcf_neq_second(row,:)).'];
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

% Multiplication by A* using NFFT software
function x = nfft_transp_mult(plan,x,transp_flag)
  if strcmp(transp_flag,'transp')
    plan.fhat = x;
    nfft_trafo(plan);
    x = plan.f;
  else
    plan.f = x;
    nfft_adjoint(plan);
    x = plan.fhat;
  end
end