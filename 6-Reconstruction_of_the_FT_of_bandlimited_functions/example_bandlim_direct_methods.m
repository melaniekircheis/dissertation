% Code file for Figures 6.1 and 6.2

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software
error('NFFT software needs to be installed on your system first!') % path above needs to be adjusted and this line to be deleted
addpath('..\3-Direct_inversion_methods_for_the_NFFT\density_comp_factors')  % density compensation methods
addpath('..\3-Direct_inversion_methods_for_the_NFFT\grids')  % grids
addpath('..\3-Direct_inversion_methods_for_the_NFFT\opt_B')  % matrix optimization method
addpath('..\3-Direct_inversion_methods_for_the_NFFT')  % window functions of NFFT

%% Setup

% Switch flag for the two experiments
flag_experiment = 1;
switch flag_experiment
    case 1 % a) first experiment
        s = 2;
        a = 5;
        R = 2*2^a;
    case 2 % b) second experiment
        s = 3;
        a = 6;
        R = 40:8:96;
end%switch
v = char('polar','mpolar','linogram','spiral');

% Init of vectors
amount = zeros(length(R),1); 
t_exact = zeros(length(R),1); t_trafo_exact = zeros(length(R),1); t_trafo_exact_tilde = zeros(length(R),1);
t_dcf_neq_first = zeros(length(R),1); t_trafo_dcf_neq_first = zeros(length(R),1); t_trafo_dcf_neq_first_tilde = zeros(length(R),1);
t_dcf_sinc_square = zeros(length(R),1); t_trafo_dcf_sinc_square = zeros(length(R),1); t_trafo_dcf_sinc_square_tilde = zeros(length(R),1);
t_opt = zeros(length(R),1); t_trafo_opt = zeros(length(R),1); t_trafo_opt_tilde = zeros(length(R),1);
err_reconstr_exact = zeros(length(R),1); err_reconstr_exact_tilde = zeros(length(R),1);
err_reconstr_dcf_neq_first = zeros(length(R),1); err_reconstr_dcf_neq_first_tilde = zeros(length(R),1);
err_reconstr_dcf_sinc_square = zeros(length(R),1); err_reconstr_dcf_sinc_square_tilde = zeros(length(R),1);
err_reconstr_opt = zeros(length(R),1); err_reconstr_opt_tilde = zeros(length(R),1);
err_reconstr_equi = zeros(length(R),1); err_reconstr_equi_tilde = zeros(length(R),1);

%% Set parameters

z = 2^a;
M = [z;z]; % Maximum bandwidth
B = ceil(M(1)*3/8); % Actual bandwidth

%% Set Fourier data

% Set equispaced nodes in Fourier domain
k = -M(1)/2:M(1)/2-1;
[K1,K2] = ndgrid(k);

% Set Fourier coefficients
g = triangularPulse(-B,B,k);
Fhat = g.*g.';
fhat = Fhat(:);

% % Visualization of Fourier coefficients
% figure(1); fig = surf(K2,K1,Fhat);
% fig.EdgeColor = 'none'; view(0,90); colorbar

%% Main loop

for l = 1:length(R)

%% Set grid

Nx = R(l);
switch s
    case 1
        % Polar grid
        x = polar(Nx,2*Nx);

    case 2
        % Modified polar grid
        x = mpolar(Nx,2*Nx);

    case 3
        % Linogram grid
        x = linogram(Nx,2*Nx);
        
    case 4
        % Spiral grid
        x = spiral(2*Nx^2);
end%switch
x = unique(x,'stable','rows');
N = length(x); % Number of nodes
amount(l) = N;
% figure(2); plot(x(:,1),x(:,2),'o')

%% Setup of spatial data

% Set corresponding function evaluations
f1 = B*my_sinc(B*pi,x(:,1)).^2;
f2 = B*my_sinc(B*pi,x(:,2)).^2;
f = f1.*f2; %F = reshape(f,Nx,Nx);

%% Precomputation

% Exact weights
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

% Density compensation factors using normal equations of first kind, cf. [Sedarat et al.], [Rosenfeld]
tic;
w_dcf_neq_first = dcf_neq_first(x,M,M,'nfft');
t_dcf_neq_first(l) = toc;

% sinc^2 weights, cf. [Greengard et al.]
tic;
w_dcf_sinc_square = dcf_sinc_square(x,M(1));
t_dcf_sinc_square(l) = toc;

% Optimal sparse matrix
tic;
O = optB_overdet_2d(N,M,M,x,'Dirichlet',4);
t_opt(l) = toc;

%% Reconstruction from original data

plan = nfft(2,M,N); % create plan of class type nfft
plan.x = -x;

% Reconstruction using exact weights
tic; plan.f = w_exact.*f;
nfft_adjoint(plan);
fcheck_exact = plan.fhat;
t_trafo_exact(l) = toc;
err_reconstr_exact(l) = norm(fcheck_exact - fhat,2)/norm(fhat,2);

% Reconstruction using DCF_neq_first
tic; plan.f = w_dcf_neq_first.*f;
nfft_adjoint(plan);
fcheck_dcf_neq_first = plan.fhat;
t_trafo_dcf_neq_first(l) = toc;
err_reconstr_dcf_neq_first(l) = norm(fcheck_dcf_neq_first - fhat,2)/norm(fhat,2);

% Reconstruction using sinc^2 weights
tic; plan.f = w_dcf_sinc_square.*f;
nfft_adjoint(plan);
fcheck_dcf_sinc_square = plan.fhat;
t_trafo_dcf_sinc_square(l) = toc;
err_reconstr_dcf_sinc_square(l) = norm(fcheck_dcf_sinc_square - fhat,2)/norm(fhat,2);

% Reconstruction using optimal sparse matrix
[fcheck_opt,t_trafo_opt(l)] = nfft_adjoint_opt_2d(f,O,M,M,'Dirichlet',4);
err_reconstr_opt(l) = norm(fcheck_opt - fhat,2)/norm(fhat,2);

% Reshaping the results for plots
Fcheck_exact = reshape(fcheck_exact,M(1),M(2));
err_exact = Fcheck_exact-Fhat; % Exact weights
Fcheck_dcf_neq_first = reshape(fcheck_dcf_neq_first,M(1),M(2));
err_dcf_neq_first = Fcheck_dcf_neq_first-Fhat; % DCF_neq_first
Fcheck_dcf_sinc_square = reshape(fcheck_dcf_sinc_square,M(1),M(2));
err_dcf_sinc_square = Fcheck_dcf_sinc_square-Fhat; % sinc^2 weights
Fcheck_opt = reshape(fcheck_opt,M(1),M(2));
err_opt = Fcheck_opt-Fhat; % Optimal sparse matrix

%% Reconstruction from artificial data

% Compute artificial sampling data
k2 = [K1(:) K2(:)];
f_tilde = exp(2*pi*1i*x*k2')*fhat;

% Visualization
if flag_experiment == 1
S = 100;
y = -1/2:1/S:1/2-1/S; [Y1,Y2] = ndgrid(y); y = [Y1(:) Y2(:)];
fy = (B*my_sinc(B*pi,y(:,1)).^2).*(B*my_sinc(B*pi,y(:,2)).^2); Fy = reshape(fy,S,S);
fy_tilde = exp(2*pi*1i*y*k2')*fhat; Fy_tilde = reshape(real(fy_tilde),S,S);
figure(3); fig3 = subplot(1,2,1);
fig3 = surf(Y1,Y2,Fy); fig3.EdgeColor = 'none'; view(315,0);
title('Function $f$');
fig4 = subplot(1,2,2); set(gca, 'Colormap', autumn, 'NextPlot', 'replacechildren'); 
fig4 = surf(Y1,Y2,Fy-Fy_tilde); fig4.EdgeColor = 'none'; 
title('$f-\tilde f$');
sgtitle('Figure 6.1: The function (6.25) and its periodization (6.26)');
end%if

% Reconstruction using exact weights
tic; plan.f = w_exact.*f_tilde;
nfft_adjoint(plan);
fcheck_exact_tilde = plan.fhat;
t_trafo_exact_tilde(l) = toc;
err_reconstr_exact_tilde(l) = norm(fcheck_exact_tilde - fhat,2)/norm(fhat,2);

% Reconstruction using DCF_neq_first
tic; plan.f = w_dcf_neq_first.*f_tilde;
nfft_adjoint(plan);
fcheck_dcf_neq_first_tilde = plan.fhat;
t_trafo_dcf_neq_first_tilde(l) = toc;
err_reconstr_dcf_neq_first_tilde(l) = norm(fcheck_dcf_neq_first_tilde - fhat,2)/norm(fhat,2);

% Reconstruction using sinc^2 weights
tic; plan.f = w_dcf_sinc_square.*f_tilde;
nfft_adjoint(plan);
fcheck_dcf_sinc_square_tilde = plan.fhat;
t_trafo_dcf_sinc_square_tilde(l) = toc;
err_reconstr_dcf_sinc_square_tilde(l) = norm(fcheck_dcf_sinc_square_tilde - fhat,2)/norm(fhat,2);

% Reconstruction using optimal sparse matrix
[fcheck_opt_tilde,t_trafo_opt_tilde(l)] = nfft_adjoint_opt_2d(f_tilde,O,M,M,'Dirichlet',4);
err_reconstr_opt_tilde(l) = norm(fcheck_opt_tilde - fhat,2)/norm(fhat,2);

% Reshaping the results for plots
Fcheck_exact_tilde = reshape(fcheck_exact_tilde,M(1),M(2));
err_exact_tilde = Fcheck_exact_tilde-Fhat; % Exact weights
Fcheck_dcf_neq_first_tilde = reshape(fcheck_dcf_neq_first_tilde,M(1),M(2));
err_dcf_neq_first_tilde = Fcheck_dcf_neq_first_tilde-Fhat; % DCF_neq_first
Fcheck_dcf_sinc_square_tilde = reshape(fcheck_dcf_sinc_square_tilde,M(1),M(2));
err_dcf_sinc_square_tilde = Fcheck_dcf_sinc_square_tilde-Fhat; % sinc^2 weights
Fcheck_opt_tilde = reshape(fcheck_opt_tilde,M(1),M(2));
err_opt_tilde = Fcheck_opt_tilde-Fhat; % Optimal sparse matrix

%% Comparison to equispaced setting

% Initialize equispaced grid
y = -0.5:1/M(1):0.5-1/M(1);
[Y1,Y2] = ndgrid(y,y);
y = [Y1(:) Y2(:)];
  
% Generate sampling data
fy1 = B*my_sinc(B*pi,y(:,1)).^2;
fy2 = B*my_sinc(B*pi,y(:,2)).^2;
fy = fy1.*fy2;
fy_tilde = exp(2*pi*1i*y*k2')*fhat;

% Reconstruction from original data
fcheck_equi = 1/prod(M)*exp(-2*pi*1i*k2*y')*fy;
err_reconstr_equi(l) = norm(fcheck_equi - fhat,2)/norm(fhat,2);
Fcheck_equi = reshape(fcheck_equi,M(1),M(2)); err_equi = Fcheck_equi-Fhat; % Reshaping for plots

% Reconstruction from artificial data
fcheck_equi_tilde = 1/prod(M)*exp(-2*pi*1i*k2*y')*fy_tilde;
err_reconstr_equi_tilde(l) = norm(fcheck_equi_tilde - fhat,2)/norm(fhat,2);
Fcheck_equi_tilde = reshape(fcheck_equi_tilde,M(1),M(2)); err_equi_tilde = Fcheck_equi_tilde-Fhat; % Reshaping for plots

%% Visualization

if flag_experiment == 1
fig5 = figure(4); set(fig5,'units','normalized','outerposition',[0 0 1 1])
subplot(2,5,1)
    imshow(abs(err_equi),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
    title('equispaced setting')
subplot(2,5,2)
    imshow(abs(err_exact),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
    title('Use of (3.16)')
subplot(2,5,3)
    imshow(abs(err_dcf_neq_first),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
    title('Use of (3.45)')
subplot(2,5,4)
    imshow(abs(err_dcf_sinc_square),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
    title('Use of (6.13)')
subplot(2,5,5)
    imshow(abs(err_opt),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
    title('Algorithm 3.21')
subplot(2,5,6)
    imshow(abs(err_equi_tilde),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
subplot(2,5,7)
    imshow(abs(err_exact_tilde),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
subplot(2,5,8)
    imshow(abs(err_dcf_neq_first_tilde),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
subplot(2,5,9)
    imshow(abs(err_dcf_sinc_square_tilde),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
subplot(2,5,10)
    imshow(abs(err_opt_tilde),[]), c = colorbar; set(c,'TickLabelInterpreter','latex'); colormap autumn, colormap(flipud(colormap))
sgtitle({'Figure 6.2: Pointwise absolute error $\big|\mathbf{\tilde h} - \mathbf{\hat f}\big|$ of the reconstruction of the samples $\hat f(\mathbf k)$, $\mathbf k\in\mathcal I_{\mathbf M}$,',' of the tensorized triangular pulse function (3.95) with $M=32$ and $b=12$, using the density compensation factors',' computed by (3.16), (3.45) and (6.13) as well as the optimization approach from Algorithm 3.21 for the ','modified polar grid, cf. Figure 3.4b, of size $R=2M$, $T=2R$, using samples $f(\mathbf x_j)$, $j=1,\dots,N$, (top)',' and artificial samples $\tilde f(\mathbf x_j)$ in (6.26) (bottom).'});
end%if

fprintf([datestr(now,'HH:MM:SS'),'   Nx=',num2str(Nx),' done\n']);
end%for
   
%% Generate Table 6.1 (in latex format)

if flag_experiment == 2
fileID = fopen(['sampta_bandlim_',strtrim(v(s,:)),'.txt'],'w');
fprintf(fileID,'Table\n\n');
fprintf(fileID,'Nx & N & Use of (3.16) & Use of (3.45) & Use of (6.13) & Algorithm 3.21 \n');
format = '%g & %g & %1.4e & %1.4e & %1.4e & %1.4e \\\\ \n';
matrix = [R',amount,err_reconstr_exact,err_reconstr_dcf_neq_first,err_reconstr_dcf_sinc_square,err_reconstr_opt];
fprintf(fileID,format,transpose(matrix));
fclose(fileID);
end%if

%% Nested functions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
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
