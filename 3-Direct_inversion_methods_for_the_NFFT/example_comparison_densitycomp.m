% Code file for Figure 3.6 and 3.7
% Compares several density compensation factors

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software -- path needs to be adjusted and next line to be deleted
error('NFFT software needs to be installed on your system first!')
addpath('density_comp_factors')  % density compensation methods

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;
if (save_results == 1), fileID = fopen(['example_comparison_densitycomp_',strtrim(v(s,:)),'.txt'],'w'); end%if

v = char('jittered','Chebyshev','random');
s = 1;
a = 1:10;

% Initialization of vectors
t_voronoi = zeros(length(a),1); t_trafo_voronoi = zeros(length(a),1);
t_exact = zeros(length(a),1); t_trafo_exact = zeros(length(a),1);
t_ls = zeros(length(a),1); t_trafo_ls = zeros(length(a),1);
t_dcf_pinv = zeros(length(a),1); t_trafo_dcf_pinv = zeros(length(a),1);
t_dcf_neq_first = zeros(length(a),1); t_trafo_dcf_neq_first = zeros(length(a),1);
t_dcf_neq_second = zeros(length(a),1); t_trafo_dcf_neq_second = zeros(length(a),1);
err_reconstr_voronoi = zeros(length(a),1);
err_reconstr_exact = zeros(length(a),1);
err_reconstr_ls = zeros(length(a),1);
err_reconstr_dcf_pinv = zeros(length(a),1);
err_reconstr_dcf_neq_first = zeros(length(a),1);
err_reconstr_dcf_neq_second = zeros(length(a),1);

%% Main
for l = 1:length(a)

% Set parameters
M = 2^a(l);
N = ceil(2*M); % Number of nodes

% Set grid
switch s
    case 1  % Jittered equispaced nodes 
        x = (-0.5+1/(2*N):1/N:0.5-1/(2*N)) + 1/(4*N)*(2*rand(1,N)-1);
    case 2  % Chebyshev nodes
        x = 1/2*cos((2*(N-(0:N-1))-1)/(2*N)*pi);
    case 3  % Randomly spaced nodes
        x = sort(rand(1,N)-0.5);
end
if l==3; figure(1); plot(x,zeros(N,1),'o'); xlabel('$x_j$'); title('Given nodes'); end

%% Problem setup 

% Initialization of the Fourier coefficients
k = (-M/2:M/2-1)'; fhat = 1-abs(k./ceil(M*3/8)); fhat(abs(k)>ceil(M*3/8)) = 0;

% Computation of corresponding values f     
plan = nfft(1,M,N); % create plan of class type nfft
plan.x = -x'; % set nodes in plan and perform precomputations
plan.fhat = fhat(:); % set Fourier coefficients
nfft_trafo(plan); % compute nonequispaced Fourier transform
f = plan.f;

%% Computation of weights

% Voronoi weights
tic;
h = abs(x-[x(2:end),x(1)+1]);
w_voronoi = 1/2*(h+[h(end),h(1:end-1)]).';
t_voronoi(l) = toc;

% Exact weights
tic;
plan_2M = nfft(1,2*M,N); % create plan of class type nfft
plan_2M.x = x'; % set nodes in plan and perform precomputations
e0 = zeros(2*M,1); e0(M+1) = 1;
[w_exact,~] = lsqr(@(x,transp_flag) nfft_transp_mult(plan_2M,x,transp_flag),e0,1e-15,10000);
t_exact(l) = toc;

% Least squares approx
tic;
[w_ls,~] = pcg(@(x) AAstar_mult(plan_2M,x,'notransp'),ones(N,1),1e-15,1000);
t_ls(l) = toc;

% Density compensation factors using pseudoinverse, cf. [Sedarat et al.]
tic;
w_dcf_pinv = dcf_pinv(x,M,M,'nfft');
t_dcf_pinv(l) = toc;
w_dcf_pinv(w_dcf_pinv==inf) = 0;

% Density compensation factors using normal equations of first kind, cf. [Sedarat et al.], [Rosenfeld]
tic;
w_dcf_neq_first = dcf_neq_first(x,M,M,'nfft');
t_dcf_neq_first(l) = toc;

% Density compensation factors using normal equations of second kind, cf. [Pipe, Menon]
tic;
w_dcf_neq_second = dcf_neq_second(x,M,M,'nfft');
t_dcf_neq_second(l) = toc;

if l==3; figure(2); 
subplot(2,3,1); stem(x,real(w_exact)); title('Use of (3.16)');
subplot(2,3,2); stem(x,real(w_ls)); title('Use of (3.17)');
subplot(2,3,3); stem(x,w_voronoi); title('Voronoi weights');
subplot(2,3,4); stem(x,real(w_dcf_pinv)); title('Use of (3.37)');
subplot(2,3,5); stem(x,w_dcf_neq_first); title('Use of (3.45)');
subplot(2,3,6); stem(x,w_dcf_neq_second); title('Use of (3.50)');
sgtitle({'Figure 3.7: Exemplary density compensation factors $w_j$ for given jittered',' nodes $x_j$, $j=1,\dots,N$, cf. (3.87), computed using (3.16), (3.17), (3.37), (3.45) and (3.50),',' as well as Voronoi weights with $d=1$, $M=8$ and $N=2M$.'}); 

if (save_results == 1)
fprintf(fileID,'Weights\n\n');
p = [w_exact,w_ls,w_voronoi,w_dcf_pinv,w_dcf_neq_first,w_dcf_neq_second];
format = '%g %1.4e \n';
for k = 1:size(p,2)
	matrix = [x.',p(:,k)];
	fprintf(fileID,format,transpose(matrix));
    fprintf(fileID,'\n------------------------------\n\n\n');
end%for
end%if
end%if

%% Reconstruction of the Fourier coefficients

% Reconstruction using Voronoi weights
tic; plan.f = w_voronoi.*f;
nfft_adjoint(plan);
fcheck_voronoi = plan.fhat;
t_trafo_voronoi(l) = toc;
err_reconstr_voronoi(l) = norm(fcheck_voronoi - fhat,2)/norm(fhat,2);

% Reconstruction using exact weights
tic; plan.f = w_exact.*f;
nfft_adjoint(plan);
fcheck_exact = plan.fhat;
t_trafo_exact(l) = toc;
err_reconstr_exact(l) = norm(fcheck_exact - fhat,2)/norm(fhat,2);

% Reconstruction using least squares approx weights
tic; plan.f = w_ls.*f;
nfft_adjoint(plan);
fcheck_ls = plan.fhat;
t_trafo_ls(l) = toc;
err_reconstr_ls(l) = norm(fcheck_ls - fhat,2)/norm(fhat,2);

% Reconstruction using DCF_pinv
tic; plan.f = w_dcf_pinv.*f;
nfft_adjoint(plan);
fcheck_dcf_pinv = plan.fhat;
t_trafo_dcf_pinv(l) = toc;
err_reconstr_dcf_pinv(l) = norm(fcheck_dcf_pinv - fhat,2)/norm(fhat,2);

% Reconstruction using DCF_neq_first
tic; plan.f = w_dcf_neq_first.*f;
nfft_adjoint(plan);
fcheck_dcf_neq_first = plan.fhat;
t_trafo_dcf_neq_first(l) = toc;
err_reconstr_dcf_neq_first(l) = norm(fcheck_dcf_neq_first - fhat,2)/norm(fhat,2);

% Reconstruction using DCF_neq_second
tic; plan.f = w_dcf_neq_second.*f;
nfft_adjoint(plan);
fcheck_dcf_neq_second = plan.fhat;
t_trafo_dcf_neq_second(l) = toc;
err_reconstr_dcf_neq_second(l) = norm(fcheck_dcf_neq_second - fhat,2)/norm(fhat,2);
  
end%for

%% Visualization

figure(3); loglog(2.^a,err_reconstr_exact,'-*',2.^a,err_reconstr_voronoi,'square-',2.^a,err_reconstr_dcf_pinv,'^-',2.^a,err_reconstr_dcf_neq_first,'diamond-',2.^a,err_reconstr_dcf_neq_second,'-o');
xlabel('$M$'); ylabel('$e_2$');
colororder(["#35B263";"#4DBEEE";"#D95319";"#EDB120";"#0072BD"]);
legend('Use of (3.16)','Use of Voronoi weights','Use of (3.37)','Use of (3.45)','Use of (3.50)');
title({'Figure 3.6: Relative errors (3.94) of the reconstruction of the Fourier coefficients',' of a trigonometric polynomial (2.8) with given $\hat f_{k}$, $k\in\mathcal I_{{M}}$, computed by',' means of Algorithm 3.2 using Voronoi weights as well as the density compensation factors',' from (3.16), (3.37), (3.45), (3.50), for one-dimensional jittered grids (3.87)',' with $N=2M$ and $M=2^c$, $c=1,\dots,10$.'}); 
   
%% Generate tables

if (save_results == 1)
fprintf(fileID,'\n\n----------------------------------------------------------------\n\n');
fprintf(fileID,'Reconstruction errors\n\n');
p = [err_reconstr_exact,err_reconstr_ls,err_reconstr_voronoi,err_reconstr_dcf_pinv,err_reconstr_dcf_neq_first,err_reconstr_dcf_neq_second];
format = '%g %1.4e \n';
for k = 1:size(p,2)
	matrix = [(2.^a)',p(:,k)];
	fprintf(fileID,format,transpose(matrix));
    fprintf(fileID,'\n------------------------------\n\n\n');
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