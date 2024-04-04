% Code file for Figures 3.5
% Visualizes the error bound of Theorem 3.10.

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software -- path needs to be adjusted and next line to be deleted
error('NFFT software needs to be installed on your system first!')
addpath('grids')  % grids

% Switch flag for saving results to txt-file
save_results = 0;

% Disable the warning that tolerance for PCG might be too small
warning('off','MATLAB:pcg:tooSmallTolerance')

%% Main
for d = 1%:3 % Dimension

% Init
a = 11-d; % Maximum degree
maxiter = 10; % Maximum number of repetitions

% Initialize error vectors
err_inf_abs = zeros(a,maxiter); err_inf_rel = zeros(a,maxiter);
err_2_abs = zeros(a,maxiter); err_2_rel = zeros(a,maxiter);

% Initialize time vectors
t_exact = zeros(a,maxiter); t_trafo_exact = zeros(a,maxiter);

% Initialize vectors for number of Fourier coefficients
amount = zeros(a,1);

% Initialize residual vector
r = zeros(a,1);

for l = 1:maxiter
    % Initialize random grid
    Nx = 2^(9-d);
    x = rand(Nx^d,d)-0.5;
    x = unique(x,'stable','rows');
    N = size(x,1); % Number of nodes
    
    for j = 1:a

    % Random Fourier coefficients in [1,10]
    z = 2^j;
    M = ones(d,1)*z;
    fhat = rand(prod(M),1)*9+1;                 
    amount(j) = prod(M); % Save amount for plotting

    % Compute evaluations of the trigonometric polynomial       
    plan = nfft(d,M,N); 
    plan.x = -x; 
    plan.fhat = fhat(:); 
    nfft_trafo(plan); 
    f = plan.f;

    %% iNFFT using density compensation approach

    % Precompute the weights
    tic;
    plan_2M = nfft(d,2*M,N); % create plan of class type nfft
    plan_2M.x = x; % set nodes in plan and perform precomputations
    % Init right side
    e0 = zeros(prod(2*M),1); e0(sum((2*M(1)).^(0:d-1)*M(1))+1) = 1;
    % Iteration using function handle
    if N<prod(2*M)
        [w_exact,~] = pcg(@(x) AAstar_mult(plan_2M,x,'notransp'),ones(N,1),1e-15,1000);
    else
        [y,~] = pcg(@(x) AAstar_mult(plan_2M,x,'transp'),e0,1e-15,1000);
        % NFFT
        plan_2M.fhat = y;
        nfft_trafo(plan_2M);
        w_exact = plan_2M.f;
    end%if
    t_exact(j,l) = toc;

    % Compute residual
    plan_2M.f = w_exact;
    nfft_adjoint(plan_2M);
    r(j,l) = max(abs(plan_2M.fhat - e0));
    
    % Reconstruction - adjoint NFFT
    tic; plan.f = w_exact.*f;
    nfft_adjoint(plan);
    fcheck_exact = plan.fhat;
    t_trafo_exact(j,l) = toc;

    % Compute errors
    err_inf_abs(j,l) = norm(fhat-fcheck_exact,Inf);                   % Absolute \ell_infty-error
    err_inf_rel(j,l) = norm(fhat-fcheck_exact,Inf)/norm(fhat,Inf);    % Relative \ell_infty-error
    err_2_abs(j,l) = norm(fhat-fcheck_exact,2);                       % Absolute \ell_2-error
    err_2_rel(j,l) = norm(fhat-fcheck_exact,2)/norm(fhat,2);          % Relative \ell_2-error

    % Display progress
    fprintf([datestr(datetime('now')),'   j=',num2str(j),', l=',num2str(l),' done\n']);

    end%for
end%for

%% Visualization

% Compute errors
max_err_inf_abs = max(err_inf_abs,[],2);
max_err_inf_rel = max(err_inf_rel,[],2);
max_err_2_abs = max(err_2_abs,[],2); 
max_err_2_rel = max(err_2_rel,[],2);
eps_max = max(r,[],2);

% Error plots
figure(1); subplot(1,3,d); 
loglog(amount,max_err_2_rel,'b-.^',amount,max_err_inf_rel,'-m*',amount,amount.*eps_max,'--ok')
xlabel('$|\mathcal I_{M}|$'); title(['$d=$ ',num2str(d)]);
legend('$e_2$','$e_{\infty}$','$\varepsilon$','Location','north west');
sgtitle({'Figure 3.5: Relative errors (3.94) of the reconstruction of the Fourier coefficients ','of a trigonometric polynomial (2.8) with given $\hat f_{k} \in [1,10]$, $k\in\mathcal I_{{M}}$,',' computed using the density compensation factors from Algorithm 3.9, ','for random grids with $N_t = 2^{9-d}$, $t=1,\dots,d$, and $M=M \cdot 1_d$, $M=2^c$ with $c=1,\dots,11-d$.'});

%% Generate tables for tikz
if (save_results == 1)
v_char = char('absolute ell_infty error','absolute ell_2 error','relative ell_infty error','relative ell_2 error','residual');
v = [max_err_inf_abs, max_err_2_abs, max_err_inf_rel, max_err_2_rel, amount.*eps_max];

fileID = fopen(['example_quality_densitycomp_',num2str(d),'d.txt'],'w');
format = '%g %1.4e \n';
for k = 1:size(v,2)
	fprintf(fileID,strtrim(v_char(k,:)));
    fprintf(fileID,'\n\n');
    matrix = [amount,v(:,k)];
	fprintf(fileID,format,transpose(matrix));
    fprintf(fileID,'\n'); 
end%for
end%if
end%for
if (save_results == 1), fclose(fileID); end%if

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