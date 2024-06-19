% Code file for Figure 5.4

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('..\nfft\matlab\nfft')  % NFFT software
addpath('..\nfft\matlab\nnfft')  % NNFFT software
error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted

%% Setup

% Set switch flag for saving results to txt-file
save_results = 0;
if (save_results == 1)
fileID = fopen('error_approx_sinc_ls.txt','w');
format = '%d %1.4e \n';
end%if

% Set vectors
maxsize = 13; % Maximum size of nonharmonic bandwidth
iter_ls_cheb = zeros(maxsize,3); iter_ls_leg = zeros(maxsize,3); iter_ls_equi = zeros(maxsize,3);
t_ls_cheb = zeros(maxsize,3); t_ls_leg = zeros(maxsize,3); t_ls_equi = zeros(maxsize,3); t_analytic = zeros(maxsize,3);
err_max_ls_cheb = zeros(maxsize,3); err_max_ls_leg = zeros(maxsize,3); err_max_ls_equi = zeros(maxsize,3); err_max_analytic = zeros(maxsize,3);

% Set parameters for iteration procedure
tol = 1e-11; % Tolerance for LSQR
maxit = 1000; % Maximum iteration number for LSQR

%% Test for several choices of the points
for switch_y = 1:3
for j = 1:maxsize
    % Set degree of sinc kernel
    N = 2^j;
    
    % Set evaluation points
    P = 2.5*(4*N); % Number of evaluation points
    switch switch_y
        case 1
            y = (-1:2/P:1)'; equi = 1; text = 'Equispaced'; % Equispaced points
        case 2
            y = legpts(P+1); equi = 0; text = 'Legendre'; % Legendre points
        case 3
            y = (cos((0:P)*pi/P)).'; equi = 0; text = 'Chebyshev'; % Chebyshev points
    end%switch

    % Evaluation of the sinc kernel at the evaluation points
    ev_y = sin(N*pi*y)./(N*pi*y); ev_y(y==0)=1; 

    %% Computation of weights

    % Set quadrature points
    n = 4*N; % Number of quadrature points
    z_cheb = 1/2*(cos((0:n)*pi/n)).'; % Chebyshev points
    z_leg = legpts(n+1,[-1/2,1/2]); % Legendre points
    z_equi = (-1/2:1/n:1/2)'; % Equispaced points
    
    % Least squares weights
    switch equi
        case 1 % NFFT can be used for iteration
            % Chebyshev points
            tic; plan = nfft(1,P+2,n+1); const = 2*N/P; plan.x = -z_cheb*const; % NFFT init
            [w_ls_cheb,~,~,iter_ls_cheb(j,switch_y)] = lsqr(@(x,flag) nfft_adj_mult(plan,x,flag,1),ev_y,tol,maxit); 
            t_ls_cheb(j,switch_y) = toc; 
            % Legendre points
            tic; plan = nfft(1,P+2,n+1); const = 2*N/P; plan.x = -z_leg*const;% NFFT init
            [w_ls_leg,~,~,iter_ls_leg(j,switch_y)] = lsqr(@(x,flag) nfft_adj_mult(plan,x,flag,1),ev_y,tol,maxit); 
            t_ls_leg(j,switch_y) = toc; 
            % Equispaced points
            tic; plan = nfft(1,P+2,n+1); const = 2*N/P; plan.x = -z_equi*const;% NFFT init
            [w_ls_equi,~,~,iter_ls_equi(j,switch_y)] = lsqr(@(x,flag) nfft_adj_mult(plan,x,flag,1),ev_y,tol,maxit); 
            t_ls_equi(j,switch_y) = toc; 
        case 0 % NNFFT has to be used for iteration
            % Set NNFFT parameter
            a = 1+2*8/N; N_star= ceil(a*N);
            % Chebyshev points
            tic; plan_nnfft = nnfft(1,n+2,P+1,2*N_star); % create plan of class type nfft
            plan_nnfft.x = y/2; plan_nnfft.v = N/N_star*[0;z_cheb]; % set nodes in plan
            plan_adjoint = nnfft(1,P+2,n+1,2*N_star); % create plan of class type nfft
            plan_adjoint.x = -z_cheb; plan_adjoint.v = N/N_star*[0;y]/2; % set nodes in plan
            nnfft_precompute_psi(plan_nnfft); nnfft_precompute_psi(plan_adjoint); % precomputations
            [w_ls_cheb,~,~,iter_ls_cheb(j,switch_y)] = lsqr(@(x,flag) nnfft_mult(plan_nnfft,plan_adjoint,x,flag,1,true),ev_y,tol,maxit);
            t_ls_cheb(j,switch_y) = toc; 
            % Legendre points
            tic; plan_nnfft = nnfft(1,n+2,P+1,2*N_star); % create plan of class type nfft
            plan_nnfft.x = y/2; plan_nnfft.v = N/N_star*[0;z_leg]; % set nodes in plan
            plan_adjoint = nnfft(1,P+2,n+1,2*N_star);  % create plan of class type nfft
            plan_adjoint.x = -z_leg; plan_adjoint.v = N/N_star*[0;y]/2; % set nodes in plan
            nnfft_precompute_psi(plan_nnfft); nnfft_precompute_psi(plan_adjoint); % precomputations
            [w_ls_leg,~,~,iter_ls_leg(j,switch_y)] = lsqr(@(x,flag) nnfft_mult(plan_nnfft,plan_adjoint,x,flag,1,true),ev_y,tol,maxit);
            t_ls_leg(j,switch_y) = toc; 
            % Equispaced points
            tic; plan = nfft(1,n+2,P+1); plan.x = N/n*y;
            [w_ls_equi,~,~,iter_ls_equi(j,switch_y)] = lsqr(@(x,flag) nfft_mult(plan,x,flag,1),ev_y,tol,maxit);
            t_ls_equi(j,switch_y) = toc;
    end%switch
    
    % Analytic Clenshaw-Curtis weights
    tic; k = (0:n)';
    eps = ones(n+1,1); eps([1,end])=1/2;
    alpha = (2./(1-4*(0:n/2).^2)).';
    beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
    w_analytic = [beta;flip(beta(2:end-1))];
    w_analytic = ifft(w_analytic);
    w_analytic = eps.*real(w_analytic(1:n+1));
    t_analytic(j,switch_y) = toc;

    %% Computation of function approximations

    % Compute evaluations of the sinc function
    R = 10^6; % Number of grid points for the comparison
    x = (-1:2/R:1-1/R)'; % Grid points for comparison
    ev_x = sin(N*pi*x)./(N*pi*x); ev_x(x==0)=1; 

    % Fast evaluation by means of an adjoint NFFT
    % Least squares weights obtained using Chebyshev points
    plan = nfft(1,R,n+1); 
    const = 2*N/R;
    plan.x = -z_cheb*const;
    plan.f = w_ls_cheb; 
    nfft_adjoint(plan); 
    approx_ls_cheb = plan.fhat; 
    err_max_ls_cheb(j,switch_y) = max(abs(ev_x-approx_ls_cheb));
    % Least squares weights obtained using Legendre points
    plan.x = -z_leg*const;
    plan.f = w_ls_leg; 
    nfft_adjoint(plan); 
    approx_ls_leg = plan.fhat; 
    err_max_ls_leg(j,switch_y) = max(abs(ev_x-approx_ls_leg));
    % Least squares weights obtained using Equispaced points
    plan.x = -z_equi*const;
    plan.f = w_ls_equi; 
    nfft_adjoint(plan); 
    approx_ls_equi = plan.fhat; 
    err_max_ls_equi(j,switch_y) = max(abs(ev_x-approx_ls_equi));
    % Analytic Clenshaw-Curtis weights
    plan = nfft(1,R,n+1); 
    const = 2*N/R;
    plan.x = -z_cheb*const;
    plan.f = w_analytic; 
    nfft_adjoint(plan); 
    approx_analytic = plan.fhat; 
    err_max_analytic(j,switch_y) = max(abs(ev_x-approx_analytic));
end%for

%% Visualization 

figure(1); subplot(1,3,switch_y); loglog(2.^(3:maxsize),err_max_ls_equi(3:end,switch_y),'-square',2.^(3:maxsize),err_max_ls_cheb(3:end,switch_y),'-o',2.^(3:maxsize),err_max_ls_leg(3:end,switch_y),'-*',2.^(3:maxsize),err_max_analytic(3:end,switch_y),'-diamond');
legend('$w^{\mathrm{equi}}$','$w^{\mathrm{cheb}}$','$w^{\mathrm{leg}}$','$w^{\mathrm{cc}}$');
xlabel('$M$'); title([text,' points $y_p$']);
colororder(["#000000";"#29CAF1";"#781296";"#EC801B"]); 
sgtitle({'Figure 5.4: Discrete maximum error (5.47) of the approximation of', 'the $\mathrm{sinc}$ function $\mathrm{sinc}(M\pi x)$, $x\in[-1,1]$, with $R=10^6$ for certain', 'bandwidths $M=2^s$, $s\in\{3,\dots,13\}$, shown for several choices', 'of $y_p\in[-1,1]$, $p\in \mathcal I_P$, with $P=2.5n$ and $n=4M$,', 'where the weights $w^{\mathrm{equi}}$, $w^{\mathrm{cheb}}$ and $w^{\mathrm{leg}}$ are computed using the least squares', 'approach (5.7) for $z_j\in [-\frac 12, \frac 12]$, $j=0,\dots,n$, in (5.48), (5.8), and (5.13),', 'compared to the analytic weights $w^{\mathrm{cc}}$ computed by Algorithm 5.5.'});

%% Generate tables for tikz

if (save_results == 1)
% Save approimation errors
fprintf(fileID,['Maximum approximation error for ',text,' points $y_p$\n\n']);
fprintf(fileID,'LS weights using Chebyshev points\n\n');
matrix = [2.^(1:maxsize)',err_max_ls_cheb(:,switch_y)];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n');
fprintf(fileID,'LS weights using Legendre points\n\n');
matrix = [2.^(1:maxsize)',err_max_ls_leg(:,switch_y)];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n');
fprintf(fileID,'LS weights using equispaced points\n\n');
matrix = [2.^(1:maxsize)',err_max_ls_equi(:,switch_y)];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n');
fprintf(fileID,'Analytic weights\n\n');
matrix = [2.^(1:maxsize)',err_max_analytic(:,switch_y)];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'-------------------------------------------------\n\n');
end%if
end%for

if (save_results == 1), fclose(fileID); end%if

%% Nested functions

function x = nfft_mult(plan,x,transp_flag,truncation)
      % Input check
      if nargin==2, transp_flag=[]; elseif nargin<2, error('Wrong number of input arguments.'); end
      if nargin<4 || isempty(truncation)
          truncation = 0;
      end%if

      % NFFT
      if not(strcmp(transp_flag,'transp'))
        plan.fhat = [zeros(truncation,1);x];
        nfft_trafo(plan);
        x = plan.f;
      else
        plan.f = x;
        nfft_adjoint(plan);
        x = plan.fhat;
        if truncation>=1
            x = x(truncation+1:end);
        end%if
      end
end

function x = nfft_adj_mult(plan,x,transp_flag,truncation,w)
  if nargin<4 || isempty(truncation)
      truncation = 0;
  end%if

  if strcmp(transp_flag,'transp')
    plan.fhat = [zeros(truncation,1);x];
    nfft_trafo(plan);
    if nargin < 5 
        x = plan.f;
    else
        x = w.*plan.f;
    end%if
  else
    if nargin < 5 
        plan.f = x;
    else
        plan.f = w.*x;
    end%if
    nfft_adjoint(plan);
    x = plan.fhat;
    if truncation>=1
        x = x(truncation+1:end);
    end%if
  end
end

function f = nnfft_mult(plan_nnfft,plan_adjoint,fhat,transp_flag,truncation,split_flag)
% Function handle for NNFFT to use in LSQR
% 
% Input:
%
% plans ... Initialized and precomputed plans for nnfft and adjoint
%     plan_nnfft ... plan for nnfft
%     plan_adjoint ... plan for the adjoint
%                  (simply M and N as well as v and -x are swapped)
% fhat ... vector
% transp_flag ... 'transp' computes adjoint product, else NNFFT

      % Input check
      if nargin==3, transp_flag=[]; elseif nargin<3, error('Wrong number of input arguments.'); end
      if nargin<5 || isempty(truncation)
        truncation = zeros(2,1);
      elseif sum(size(truncation))>3
          error('Invalid truncation input. Length must not be larger than 2.')
      elseif length(truncation)==1
          truncation = ones(2,1)*truncation;
      end%if
      if nargin<6 || isempty(split_flag)
        split_flag = false;
      end%if

      % NNFFT
      if not(strcmp(transp_flag,'transp'))
        plan_nnfft.fhat = [zeros(truncation(2),1);fhat];
        nnfft_trafo(plan_nnfft);
        f = plan_nnfft.f;
      else
        plan_adjoint.fhat = [zeros(truncation(1),1);fhat];
        nnfft_trafo(plan_adjoint);
        f = plan_adjoint.f;
        if split_flag
            f(1) = conj(f(1));
        end%if
      end
end