% Code file for Figure 5.7

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('../4-Regularized_Shannon_sampling_formulas'); % cardinal_bspline.m

% addpath('..\nfft\matlab\nfft')  % NFFT software
% addpath('..\nfft\matlab\nnfft')  % NNFFT software
% error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted

%% Setup

% Set switch flag for saving results to txt-file
save_results = 0;

% Set parameters
N = 64;
lambda = 1;
L = (1+lambda)*N;
mmax = 10;

% Set test function
% f = @(x) sqrt(N)*my_sinc(N*pi,x);
% f = @(x) sqrt(3*N/4)*(my_sinc(N*pi/2,x)).^2;
f = @(x) sqrt(4*N/5)*(my_sinc(N*pi,x)+my_sinc(N*pi,(x-1))/2);

% Initialization of a fine grid for evaluation of reconstruction error
S = 1e5;
s = (-S:S)';
t = s/S;
ft = f(t); % Set function evaluations for comparison

% Initialization of vectors for reconstructions
Rm_Gauss = zeros(length(t),mmax); t_reg_Gauss_naive = zeros(mmax,1);
Rm_B = zeros(length(t),mmax); t_reg_B_naive = zeros(mmax,1);
Rm_sinh = zeros(length(t),mmax); t_reg_sinh_naive = zeros(mmax,1);
Rm_cKB = zeros(length(t),mmax); t_reg_cKB_naive = zeros(mmax,1);
t_step1 = zeros(mmax,1);
Rm_Gauss_sinctrafo = zeros(length(t),mmax); t_reg_Gauss_sinctrafo = zeros(mmax,1);
Rm_B_sinctrafo = zeros(length(t),mmax); t_reg_B_sinctrafo = zeros(mmax,1);
Rm_sinh_sinctrafo = zeros(length(t),mmax); t_reg_sinh_sinctrafo = zeros(mmax,1);
Rm_cKB_sinctrafo = zeros(length(t),mmax); t_reg_cKB_sinctrafo = zeros(mmax,1);

% Initialization of error vectors
err_reg_Gauss = zeros(mmax,1); err_reg_Gauss_sinctrafo = zeros(mmax,1);
err_reg_B = zeros(mmax,1); err_reg_B_sinctrafo = zeros(mmax,1);
err_reg_sinh = zeros(mmax,1); err_reg_sinh_sinctrafo = zeros(mmax,1);
err_reg_cKB = zeros(mmax,1); err_reg_cKB_sinctrafo = zeros(mmax,1);
err_max_Gauss_sinctrafo = zeros(mmax,1); err_max_B_sinctrafo = zeros(mmax,1);
err_max_sinh_sinctrafo = zeros(mmax,1); err_max_cKB_sinctrafo = zeros(mmax,1);

%%
% Setup
for m = 1:mmax

%% Precomputation of regularized sinc functions and their Fourier transforms

% Gaussian regularized sinc function
sigma = sqrt(m./(pi*L*(L-N))); % Set variance of Gaussian
psi_Gauss = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*exp(-x.^2/(2*sigma^2));
psihat_Gauss = @(v) (erf(sqrt(2)*pi*sigma*(v+L/2))-erf(sqrt(2)*pi*sigma*(v-L/2)))/(2*L);

% B-spline regularized sinc function
s = ceil((m+1)/2);
psi_B = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*(cardinal_bspline(s*L*x/m,2*s)./cardinal_bspline(0,2*s));
func = @(u) my_sinc(pi*m/(s*L),u).^(2*s);
psihat_B = @(v) arrayfun(@(j) m*integral(func,v(j)-L/2,v(j)+L/2)./(s*L^2*cardinal_bspline(0,2*s)),(1:length(v)).');

% sinh-type regularized sinc function
beta = m*pi*lambda./(1+lambda); % Set parameter of sinh
psi_sinh = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*(sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta));
rad = @(u) sqrt(u.^2/L^2-(lambda)^2/((2+2*lambda)^2));
func = @(u) besselj(1,2*pi*m*rad(u))./rad(u);
psihat_sinh = @(v) arrayfun(@(j) beta.*integral(func,v(j)-L/2,v(j)+L/2)./(2*L^2*sinh(beta)),(1:length(v)).');

% cKB regularized sinc function
beta = m*pi*lambda./(1+lambda); % Set parameter of cKB
psi_cKB = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*((besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1));
subst = @(u) 2*pi*m/(beta*L)*u;
func = @(u) sinh(beta*sqrt(1-u.^2))./(beta*sqrt(1-u.^2))-my_sinc(beta,u);
psihat_cKB = @(v) arrayfun(@(j) beta./((besseli(0,beta)-1)*L*pi).*integral(func,subst(1)*(v(j)-L/2),subst(1)*(v(j)+L/2)),(1:length(v)).');

%% Naive computation of the regularized Shannon sums

% Set equispaced function evaluations
T = m+L;
j = (-T:T)'; fj = f(j/L);

% Loop for computation of the error
for i3 = 1:length(t)
    % Set evaluation points
    ind = (j>=-m+L*t(i3,:)) & (j<=m+L*t(i3,:)); % Find indices of nonzeros
    x = t(i3,:)-j(ind).'/L; 

    % Evaluation of Gaussian regularization
    tic; 
    Rm_Gauss(i3,m) = psi_Gauss(x)*fj(ind); 
    t_reg_Gauss_naive(m) = t_reg_Gauss_naive(m) + toc;

    % Evaluation of B-Spline regularization
    tic; 
    Rm_B(i3,m) = psi_B(x)*fj(ind); 
    t_reg_B_naive(m) = t_reg_B_naive(m) + toc;

    % Evaluation of sinh-type regularization
    tic; 
    Rm_sinh(i3,m) = psi_sinh(x)*fj(ind); 
    t_reg_sinh_naive(m) = t_reg_sinh_naive(m) + toc;

    % Evaluation of cKB regularizations
    tic; 
    Rm_cKB(i3,m) = psi_cKB(x)*fj(ind);
    t_reg_cKB_naive(m) = t_reg_cKB_naive(m) + toc;
end%for

% Computation of reconstruction errors
err_reg_Gauss(m) = norm(Rm_Gauss(:,m)-ft,inf); 
err_reg_B(m) = norm(Rm_B(:,m)-ft,inf); 
err_reg_sinh(m) = norm(Rm_sinh(:,m)-ft,inf); 
err_reg_cKB(m) = norm(Rm_cKB(:,m)-ft,inf); 

%% Approximation by means of the fast sinc transform    

% % Precomputation of the weights
nu = 2*4*(1+T/L);
n = nu*N;
k = (0:n)';
epsilon = ones(n+1,1); epsilon([1,end])=1/2;
alpha = (2./(1-4*(0:n/2).^2)).';
% Fast computation of the weights using FFT
beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
y = [beta;flip(beta(2:end-1))];
w_cc = ifft(y);
w_cc = epsilon.*real(w_cc(1:n+1));
z_cc = 1/2*(cos(k*pi/n));

% % Step 1 -- NFFT (points a equispaced)
tic; plan_nfft = nfft(1,2*(T+1),n+1);
plan_nfft.x = -2*z_cc;
plan_nfft.fhat = [0;fj];
nfft_trafo(plan_nfft)
g = plan_nfft.f; t_step1(m) = toc;

% intermed = exp(2*pi*1i*(2*L)*z_cc*j'/L)*fj;
% err_step1 = norm(intermed-g,'Inf');

% % Step 2 -- Multiplication
% Gaussian window function
tau_psihat_Gauss = psihat_Gauss(2*L*z_cc); % precompute the coefficients
tic; tau_psihat_Gauss = g.*(2*L*w_cc.*tau_psihat_Gauss);
t_reg_Gauss_sinctrafo(m) = t_step1(m) + toc;
% B-spline window function
tau_psihat_B = psihat_B(2*L*z_cc); % precompute the coefficients
tic; tau_psihat_B = g.*(2*L*w_cc.*tau_psihat_B);
t_reg_B_sinctrafo(m) = t_step1(m) + toc;
% sinh-type window function
tau_psihat_sinh = psihat_sinh(2*L*z_cc); % precompute the coefficients
tic; tau_psihat_sinh = g.*(2*L*w_cc.*tau_psihat_sinh);
t_reg_sinh_sinctrafo(m) = t_step1(m) + toc;
% cKB window function
tau_psihat_cKB = psihat_cKB(2*L*z_cc); % precompute the coefficients
tic; tau_psihat_cKB = g.*(2*L*w_cc.*tau_psihat_cKB);
t_reg_cKB_sinctrafo(m) = t_step1(m) + toc;

% % Step 3 -- NNFFT (points b nonequispaced)
% Gaussian window function
tic; L_star = ceil((1+2*8/(2*L))*(2*L));
plan_nnfft = nnfft(1,n+2,length(t),2.5*L_star); 
plan_nnfft.x = t/2.5; plan_nnfft.v = 2*L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_Gauss];
nnfft_trafo(plan_nnfft);
Rm_Gauss_sinctrafo(:,m) = plan_nnfft.f;
t_reg_Gauss_sinctrafo(m) = t_reg_Gauss_sinctrafo(m) + toc;
% B-spline window function
tic; L_star = ceil((1+2*8/(2*L))*(2*L));
plan_nnfft = nnfft(1,n+2,length(t),2.5*L_star); 
plan_nnfft.x = t/2.5; plan_nnfft.v = 2*L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_B];
nnfft_trafo(plan_nnfft);
Rm_B_sinctrafo(:,m) = plan_nnfft.f;
t_reg_B_sinctrafo(m) = t_reg_B_sinctrafo(m) + toc;
% sinh-type window function
tic; L_star = ceil((1+2*8/(2*L))*(2*L));
plan_nnfft = nnfft(1,n+2,length(t),2.5*L_star); 
plan_nnfft.x = t/2.5; plan_nnfft.v = 2*L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_sinh];
nnfft_trafo(plan_nnfft);
Rm_sinh_sinctrafo(:,m) = plan_nnfft.f;
t_reg_sinh_sinctrafo(m) = t_reg_sinh_sinctrafo(m) + toc;
% cKB window function
tic; L_star = ceil((1+2*8/(2*L))*(2*L));
plan_nnfft = nnfft(1,n+2,length(t),2.5*L_star); 
plan_nnfft.x = t/2.5; plan_nnfft.v = 2*L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_cKB];
nnfft_trafo(plan_nnfft);
Rm_cKB_sinctrafo(:,m) = plan_nnfft.f;
t_reg_cKB_sinctrafo(m) = t_reg_cKB_sinctrafo(m) + toc;

% % Error computation
% Gaussian window function
err_reg_Gauss_sinctrafo(m) = norm(Rm_Gauss_sinctrafo(:,m)-ft,inf); 
err_max_Gauss_sinctrafo(m) = norm(Rm_Gauss(:,m)-Rm_Gauss_sinctrafo(:,m),'Inf');
% B-spline window function
err_reg_B_sinctrafo(m) = norm(Rm_B_sinctrafo(:,m)-ft,inf); 
err_max_B_sinctrafo(m) = norm(Rm_B(:,m)-Rm_B_sinctrafo(:,m),'Inf');
% sinh-type window function
err_reg_sinh_sinctrafo(m) = norm(Rm_sinh_sinctrafo(:,m)-ft,inf); 
err_max_sinh_sinctrafo(m) = norm(Rm_sinh(:,m)-Rm_sinh_sinctrafo(:,m),'Inf');
% cKB window function
err_reg_cKB_sinctrafo(m) = norm(Rm_cKB_sinctrafo(:,m)-ft,inf); 
err_max_cKB_sinctrafo(m) = norm(Rm_cKB(:,m)-Rm_cKB_sinctrafo(:,m),'Inf');

% help = exp(-2*pi*1i*(2*L)*t*z_cc')*tau_psihat_Gauss; 
% err_step3 = norm(help-Rm_Gauss_sinctrafo,'Inf');
% figure(2);plot(t,real(Rm_Gauss_sinctrafo))

end%for

%% Visualization

% Figure 5.7
figure(1); subplot(1,2,1); semilogy(1:mmax,err_max_Gauss_sinctrafo,'-o',1:mmax,err_max_B_sinctrafo,'-*',1:mmax,err_max_sinh_sinctrafo,'-diamond',1:mmax,err_max_cKB_sinctrafo,'-square');
legend('$\psi_{\mathrm{Gauss}}$','$\psi_{\mathrm{B}}$','$\psi_{\mathrm{sinh}}$','$\psi_{\mathrm{cKB}}$');
xlabel('$m$'); xlim([1,mmax]); xticks(1:mmax); title({'Maximum approximation', 'error (5.50)'});
subplot(1,2,2); semilogy(1:mmax,t_reg_Gauss_naive,'--o',1:mmax,t_reg_B_naive,'--*',1:mmax,t_reg_sinh_naive,'--diamond',1:mmax,t_reg_cKB_naive,'--square',1:mmax,t_reg_Gauss_sinctrafo,':o',1:mmax,t_reg_B_sinctrafo,':*',1:mmax,t_reg_sinh_sinctrafo,':diamond',1:mmax,t_reg_cKB_sinctrafo,':square');
legend('$\psi_{\mathrm{Gauss}}$','$\psi_{\mathrm{B}}$','$\psi_{\mathrm{sinh}}$','$\psi_{\mathrm{cKB}}$','Location','east');
xlabel('$m$'); xlim([1,mmax]); xticks(1:mmax); title({'Time needed for the', 'computation (in seconds)'});
colororder(["#801919";"#10520e";"#77AC30";"#00008B"])
sgtitle({'Figure 5.7: Comparison of the regularization $R_{\varphi,m} f$ in (4.74) in', 'spatial domain (dashed) and its approximation $\tilde R_{\varphi,m} f$ in (5.32) by means', 'of the fast $\mathrm{sinc}$ transform in Algorithm 5.10 (dotted) using the scaled', 'weights (5.31), for several $m\in\{1,\dots,10\}$ and $\varphi \in \{\varphi_{\mathrm{Gauss}},\, \varphi_{\mathrm{B}},\, \varphi_{\sinh}, \, {\varphi}_{\mathrm{cKB}}\}$,', 'in (4.60), (4.61), (4.62), and (4.63), for the bandlimited', 'function (4.151) with $M=64$, $\lambda=1$, and $L=(1+\lambda)M$.'});

% Visualization of the quality of the regularized Shannon formulas
figure(2); semilogy(1:mmax,err_reg_Gauss,'-',1:mmax,err_reg_Gauss_sinctrafo,'o',1:mmax,err_reg_B_sinctrafo,'-',1:mmax,err_reg_B,'*',1:mmax,err_reg_sinh,'-',1:mmax,err_reg_sinh_sinctrafo,'diamond',1:mmax,err_reg_cKB_sinctrafo,'-',1:mmax,err_reg_cKB_sinctrafo,'square');
legend('$\psi_{\mathrm{Gauss}}$ naive','$\psi_{\mathrm{Gauss}}$ sinctrafo','$\psi_{\mathrm{B}}$ naive','$\psi_{\mathrm{B}}$ sinctrafo','$\psi_{\mathrm{sinh}}$ naive','$\psi_{\mathrm{sinh}}$ sinctrafo','$\psi_{\mathrm{cKB}}$ naive','$\psi_{\mathrm{cKB}}$ sinctrafo');
xlabel('$m$'); xlim([1,mmax]); xticks(1:mmax); 
colororder(["#801919";"#801919";"#10520e";"#10520e";"#77AC30";"#77AC30";"#00008B";"#00008B"]);
title('Approximation error of regularized Shannon formula');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('test_sinc_trafo_spatial_reg.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'Error of sinc transform\n\n');
    fprintf(fileID,'psi_Gauss\n\n');
    matrix = [(1:mmax)',err_max_Gauss_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_B\n\n');
    matrix = [(1:mmax)',err_max_B_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_sinh\n\n');
    matrix = [(1:mmax)',err_max_sinh_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cKB\n\n');
    matrix = [(1:mmax)',err_max_cKB_sinctrafo];
    fprintf(fileID,format,matrix');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Computation time\n\n');
    fprintf(fileID,'psi_Gauss naive\n\n');
    matrix = [(1:mmax)',t_reg_Gauss_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_Gauss sinc trafo\n\n');
    matrix = [(1:mmax)',t_reg_Gauss_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_B naive\n\n');
    matrix = [(1:mmax)',t_reg_B_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_B sinc trafo\n\n');
    matrix = [(1:mmax)',t_reg_B_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_sinh naive\n\n');
    matrix = [(1:mmax)',t_reg_sinh_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_sinh sinc trafo\n\n');
    matrix = [(1:mmax)',t_reg_sinh_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cKB naive\n\n');
    matrix = [(1:mmax)',t_reg_cKB_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cKB sinc trafo\n\n');
    matrix = [(1:mmax)',t_reg_cKB_sinctrafo];
    fprintf(fileID,format,matrix');
fclose(fileID);
end%if

%% Nested functions

function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end