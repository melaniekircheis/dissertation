% Code file for Figure 5.6

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('../4-Regularized_Shannon_sampling_formulas'); % cardinal_bspline.m

addpath('..\nfft\matlab\nfft')  % NFFT software
addpath('..\nfft\matlab\nnfft')  % NNFFT software
error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted

%% Setup

% Set switch flag for saving results to txt-file
save_results = 0;

% Set parameters
N = 64;
lambda = 1;
L = (1+lambda)*N;
cmax = 10;

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
Rm_lin = zeros(length(t),cmax); t_reg_lin_naive = zeros(cmax,1);
Rm_cub = zeros(length(t),cmax); t_reg_cub_naive = zeros(cmax,1);
Rm_cos = zeros(length(t),cmax); t_reg_cos_naive = zeros(cmax,1);
Rm_conv2 = zeros(length(t),cmax); t_reg_conv2_naive = zeros(cmax,1);
t_step1 = zeros(cmax,1);
Rm_lin_sinctrafo = zeros(length(t),cmax); t_reg_lin_sinctrafo = zeros(cmax,1);
Rm_cub_sinctrafo = zeros(length(t),cmax); t_reg_cub_sinctrafo = zeros(cmax,1);
Rm_cos_sinctrafo = zeros(length(t),cmax); t_reg_cos_sinctrafo = zeros(cmax,1);
Rm_conv2_sinctrafo = zeros(length(t),cmax); t_reg_conv2_sinctrafo = zeros(cmax,1);

% Initialization of error vectors
err_reg_lin = zeros(cmax,1); err_reg_lin_sinctrafo = zeros(cmax,1);
err_reg_cub = zeros(cmax,1); err_reg_cub_sinctrafo = zeros(cmax,1);
err_reg_cos = zeros(cmax,1); err_reg_cos_sinctrafo = zeros(cmax,1);
err_reg_conv2 = zeros(cmax,1); err_reg_conv2_sinctrafo = zeros(cmax,1);
err_max_lin_sinctrafo = zeros(cmax,1); err_max_cub_sinctrafo = zeros(cmax,1);
err_max_cos_sinctrafo = zeros(cmax,1); err_max_conv2_sinctrafo = zeros(cmax,1);

%% Set frequency window functions and their inverse Fourier transforms

% Linear frequency window function
psi_lin = @(x) (N+L)/2.*my_sinc((N+L)/2*pi,x).*my_sinc((L-N)/2*pi,x);
psihat_lin = @(v) (abs(v)<=N/2) + (N/2<abs(v) & abs(v)<=L/2).*(1-(2*abs(v)-N)/(L-N));

% Cubic frequency window function
psi_cub = @(x) (abs(x)<=eps).*((L+N)/2) + (abs(x)>eps).*fillmissing((L+N)/2.*my_sinc((L+N)/2*pi,x).*(12./((L-N)^2*pi^2*x.^2).*(my_sinc((L-N)/2*pi,x)-cos((L-N)/2*pi*x))),'constant',(L+N)/2);
psihat_cub = @(v) (abs(v)<=N/2) + (N/2<abs(v) & abs(v)<=L/2).*(48/(3*(L-N)^3)*(abs(v)-L/2).^2.*(abs(v)-(3*N-L)/4));

% Raised cosine frequency window function
psi_cos = @(x) fillmissing((abs(abs(x)-1/(L-N))>eps).*((L+N)/2.*my_sinc((L+N)/2*pi,x)./(1-x.^2*(L-N)^2).*cos((L-N)/2*pi*x)),'constant',(L-N)/4*cos(N*pi/(L-N)));
psihat_cos = @(v) (abs(v)<=N/2) + (N/2<abs(v) & abs(v)<=L/2).*(1/2.*(1+cos(pi*(2*abs(v)-N)/(L-N))));

% Convolutional frequency window function
psi_conv2 = @(x) (N+L)/2.*my_sinc((N+L)/2*pi,x).*(my_sinc((L-N)/(2*2)*pi,x)).^2;
% psihat_conv = @(v) (abs(v)<=N/2) + (N/2<abs(v) & abs(v)<=L/2).*cardinal_bspline((2+1)/(L-N)*(abs(v)-N/2),2+1)./cardinal_bspline(0,2+1);
func = @(w) cardinal_bspline((2*2)/(L-N)*w,2);

%% Loop for several truncation parameters T
for c = 1:cmax

T = c*L;

% Set grid points for the comparison
R = 10^6;
x = (-1:2/R:1-1/R)'*(1*T/L);
maxnu = 2*4*(1+T/L); % Maximum size for scaling parameter

% Init of vectors
err_psi_lin = zeros(maxnu,1);
err_psi_cub = zeros(maxnu,1);
err_psi_cos = zeros(maxnu,1);
err_psi_conv2 = zeros(maxnu,1);

%% Approximation of the function psi

% Error constant
const_cc = 48/35*(1+lambda)*N.*(exp(N*(3/4*pi*(1+lambda)-(1:maxnu)*log(2)))+exp(-N*(3/4*pi*(1+lambda)+(1:maxnu)*log(2))))/2; %% = 48/35*2.^(-n).*(1+lambda)*N.*cosh(3/4*pi*(1+lambda)*N);

% Set evaluations
ev_psi_lin = 1/L*psi_lin(x);
ev_psi_cub = 1/L*psi_cub(x);
ev_psi_cos = 1/L*psi_cos(x);
ev_psi_conv2 = 1/L*psi_conv2(x);

% Compute approximation
for nu = 1:maxnu
    n = nu*N;
    k = (0:n)';
    epsilon = ones(n+1,1); epsilon([1,end])=1/2;
    alpha = (2./(1-4*(0:n/2).^2)).';
    % Fast compuation of the weights using FFT
    beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
    y = [beta;flip(beta(2:end-1))];
    w_cc = ifft(y);
    w_cc = epsilon.*real(w_cc(1:n+1));
    z_cc = 1/2*(cos(k*pi/n));

    % Fast evaluation by means of an adjoint NFFT
    plan = nfft(1,R,n+1); 
    plan.x = -z_cc*2*T/R;
    % Linear frequency window function
    plan.f = w_cc.*psihat_lin(L*z_cc); 
    nfft_adjoint(plan); 
    approx_psi_lin = plan.fhat; 
    err_psi_lin(nu) = max(abs(ev_psi_lin-approx_psi_lin));
    % Cubic frequency window function
    plan.f = w_cc.*psihat_cub(L*z_cc); 
    nfft_adjoint(plan); 
    approx_psi_cub = plan.fhat; 
    err_psi_cub(nu) = max(abs(ev_psi_cub-approx_psi_cub));
    % Raised cosine frequency window function
    plan.f = w_cc.*psihat_cos(L*z_cc); 
    nfft_adjoint(plan); 
    approx_psi_cos = plan.fhat; 
    err_psi_cos(nu) = max(abs(ev_psi_cos-approx_psi_cos));
    % Convolutional frequency window function
    psihat_conv2 = zeros(n+1,1);
    for j = 1:length(z_cc)
        psihat_conv2(j) = integral(func,L*z_cc(j)-(L+N)/4,L*z_cc(j)+(L+N)/4).*(2*2)/(L-N);
    end%for
    plan.f = w_cc.*psihat_conv2; 
    nfft_adjoint(plan); 
    approx_psi_conv2 = plan.fhat; 
    err_psi_conv2(nu) = max(abs(ev_psi_conv2-approx_psi_conv2));
end%for

% Visualization
figure(c); semilogy(1:maxnu,err_psi_lin,1:maxnu,err_psi_cub,1:maxnu,err_psi_cos,1:maxnu,err_psi_conv2,4*(1+T/L)*[1,1],[1e-10,1],'-.');
legend('$\psi_{\mathrm{lin}}$','$\psi_{\mathrm{cub}}$','$\psi_{\mathrm{cos}}$','$\psi_{\mathrm{conv,2}}$','$4(1+\frac TL)$');
xlabel('$c$'); xlim([1,maxnu]);
colororder(["#EDB120";"#4DBEEE";"#0072BD";"#D95319";"#000000"]);
title({'Error of the approximation of the regularized $\mathrm{sinc}$ functions', ['by means of exponential sums for $T=$ ',num2str(T)]});

%% Naive computation of the regularized Shannon sums

% Set equispaced function evaluations
j = (-T:T)'; fj = f(j/L);

% Loop for computation of the error
for i3 = 1:length(t)
    % Set evaluation points
    x = t(i3,:)-j.'/L; 

    % Evaluation of linear frequency regularization
    tic; 
    Rm_lin(i3,c) = psi_lin(x)/L*fj; 
    t_reg_lin_naive(c) = t_reg_lin_naive(c) + toc;

    % Evaluation of cubic frequency regularization
    tic; 
    Rm_cub(i3,c) = psi_cub(x)/L*fj; 
    t_reg_cub_naive(c) = t_reg_cub_naive(c) + toc;

    % Evaluation of raised cosine regularization
    tic; 
    Rm_cos(i3,c) = psi_cos(x)/L*fj; 
    t_reg_cos_naive(c) = t_reg_cos_naive(c) + toc;

    % Evaluation of convolution based regularizations
    tic; 
    Rm_conv2(i3,c) = psi_conv2(x)/L*fj; 
    t_reg_conv2_naive(c) = t_reg_conv2_naive(c) + toc;
end%for

% Computation of reconstruction errors
err_reg_lin(c) = norm(Rm_lin(:,c)-ft,inf); 
err_reg_cub(c) = norm(Rm_cub(:,c)-ft,inf); 
err_reg_cos(c) = norm(Rm_cos(:,c)-ft,inf); 
err_reg_conv2(c) = norm(Rm_conv2(:,c)-ft,inf); 

%% Approximation by means of the fast sinc transform    

% % Precomputation of the weights
nu = 4*(1+T/L);
n = nu*N;
k = (0:n)';
epsilon = ones(n+1,1); epsilon([1,end])=1/2;
alpha = (2./(1-4*(0:n/2).^2)).';
% Fast compuation of the weights using FFT
beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
y = [beta;flip(beta(2:end-1))];
w_cc = ifft(y);
w_cc = epsilon.*real(w_cc(1:n+1));
z_cc = 1/2*(cos(k*pi/n));

% % Step 1 -- NFFT (points a equispaced)
tic; plan_nfft = nfft(1,2*(T+1),n+1);
plan_nfft.x = -z_cc;
plan_nfft.fhat = [0;fj];
nfft_trafo(plan_nfft)
g = plan_nfft.f; t_step1(c) = toc;

% % Step 2 -- Multiplication
% Linear frequency window function
tau_psihat_lin = psihat_lin(L*z_cc); % precompute the coefficients
tic; tau_psihat_lin = g.*(w_cc.*tau_psihat_lin);
t_reg_lin_sinctrafo(c) = t_step1(c) + toc;
% Cubic frequency window function
tau_psihat_cub = psihat_cub(L*z_cc); % precompute the coefficients
tic; tau_psihat_cub = g.*(w_cc.*tau_psihat_cub);
t_reg_cub_sinctrafo(c) = t_step1(c) + toc;
% Raised cosine frequency window function
tau_psihat_cos = psihat_cos(L*z_cc); % precompute the coefficients
tic; tau_psihat_cos = g.*(w_cc.*tau_psihat_cos);
t_reg_cos_sinctrafo(c) = t_step1(c) + toc;
% Convolutional frequency window function
tau_psihat_conv2 = zeros(n+1,1);
for ind = 1:length(z_cc)
    tau_psihat_conv2(ind) = integral(func,L*z_cc(ind)-(L+N)/4,L*z_cc(ind)+(L+N)/4).*(2*2)/(L-N);
end%for
tic; tau_psihat_conv2 = g.*(w_cc.*tau_psihat_conv2);
t_reg_conv2_sinctrafo(c) = t_step1(c) + toc;

% % Step 3 -- NNFFT (points b nonequispaced)
% Linear frequency window function
tic; L_star = ceil((1+2*8/L)*L);
plan_nnfft = nnfft(1,n+2,length(t),3*L_star); 
plan_nnfft.x = t/3; plan_nnfft.v = L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_lin];
nnfft_trafo(plan_nnfft);
Rm_lin_sinctrafo(:,c) = plan_nnfft.f;
t_reg_lin_sinctrafo(c) = t_reg_lin_sinctrafo(c) + toc;
% Cubic frequency window function
tic; L_star = ceil((1+2*8/L)*L);
plan_nnfft = nnfft(1,n+2,length(t),3*L_star); 
plan_nnfft.x = t/3; plan_nnfft.v = L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_cub];
nnfft_trafo(plan_nnfft);
Rm_cub_sinctrafo(:,c) = plan_nnfft.f;
t_reg_cub_sinctrafo(c) = t_reg_cub_sinctrafo(c) + toc;
% Raised cosine frequency window function
tic; L_star = ceil((1+2*8/L)*L);
plan_nnfft = nnfft(1,n+2,length(t),3*L_star); 
plan_nnfft.x = t/3; plan_nnfft.v = L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_cos];
nnfft_trafo(plan_nnfft);
Rm_cos_sinctrafo(:,c) = plan_nnfft.f;
t_reg_cos_sinctrafo(c) = t_reg_cos_sinctrafo(c) + toc;
% Convolutional frequency window function
tic; L_star = ceil((1+2*8/L)*L);
plan_nnfft = nnfft(1,n+2,length(t),3*L_star); 
plan_nnfft.x = t/3; plan_nnfft.v = L/L_star*[0;z_cc]; % set nodes in plan (scaling necessary!)
nnfft_precompute_psi(plan_nnfft); % precomputations
plan_nnfft.fhat = [0;tau_psihat_conv2];
nnfft_trafo(plan_nnfft);
Rm_conv2_sinctrafo(:,c) = plan_nnfft.f;
t_reg_conv2_sinctrafo(c) = t_reg_conv2_sinctrafo(c) + toc;

% % Error computation
% Linear frequency window function
err_reg_lin_sinctrafo(c) = norm(Rm_lin_sinctrafo(:,c)-ft,inf); 
err_max_lin_sinctrafo(c) = norm(Rm_lin(:,c)-Rm_lin_sinctrafo(:,c),'Inf');
% Cubic frequency window function
err_reg_cub_sinctrafo(c) = norm(Rm_cub_sinctrafo(:,c)-ft,inf); 
err_max_cub_sinctrafo(c) = norm(Rm_cub(:,c)-Rm_cub_sinctrafo(:,c),'Inf');
% Raised cosine frequency window function
err_reg_cos_sinctrafo(c) = norm(Rm_cos_sinctrafo(:,c)-ft,inf); 
err_max_cos_sinctrafo(c) = norm(Rm_cos(:,c)-Rm_cos_sinctrafo(:,c),'Inf');
% Convolutional frequency window function
err_reg_conv2_sinctrafo(c) = norm(Rm_conv2_sinctrafo(:,c)-ft,inf); 
err_max_conv2_sinctrafo(c) = norm(Rm_conv2(:,c)-Rm_conv2_sinctrafo(:,c),'Inf');

% intermed = exp(2*pi*1i*z_cc*j')*fj;
% err_step1 = norm(intermed-g,'Inf');
% help = exp(-2*pi*1i*L*t*z_cc')*tau_psihat_lin; 
% err_step3 = norm(help-Rm_lin_sinctrafo,'Inf');

end%for

%% Visualization
 
% Figure 5.6
figure(cmax+1); subplot(1,2,1); semilogy((1:cmax)*L,err_max_lin_sinctrafo,'-o',(1:cmax)*L,err_max_cub_sinctrafo,'-square',(1:cmax)*L,err_max_cos_sinctrafo,'-diamond',(1:cmax)*L,err_max_conv2_sinctrafo,'-^');
legend('$\psi_{\mathrm{lin}}$','$\psi_{\mathrm{cub}}$','$\psi_{\mathrm{cos}}$','$\psi_{\mathrm{conv,2}}$');
xlabel('$T$'); xlim([L,cmax*L]); xticks((2:2:cmax)*L); xticklabels({'$2L$','$4L$','$6L$','$8L$','$10L$'});
colororder(["#EDB120";"#4DBEEE";"#0072BD";"#D95319"]);
title({'Maximum approximation', 'error (5.49)'});
subplot(1,2,2); semilogy((1:cmax)*L,t_reg_lin_naive,'--o',(1:cmax)*L,t_reg_cub_naive,'--square',(1:cmax)*L,t_reg_cos_naive,'--diamond',(1:cmax)*L,t_reg_conv2_naive,'--^',(1:cmax)*L,t_reg_lin_sinctrafo,':o',(1:cmax)*L,t_reg_cub_sinctrafo,':square',(1:cmax)*L,t_reg_cos_sinctrafo,':diamond',(1:cmax)*L,t_reg_conv2_sinctrafo,':^');
legend('$\psi_{\mathrm{lin}}$','$\psi_{\mathrm{cub}}$','$\psi_{\mathrm{cos}}$','$\psi_{\mathrm{conv,2}}$','Location','east');
xlabel('$T$'); xlim([L,cmax*L]); xticks((2:2:cmax)*L); xticklabels({'$2L$','$4L$','$6L$','$8L$','$10L$'});
title({'Time needed for the', 'computation (in seconds)'});
sgtitle({'Figure 5.6: Comparison of the regularization $P_{\psi,T} f$ in (4.38) in', 'frequency domain (dashed) and its approximation $\tilde P_{\psi,T} f$ in (5.30)', 'by means of the fast $\mathrm{sinc}$ transform in Algorithm 5.10 (dotted)', 'using the scaled weights (5.29), for several $T=cL$ with $c\in\{1,\dots,10\}$ and', '$\psi \in \{\psi_{\mathrm{lin}},\, \psi_{\mathrm{cub}},\, \psi_{\cos}, \, {\psi}_{\mathrm{conv},2}\}$, in (4.45), (4.48), (4.51), and (4.57),', 'for the bandlimited function (4.151) with $M=64$, $\lambda=1$, and $L=(1+\lambda)M$.'});

% Visualization of the quality of the regularized Shannon formulas
figure(cmax+2); semilogy((1:cmax)*L,err_reg_lin,'-',(1:cmax)*L,err_reg_lin_sinctrafo,'o',(1:cmax)*L,err_reg_cub,'-',(1:cmax)*L,err_reg_cub_sinctrafo,'square',(1:cmax)*L,err_reg_cos,'-',(1:cmax)*L,err_reg_cos_sinctrafo,'diamond',(1:cmax)*L,err_reg_conv2,'-',(1:cmax)*L,err_reg_conv2_sinctrafo,'^');
legend('$\psi_{\mathrm{lin}}$ naive','$\psi_{\mathrm{lin}}$ sinctrafo','$\psi_{\mathrm{cub}}$ naive','$\psi_{\mathrm{cub}}$ sinctrafo','$\psi_{\mathrm{cos}}$ naive','$\psi_{\mathrm{cos}}$ sinctrafo','$\psi_{\mathrm{conv,2}}$ naive','$\psi_{\mathrm{conv,2}}$ sinctrafo');
xlabel('$T$'); xlim([L,cmax*L]); xticks((2:2:cmax)*L); xticklabels({'$2L$','$4L$','$6L$','$8L$','$10L$'});
colororder(["#EDB120";"#EDB120";"#4DBEEE";"#4DBEEE";"#0072BD";"#0072BD";"#D95319";"#D95319"]);
title('Approximation error of regularized Shannon formula');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('test_sinc_trafo_freq_reg.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'Error of sinc transform\n\n');
    fprintf(fileID,'psi_lin\n\n');
    matrix = [(1:cmax)'*L,err_max_lin_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cub\n\n');
    matrix = [(1:cmax)'*L,err_max_cub_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cos\n\n');
    matrix = [(1:cmax)'*L,err_max_cos_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_conv2\n\n');
    matrix = [(1:cmax)'*L,err_max_conv2_sinctrafo];
    fprintf(fileID,format,matrix');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Computation time\n\n');
    fprintf(fileID,'psi_lin naive\n\n');
    matrix = [(1:cmax)'*L,t_reg_lin_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_lin sinc trafo\n\n');
    matrix = [(1:cmax)'*L,t_reg_lin_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cub naive\n\n');
    matrix = [(1:cmax)'*L,t_reg_cub_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cub sinc trafo\n\n');
    matrix = [(1:cmax)'*L,t_reg_cub_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cos naive\n\n');
    matrix = [(1:cmax)'*L,t_reg_cos_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_cos sinc trafo\n\n');
    matrix = [(1:cmax)'*L,t_reg_cos_sinctrafo];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_conv2 naive\n\n');
    matrix = [(1:cmax)'*L,t_reg_conv2_naive];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\npsi_conv2 sinc trafo\n\n');
    matrix = [(1:cmax)'*L,t_reg_conv2_sinctrafo];
    fprintf(fileID,format,matrix');
fclose(fileID);
end%if

%% Nested functions

function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end