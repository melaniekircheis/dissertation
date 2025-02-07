% Code file for Figure 3.10
% Compares several matrix optimization methods

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software
addpath('..\nfft\applications\polarFFT')  % Shepp Logan phantom
error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted
addpath('opt_B')  % optimization methods for sparse matrix B

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Set parameters
M = 64; % number of Fourier coefficients
sigma = 1.0;    % Oversampling factor
m = ceil(sigma.*M);
ind = (mod(m,2)==1);  % if it's odd, make it even
m(ind) = m(ind)+1;
n = 2;  % Truncation parameter 

% Initialization of nodes
N = 2*M; % number of nodes
s = 1; % switch flag
switch s
    case 1
    % Jittered equispaced nodes 
        x = (-0.5+1/(2*N):1/N:0.5-1/(2*N)) + 1/(4*N)*(2*rand(1,N)-1);
    case 2
    % Chebyshev nodes
        x = 1/2*cos((2*(N-(0:N-1))-1)/(2*N)*pi);
    case 3
    % Randomly spaced nodes
        x = sort(rand(1,N)-0.5);
end

%% Generate original matrices

% Initialize identity matrices 
IM = eye(prod(M));
Im = eye(m);

% Generate matrix A
A = zeros(N,M);
j = 1:N;
for k = -M/2 : M/2-1
    A(:,k + M/2 + 1) = exp(2*pi*1i*k*x(j));
end

% Generate matrix F
F = zeros(m,M);
l = -m/2 : m/2-1;
for k = -M/2 : M/2-1
    F(:,k + M/2 + 1) = exp(2*pi*1i*k*l/m);
end

% Generate matrix D
window = 'BSpline';
D = 1/m * sparse(1:M,1:M,1./phi_hat(-M/2:M/2-1,M,m,window,n),M,M);
D_inv = inv(D); % Compute the inverse matrix
D_Dirichlet = 1/m * sparse(1:M,1:M,1./phi_hat(-M/2:M/2-1,M,m,'Dirichlet',n),M,M);

% Generate matrix B
B = sparse(N,prod(m));
for j = 1:N
  y = x(j)-l./m';
  y = mod(y+0.5,1)-0.5;
  phi1 = phi(y,M(1),m(1),window,n);
  B = B + sparse(j,1:prod(m),phi1,N,prod(m));
end

% Compute the pseudoinverse (of adjoint matrix)
B_dagger = pinv(full(B'));

% % Visualization of sparsity
% figure(1); subplot(1,2,1); spy(B); title('$B$'); subplot(1,2,2); spy(B_dagger); title('$(B^*)^\dagger$'); sgtitle('Sparsity of pseudoinverse')

% Sanity check - compute Frobenius norms using dense pseudoinverse
accuracy_pinv = norm(B'*B_dagger-Im,'fro'); 
prod_dense = 1/prod(m)*D_inv'*F'*B_dagger'*B*F*D; norm_dense = norm(prod_dense-IM,'fro');

%% Compute simple sparsification of dense pseudoinverse

% Generate sparse approximation -- simply cut-out of dense matrix
B_dagger_sub = sparse(N,prod(m));
for j = 1:N
  ind = find(B(j,1:prod(m))>0); % Find nonzero indices of matrix B
  phi1 = B_dagger(j,1:prod(m)); % Entries of appropriate row of the dense matrix
  B_dagger_sub = B_dagger_sub + sparse(j,ind,phi1(ind),N,prod(m)); % Compose sparse approximation
end

% % Visualisation of structure of the dense pseudoinverse
% figure(2); subplot(1,4,1); spy(B); title('$B$'); subplot(1,4,2); spy(B_dagger>1e-5); title('$(B^*)^\dagger>1e-5$'); 
% subplot(1,4,3); spy(B_dagger>1e-3); title('$(B^*)^\dagger>1e-3$'); subplot(1,4,4); spy(B_dagger_sub); title('$(B^*)^\dagger_{\mathrm{sub}}$'); 
% sgtitle('Sparsity of reconstructed matrix')

% Compute Frobenius norm using the sparse pseudoinverse
prod_sparse = 1/prod(m)*D_inv'*F'*B_dagger_sub'*B*F*D; norm_sparse = norm(prod_sparse-IM,'fro');

%% Examine reconstruction properties

% Define test case
ex = 'deterministic';

% Generate test data
k = -M/2:M/2-1;
switch ex
    case 'constant'
        fhat = ones(M,1);
    case 'random'
        fhat = rand(M,1)+1;
    case 'deterministic'
%         P = phantom(M); fhat = P(13/16*M,:).';
        b = ceil(3/8*M); fhat = 1-abs(k'/b); fhat(abs(k)>=b) = 0;
end%switch
f = A*fhat;

% Compute reconstructions
fcheck_dense = 1/prod(m)*D_inv'*F'*B_dagger'*f;
fcheck_sparse = 1/prod(m)*D_inv'*F'*B_dagger_sub'*f;

% Compute reconstruction errors
switch ex
    case 'random'
        err_dense = abs(fhat-fcheck_dense)./fcheck_dense; maxerr_dense = max(err_dense);
        err_sparse = abs(fhat-fcheck_sparse)./fcheck_sparse; maxerr_sparse = max(err_sparse);
    otherwise
        err_dense = abs(fhat-fcheck_dense); maxerr_dense = max(err_dense);
        err_sparse = abs(fhat-fcheck_sparse); maxerr_sparse = max(err_sparse);
end%switch
accuracy_reconstr = [norm(err_dense,inf), norm(err_sparse,inf)];

% Visualization of reconstruction
figure(3); subplot(1,2,1); plot(k,fhat,'-ko',k,real(fcheck_dense),'-*',k,real(fcheck_sparse),'-*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('given coeffs','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$');
subplot(1,2,2); semilogy(k,zeros(M,1),k,err_dense,k,err_sparse); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$');
sgtitle('Reconstruction and error');

%% Examine optimization procedure

% Compute optimized matrix
[S_opt,time_opt,~] = optB_pinv_1d(N,M,m,x,window,n,[],'pcg');
% figure(4); subplot(1,4,1); spy(B); title('$B$'); subplot(1,4,2); spy(B_dagger); title('$(B^*)^\dagger$'); 
% subplot(1,4,3); spy(B_dagger_sub); title('$(B^*)^\dagger_{\mathrm{sub}}$'); subplot(1,4,4); spy(S_opt); title('$S$'); 
% sgtitle('Sparsity of reconstructed matrix')

% Compute Frobenius norm using the optimized pseudoinverse
prod_opt = D_inv'*F'*S_opt'*B*F*D; norm_opt = norm(prod_opt-IM,'fro');

% Compute reconstruction and reconstruction error
fcheck_opt = D_inv'*F'*S_opt'*f;
err_opt = abs(fhat-fcheck_opt); maxerr_opt = max(err_opt);

% Visualization of reconstruction
figure(3); subplot(1,2,1); plot(k,fhat,'-ko',k,real(fcheck_dense),'-*',k,real(fcheck_sparse),'-*',k,real(fcheck_opt),'-*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('given coeffs','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$','reconstr $S$');
subplot(1,2,2); semilogy(k,zeros(M,1),k,err_dense,k,err_sparse,k,err_opt); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$','reconstr $S$');
sgtitle('Reconstruction and error');

%% Comparison with previously used optimization using A'BFD

% Computation of corresponding reconstruction and error
[B_opt,~] = optB_overdet_1d(N,M,m,x,'Dirichlet',n);
fcheck_opt_overdet = D_Dirichlet'*F'*B_opt'*f;
err_opt_overdet = abs(fhat-fcheck_opt_overdet); maxerr_opt_overdet = max(err_opt_overdet);

% Visualization of reconstruction
figure(3); subplot(1,2,1); plot(k,fhat,'-ko',k,real(fcheck_dense),'-*',k,real(fcheck_sparse),'-*',k,real(fcheck_opt),'-*',k,real(fcheck_opt_overdet),'-*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('given coeffs','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$','reconstr $S$','reconstr $B_{\mathrm{opt}}$');
subplot(1,2,2); semilogy(k,zeros(M,1),k,err_dense,k,err_sparse,k,err_opt,k,err_opt_overdet); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$','reconstr $S$','reconstr $B_{\mathrm{opt}}$');
sgtitle('Reconstruction and error');

%% Examine truncated optimization procedure

% Compute optimized matrices
S_opt_half = optB_pinv_1d(N,M,m,x,window,n,m/4,'pcg');
S_opt_three = optB_pinv_1d(N,M,m,x,window,n,3,'pcg');
S_opt_two = optB_pinv_1d(N,M,m,x,window,n,2,'pcg');
S_opt_one = optB_pinv_1d(N,M,m,x,window,n,1,'pcg');
S_opt_zero = optB_pinv_1d(N,M,m,x,window,n,0,'pcg');

% Compute Frobenius norms using the optimized pseudoinverse
prod_opt_half = D_inv'*F'*S_opt_half'*B*F*D; norm_opt_half = norm(prod_opt_half-IM,'fro');
prod_opt_three = D_inv'*F'*S_opt_three'*B*F*D; norm_opt_three = norm(prod_opt_three-IM,'fro');
prod_opt_two = D_inv'*F'*S_opt_two'*B*F*D; norm_opt_two = norm(prod_opt_two-IM,'fro');
prod_opt_one = D_inv'*F'*S_opt_one'*B*F*D; norm_opt_one = norm(prod_opt_one-IM,'fro');
prod_opt_zero = D_inv'*F'*S_opt_zero'*B*F*D; norm_opt_zero = norm(prod_opt_zero-IM,'fro');

% Compute reconstructions and reconstruction errors
fcheck_opt_half = D_inv'*F'*S_opt_half'*f;
err_opt_half = abs(fhat-fcheck_opt_half); maxerr_opt_half = max(err_opt_half);
fcheck_opt_three = D_inv'*F'*S_opt_three'*f;
err_opt_three = abs(fhat-fcheck_opt_three); maxerr_opt_three = max(err_opt_three);
fcheck_opt_two = D_inv'*F'*S_opt_two'*f;
err_opt_two = abs(fhat-fcheck_opt_two); maxerr_opt_two = max(err_opt_two);
fcheck_opt_one = D_inv'*F'*S_opt_one'*f;
err_opt_one = abs(fhat-fcheck_opt_one); maxerr_opt_one = max(err_opt_one);
fcheck_opt_zero = D_inv'*F'*S_opt_zero'*f;
err_opt_zero = abs(fhat-fcheck_opt_zero); maxerr_opt_zero = max(err_opt_zero);

% Visualization of reconstruction
figure(4); subplot(1,2,1); plot(k,fhat,'-ko',k,real(fcheck_opt),'-*',k,real(fcheck_opt_half),'-*',k,real(fcheck_opt_three),'-*',k,real(fcheck_opt_two),'-*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); title('Reconstruction');
legend('given coeffs','$r = 32$','$r = 16$','$r=3$','$r=2$');
subplot(1,2,2); semilogy(k,zeros(M,1),k,err_opt,k,err_opt_half,k,err_opt_three,k,err_opt_two); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); title('Pointwise reconstruction error');
legend('','$r = 32$','$r = 16$','$r=3$','$r=2$');
sgtitle('The reconstruction $D^{-1} F^* S_{\mathrm{opt}} f$ with several choices of $r \leq \frac{|\mathcal I_{M_\sigma}|}{2} = 32$.');

%% Comparison to [GeSo14] -  sparse version only

% Suppress warning when matrices are rank decifient
warning('off','MATLAB:rankDeficientMatrix')

% Setup corresponding matrices
IN = eye(N);
T = A*D'*F'*B';
Lambda = sparse(N,N);
s = ceil(log(N)); % bandwidth parameter

% Columnwise computation of the sparse matrix
for j = 1:N
    ind = max(1,j-s+1):min(N,j+s-1);
    Tj = T(:,ind);
    ej = IN(:,j);
    lambda = Tj\ej;
    Lambda = Lambda + sparse(ind,j,lambda,N,N);
end%for

% Computation of corresponding reconstruction and error
fcheck_opt_Lambda = D'*F'*B'*Lambda*f;
err_opt_Lambda = abs(fhat-fcheck_opt_Lambda); maxerr_opt_Lambda = max(err_opt_Lambda);

% Visualization of reconstruction
figure(3); subplot(1,2,1); plot(k,fhat,'-ko',k,real(fcheck_dense),'-*',k,real(fcheck_sparse),'-*',k,real(fcheck_opt),'-*',k,real(fcheck_opt_overdet),'-*',k,real(fcheck_opt_Lambda),'-*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('given coeffs','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$','reconstr $S$','reconstr $B_{\mathrm{opt}}$','reconstr $\Lambda$');
subplot(1,2,2); semilogy(k,zeros(M,1),k,err_dense,k,err_sparse,k,err_opt,k,err_opt_overdet,k,err_opt_Lambda); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); legend('','reconstr $(B^*)^\dagger$','reconstr $(B^*)^\dagger_{\mathrm{sub}}$','reconstr $S$','reconstr $B_{\mathrm{opt}}$','reconstr $\Lambda$');
sgtitle('Reconstruction and error');

%% Comparison to [KiPo19] - optimization using FDA'B

% Computation of corresponding reconstruction and error
plan = infft(x,M,'m',n,'sigma',sigma); % Initialization and node-dependent precomputations
plan.f = f; % Set function values
infft_trafo(plan); % Compute inverse nonequispaced Fourier transform
fcheck_opt_old = plan.fcheck;
err_opt_old = abs(fhat-fcheck_opt_old); maxerr_opt_old = max(err_opt_old);

% Visualization of reconstruction
figure(3); subplot(1,2,1); plot(k,fhat,'-',k,real(fcheck_dense),'-square',k,real(fcheck_sparse),'--x',k,real(fcheck_opt),'-o',k,real(fcheck_opt_Lambda),':diamond',k,real(fcheck_opt_old),'-^',k,real(fcheck_opt_overdet),'-.*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); title('Reconstruction');
colororder(["#000000";"#0072BD";"#4DBEEE";"#35B263";"#D95319";"#EDB120";"#A2142F"]);
legend('$\hat f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \, D^{-1} F^* (B^*)^\dagger f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \, D^{-1} F^* (B^*)^\dagger_{\mathrm{sub}} f$','$D^{-1} F^* S_{\mathrm{opt}} f$','$D^* F^* B^* \Lambda f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \,I_{|\mathcal I_{M}|} F^* \tilde B_{\mathrm{opt}}^* f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \,I_{|\mathcal I_{M}|} F^* B_{\mathrm{opt}}^* f$','Location','northeastoutside');
subplot(1,2,2); semilogy(k,zeros(M,1),k,err_dense,'-square',k,err_sparse,'--x',k,err_opt,'-o',k,err_opt_Lambda,':diamond',k,err_opt_old,'-^',k,err_opt_overdet,'-.*'); xlim([min(k),max(k)]); 
xticks([-M/2,M/2-1]); xticklabels({'$-\frac M2$','$\frac M2-1$'}); title('Poinwise reconstruction error');
legend('$\hat f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \, D^{-1} F^* (B^*)^\dagger f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \, D^{-1} F^* (B^*)^\dagger_{\mathrm{sub}} f$','$D^{-1} F^* S_{\mathrm{opt}} f$','$D^* F^* B^* \Lambda f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \,I_{|\mathcal I_{M}|} F^* \tilde B_{\mathrm{opt}}^* f$','$\frac{1}{|\mathcal I_{{M_{\sigma}}}|} \,I_{|\mathcal I_{M}|} F^* B_{\mathrm{opt}}^* f$','Location','northeastoutside');
sgtitle({'Figure 3.10: Reconstructions of the triangular pulse (3.95) and ','pointwise errors for the different matrix optimization methods ','from Section 3.3 with $d=1$, $M=64$, $\sigma=1.0$, $m=2$',' and $N=2M$ jittered points (3.87).'});

%% Generate tables for tikz

if ( save_results == 1 )
fileID = fopen('example_comparison_optB.txt','w');
format = '%d %1.4e \n';
% Save reconstructions
fprintf(fileID,'Reconstructions');
fprintf(fileID,'\n\n--------------------');
fprintf(fileID,'\n\nOriginal\n\n');
matrix = [k',fhat];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nreconstr pinv\n\n');
matrix = [k',fcheck_dense];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nreconstr cutout pinv\n\n');
matrix = [k',fcheck_sparse];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nreconstr opt S\n\n');
matrix = [k',fcheck_opt];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nreconstr opt B\n\n');
matrix = [k',fcheck_opt_overdet];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nreconstr opt Lambda\n\n');
matrix = [k',fcheck_opt_Lambda];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nreconstr opt B (old)\n\n');
matrix = [k',fcheck_opt_old];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n-----------------------------------------------------\n\n');
% Save errors
fprintf(fileID,'Pointwise errors');
fprintf(fileID,'\n\n--------------------');
fprintf(fileID,'\n\nerr pinv\n\n');
matrix = [k',err_dense];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nerr cutout pinv\n\n');
matrix = [k',err_sparse];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nerr opt S\n\n');
matrix = [k',err_opt];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nerr opt B\n\n');
matrix = [k',err_opt_overdet];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nerr opt Lambda\n\n');
matrix = [k',err_opt_Lambda];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\nerr opt B (old)\n\n');
matrix = [k',err_opt_old];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n-----------------------------------------------------\n\n');
% Save errors for subproblems
fprintf(fileID,'Maximum errors for solving subproblems');
fprintf(fileID,'\n\n--------------------\n\n');
r = [m/2; m/4; 3; 2; 1; 0];
maxerr = [maxerr_opt; maxerr_opt_half; maxerr_opt_three; maxerr_opt_two; maxerr_opt_one; maxerr_opt_zero];
matrix = [r,maxerr];
fprintf(fileID,format,transpose(matrix));
fclose(fileID);
end%if
