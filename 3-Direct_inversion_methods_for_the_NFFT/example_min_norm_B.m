% Code file for Figure 3.9
% Minimization of the Frobenius norm of A'BFD

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Add necessary folders to search path
addpath('..\nfft\matlab\nfft')  % NFFT software -- path needs to be adjusted and next line to be deleted
error('NFFT software needs to be installed on your system first!')
addpath('grids')  % grids
addpath('opt_B')  % optimization methods for sparse matrix B

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;
if (save_results == 1)
filename = 'min_norm.txt';
fileID = fopen(filename,'w');
format = '%g %1.4e \n';
end%if

% Set char vector of possible grids
v = char('Equispaced x jittered','1d-jittered','2d-jittered','Polar','Modified polar','Linogram','Golden angle polar','Golden angle linogram');

% Set initial parameters
z = 2.^(2:7);
M = [12;12];
mode = 'normal_equations';

% Set vector for number of points
NN = zeros(length(z),1);

% Set vectors of norms
norm_fro_B_BSpline = zeros(4,length(z),size(v,1)); norm_2_B_BSpline = zeros(4,length(z),size(v,1));
norm_fro_opt_BSpline = zeros(4,length(z),size(v,1)); norm_2_opt_BSpline = zeros(4,length(z),size(v,1));
norm_fro_opt_Dirichlet = zeros(4,length(z),size(v,1)); norm_2_opt_Dirichlet = zeros(4,length(z),size(v,1));
norm_fro_opt_BSpline_D = zeros(4,length(z),size(v,1)); norm_fro_opt_Dirichlet_D = zeros(4,length(z),size(v,1));
deltaD_BSpline = zeros(4,length(z),size(v,1)); deltaD_Dirichlet = zeros(4,length(z),size(v,1));

% Set vectors of computation time
t_BSpline = zeros(4,length(z),size(v,1));
t_Dirichlet = zeros(4,length(z),size(v,1));

%% Main computations
for s = 4:5 %1:size(v,1) %% test all desired grids
for r = 1:4

% Set parameters
sigma = [ceil(r/2);ceil(r/2)];    % Oversampling factor
m = ceil(sigma.*M);
ind = (mod(m,2)==1);  % if it's odd, make it even
m(ind) = m(ind)+1;
n = (mod(r-1,2)+1)*2;  % Truncation parameter 

for c = 1:length(z)
    %% Set grid
    R = z(c);

    switch s
        case 1
            % Equispaced x jittered grid
            x1 = -0.5:1/R:0.5-1/R;
            x2 = (-0.5+1/(2*R):1/R:0.5-1/(2*R)) + 1/(4*R)*(2*rand(1,R)-1);
            [X1,X2] = ndgrid(x1,x2);
            x = [X1(:) X2(:)];
            clear x1 x2 X1 X2;

        case 2
            % 1d-jittered grid
            x1 = (-0.5+1/(2*R):1/R:0.5-1/(2*R)) + 1/(4*R)*(2*rand(1,R)-1);
            x2 = (-0.5+1/(2*R):1/R:0.5-1/(2*R)) + 1/(4*R)*(2*rand(1,R)-1);
            [X1,X2] = ndgrid(x1,x2);
            x = [X1(:) X2(:)];
            clear x1 x2 X1 X2;

        case 3
            % 2d-jittered grid
            x = (-0.5+1/(2*R):1/R:0.5-1/(2*R));
            [X1,X2] = ndgrid(x,x);
            x = [X1(:) X2(:)];
            x = x + 1/(4*R)*(2*rand(R^2,2)-1);
            clear X1 X2;

        case 4
            % Polar grid
            x = polar(R,2*R);

        case 5
            % Modified polar grid
            x = mpolar(R,2*R);

        case 6
            % Linogram grid
            x = linogram(R,2*R);

        case 7
            % Golden angle polar grid
            x = golden_angle_polar(R,2*R);

        case 8
            % Golden angle linogram grid
            x = golden_angle_linogram(R,2*R);
    end%switch

    x = unique(x,'stable','rows');
    N = size(x,1); % Number of nodes
    NN(c) = N; % Save amount for plotting

    %% Generate original matrices

    % Initialize identity matrix 
    IM = eye(prod(M));

    % Generate matrix A
    k1 = -M(1)/2:M(1)/2-1;
    k2 = -M(2)/2:M(2)/2-1;
    [K1,K2] = ndgrid(k1,k2);
    k = [K1(:) K2(:)];
    clear k1 k2 K1 K2;
    A = exp(2*pi*1i*(x*k.'));

    % Generate matrix F
    l1 = -m(1)/2:m(1)/2-1;
    l2 = -m(2)/2:m(2)/2-1;
    [L1,L2] = ndgrid(l1,l2);
    l = [L1(:) L2(:)];
    clear l1 l2 L1 L2;
    F = exp(2*pi*1i*((l./m')*k.'));
    
    % ------------------------------------
    % Generate matrix D using BSpline
    phihat_1 = phi_hat(k(:,1),M(1),m(1),'BSpline',n);
    phihat_2 = phi_hat(k(:,2),M(2),m(2),'BSpline',n);
    D_BSpline = 1/prod(m) * sparse(1:prod(M),1:prod(M),1./(phihat_1.*phihat_2),prod(M),prod(M));

    % Generate matrix B using BSpline
    B_BSpline = sparse(N,prod(m));
    for j = 1:N
      y = x(j,:)-l./m';
      y = mod(y+0.5,1)-0.5;
      phi_1 = phi(y(:,1),M(1),m(1),'BSpline',n);
      phi_2 = phi(y(:,2),M(2),m(2),'BSpline',n);
      B_BSpline = B_BSpline + sparse(j,1:prod(m),phi_1.*phi_2,N,prod(m));
    end
    
    % ------------------------------------
    % Generate matrix D using the Dirichlet kernel
    D_Dirichlet = 1/prod(m) * sparse(1:prod(M),1:prod(M),1.0,prod(M),prod(M));

    %% Fast computation of optimized matrix B
    [O_BSpline,t_BSpline(r,c,s)] = optB_overdet_2d(N,M,m,x,'BSpline',n,mode);
    [O_Dirichlet,t_Dirichlet(r,c,s)] = optB_overdet_2d(N,M,m,x,'Dirichlet',n,mode);

    %% Computation of Frobenius norms
    norm_fro_B_BSpline(r,c,s) = norm(A'*B_BSpline*F*D_BSpline-IM,'fro');  % using the original matrix with BSpline
    norm_fro_opt_BSpline(r,c,s) = norm(A'*O_BSpline*F*D_BSpline-IM,'fro');  % using the optimized matrix with BSpline
    norm_fro_opt_Dirichlet(r,c,s) = norm(A'*O_Dirichlet*F*D_Dirichlet-IM,'fro');  % using the optimized matrix with Dirichlet

    %% Computation of spectral norms
    norm_2_B_BSpline(r,c,s) = norm(A'*B_BSpline*F*D_BSpline-IM,2);  % using the original matrix with BSpline
    norm_2_opt_BSpline(r,c,s) = norm(A'*O_BSpline*F*D_BSpline-IM,2);  % using the optimized matrix with BSpline
    norm_2_opt_Dirichlet(r,c,s) = norm(A'*O_Dirichlet*F*D_Dirichlet-IM,2);  % using the optimized matrix with Dirichlet
    %% Optimization of the matrix D

    % Set optimized matrix product
    beta_BSpline = A'*O_BSpline*F;
    % Computation of the corresponding optimal diagonal matrix D
    d_BSpline = zeros(prod(M),1);
    for k = 1:prod(M)
        d_BSpline(k) = conj(beta_BSpline(k,k))./(beta_BSpline(:,k)'*beta_BSpline(:,k));
    end%for
    Dtilde_BSpline = sparse(1:prod(M),1:prod(M),d_BSpline,prod(M),prod(M));
    % Computation of Frobenius norm
    norm_fro_opt_BSpline_D(r,c) = norm(A'*O_BSpline*F*Dtilde_BSpline-IM,'fro');
    % Measure the difference between original and optimized D
    deltaD_BSpline(r,c) = max(max(abs(D_BSpline-Dtilde_BSpline)));

    % ------------------------------------
    % Set optimized matrix product
    beta_Dirichlet = A'*O_Dirichlet*F;
    % Computation of the corresponding optimal diagonal matrix D
    d_Dirichlet = zeros(prod(M),1);
    for k = 1:prod(M)
        d_Dirichlet(k) = conj(beta_Dirichlet(k,k))./(beta_Dirichlet(:,k)'*beta_Dirichlet(:,k));
    end%for
    Dtilde_Dirichlet = sparse(1:prod(M),1:prod(M),d_Dirichlet,prod(M),prod(M));
    % Computation of Frobenius norm
    norm_fro_opt_Dirichlet_D(r,c) = norm(A'*O_Dirichlet*F*Dtilde_Dirichlet-IM,'fro');  % using the optimized matrices with BSpline
    % Measure the difference between original and optimized D
    deltaD_Dirichlet(r,c) = max(max(abs(D_Dirichlet-Dtilde_Dirichlet)));

end%for
fprintf(['Done m = ',num2str((mod(r-1,2)+1)*2),', sigma = ',num2str(ceil(r/2)),' %s\n'], datestr(datetime('now'))); 

%% Generate tables for tikz

if (save_results == 1)
fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,[strtrim(v(s,:)),'\n \n']); 
fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,'------------------\n'); 
fprintf(fileID,[' m = ',num2str((mod(r-1,2)+1)*2),' sigma = ',num2str(ceil(r/2)),'\n \n']); 
fprintf(fileID,'------------------\n'); 

fprintf(fileID,'\n Original Frobenius norm with BSpline\n\n'); 
matrix = [NN, norm_fro_B_BSpline(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Optimized Frobenius norm with BSpline\n\n'); 
matrix = [NN, norm_fro_opt_BSpline(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Original spectral norm with BSpline\n\n'); 
matrix = [NN, norm_2_B_BSpline(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Optimized spectral norm with BSpline\n\n'); 
matrix = [NN, norm_2_opt_BSpline(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Computation times for BSpline\n\n'); 
matrix = [NN, t_BSpline(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n-------------------------------\n'); 
fprintf(fileID,'\n Optimized Frobenius norm with Dirichlet kernel\n\n'); 
matrix = [NN, norm_fro_opt_Dirichlet(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Optimized spectral norm with Dirichlet kernel\n\n'); 
matrix = [NN, norm_2_opt_Dirichlet(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Computation times for Dirichlet\n\n'); 
matrix = [NN, t_Dirichlet(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n-------------------------------\n'); 
fprintf(fileID,'\n Frobenius norm using optimized D with BSpline\n\n'); 
matrix = [NN, norm_fro_opt_BSpline_D(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Frobenius norm using optimized D with Dirichlet kernel\n\n'); 
matrix = [NN, norm_fro_opt_Dirichlet_D(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n-------------------------------\n'); 
fprintf(fileID,'\n Difference of norms for BSpline\n\n'); 
matrix = [NN, deltaD_BSpline(r,:,s).'];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'\n Difference of norms for Dirichlet kernel\n\n'); 
matrix = [NN, deltaD_Dirichlet(r,:,s).'];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n');
fprintf(fileID,'--------------------------------------------------------\n'); 
end%if
end%for 

%% Visualization

figure(s)
loglog(NN,norm_fro_B_BSpline(1,:,s),'diamond-',NN,norm_fro_opt_BSpline(1,:,s),'o-',NN,norm_fro_opt_Dirichlet(1,:,s),'square-',NN,norm_fro_B_BSpline(4,:,s),'diamond--',NN,norm_fro_opt_BSpline(4,:,s),'o--',NN,norm_fro_opt_Dirichlet(4,:,s),'square--',NN,norm_2_B_BSpline(1,:,s),'^-',NN,norm_2_opt_BSpline(1,:,s),'-x',NN,norm_2_opt_Dirichlet(1,:,s),'-+', NN,norm_2_B_BSpline(4,:,s),'^--',NN,norm_2_opt_BSpline(4,:,s),'x--',NN,norm_2_opt_Dirichlet(4,:,s),'+--')
xlabel(['$N$']); colororder(["#7E2F8E";"#F98E36";"#4DBEEE";"#7E2F8E";"#F98E36";"#4DBEEE";"#744219";"#D40F0F";"#000000";"#744219";"#D40F0F";"#000000"])
legend('$n_{\mathrm{F}}(\varphi_{\mathrm{B}},2,1)$','$n_{\mathrm{F}}^{\textrm{opt}}(\varphi_{\mathrm{B}},2,1)$','$n_{\mathrm{F}}^{\textrm{opt}}(\varphi_{\mathrm{D}},2,1)$', '$n_{\mathrm{F}}(\varphi_{\mathrm{B}},4,2)$','$n_{\mathrm{F}}^{\textrm{opt}}(\varphi_{\mathrm{B}},4,2)$','$n_{\mathrm{F}}^{\textrm{opt}}(\varphi_{\mathrm{D}},4,2)$','$n_{2}(\varphi_{\mathrm{B}},2,1)$','$n_{2}^{\textrm{opt}}(\varphi_{\mathrm{B}},2,1)$','$n_{2}^{\textrm{opt}}(\varphi_{\mathrm{D}},2,1)$', '$n_{2}(\varphi_{\mathrm{B}},4,2)$','$n_{2}^{\textrm{opt}}(\varphi_{\mathrm{B}},4,2)$','$n_{2}^{\textrm{opt}}(\varphi_{\mathrm{D}},4,2)$','Location','southoutside','NumColumns',2)
title([strtrim(v(s,:)),' grid']);
sgtitle({'Figure 3.9: Frobenius norms (3.96) as well as spectral norms (3.97)',' of the original matrix $B$ (violet) and the optimized matrix $B_{\mathrm{opt}}$ ','generated by Algorithm 3.24 using the $\mathrm B$-spline window function $\varphi_{\mathrm{B}}$ (orange) ','as well as the Dirichlet window function $\varphi_{\mathrm{D}}$ (cyan) from (3.66) ','with $R=2^\mu$, $\mu\in\{2,\dots,7\}$, and $T=2R$ as well as $\b M=(12,12)^T$, ', '$m\in\{2,4\}$, and $\sigma\in\{1,2\}$.'});

end%for
if (save_results == 1), fclose(fileID); end%if