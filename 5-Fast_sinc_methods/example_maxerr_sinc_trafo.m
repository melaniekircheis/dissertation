% Code file for Figure 5.5

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('..\nfft\matlab\nfft')  % NFFT software
addpath('..\nfft\matlab\nnfft')  % NNFFT software
error('NFFT software needs to be installed on your system first!') % paths above need to be adjusted and this line to be deleted

%% Init

% Set parameters
d = 1; % Dimension
maxsize = 13; % Maximum size of nonharmonic bandwidth
rep = 100; % Repetition number for random nodes
eta = 2.5;
Nu = [2;4;6;8];

% Setup test setting
equi_a = false; % Random points - for equispaced points : equi_a = true;
equi_b = false; % Random points - for equispaced points : equi_b = true;

% Set switch flag for saving results to txt-file
save_results = 0;
if (save_results == 1)
fileID = fopen('example_maxerr_sinc_trafo.txt','w');
format = '%d %1.4e \n';
end%if

% Set vectors
err_max_ls = zeros(length(Nu),maxsize-1,rep,length(eta)); err_max_cc = zeros(length(Nu),maxsize-1,rep);
t_naive = zeros(length(Nu),maxsize-1,rep);
t_fast_ls = zeros(length(Nu),maxsize-1,rep); t_fast_cc = zeros(length(Nu),maxsize-1,rep);

%% Fast sinc transform

for o = 1:length(Nu)
    nu = Nu(o);
    for j = 2:maxsize
        N = 2^j; % Nonharmonic bandwidth
        N_star = ceil((1+2*8/N)*N);
        n = nu*N;

        % Precompute least square weights
        plan_ls = fastsinc(d,N,'n',n,'P',eta*n);
        fastsinc_precompute_least_squares(plan_ls);
        
        % Precompute Clenshaw-Curtis weights
        plan_cc = fastsinc(d,N,'n',n);
        fastsinc_precompute_analytic(plan_cc);

        %% Repetitions for random nodes
        for i = 1:rep
            % Set data of the discrete sinc transform
            K = N;
            switch equi_a
                case false % Random points
                    a = rand(K,d)-1/2;
                case true % Equispaced points
                    a = (-K/2:K/2-1)'/K;
            end%switch
            c = rand(K,1); % Random coefficients
            L = N;
            switch equi_b
                case false % Random points
                    b = rand(L,d)-1/2;
                case true % Equispaced points
                    b = (-L/2:L/2-1)'/L;
            end%switch

            % Naive evaluation of the discrete sinc transform
            tic; S = ones(L,K);
            for dim = 1:d % Set up sinc matrix
                C = sin(N*pi*(b(:,dim)-a(:,dim)'))./(N*pi*(b(:,dim)-a(:,dim)'));
                C(isnan(C)) = 1;
                S = S.*C;
            end%for
            prod_naive = S*c; t_naive(o,j-1,i) = toc;

            % Fast approximation using empirical least squares weights
            plan_ls.a = a; plan_ls.b = b; % Set points of the discrete sinc tranform
            plan_ls.c = c; % Set coefficients of the discrete sinc tranform
            prod_fast_ls = fastsinc_trafo(plan_ls); % fast sinc transform
            t_fast_ls(o,j-1,i) = plan_ls.time.t_precompute+plan_ls.time.t_trafo;
            err_max_ls(o,j-1,i) = norm(prod_naive-prod_fast_ls,'Inf');

            % Fast approximation using analytic Clenshaw-Curtis weights
            plan_cc.a = a; plan_cc.b = b; % Set points of the discrete sinc tranform
            plan_cc.c = c; % Set coefficients of the discrete sinc tranform
            prod_fast_cc = fastsinc_trafo(plan_cc); % fast sinc transform
            t_fast_cc(o,j-1,i) = plan_cc.time.t_precompute+plan_cc.time.t_trafo;
            err_max_cc(o,j-1,i) = norm(prod_naive-prod_fast_cc,'Inf');
        end%for

        err_max_ls = max(err_max_ls,[],3); err_max_cc = max(err_max_cc,[],3);
        t_naive = max(t_naive,[],3); t_fast_ls = max(t_fast_ls,[],3); t_fast_cc = max(t_fast_cc,[],3);
    end%for
end%for

%% Visualization

figure(1); subplot(1,2,1); loglog(2.^(2:maxsize),err_max_cc(1,:),'-*',2.^(2:maxsize),err_max_cc(2,:),'-o',2.^(2:maxsize),err_max_cc(3,:),'-^',2.^(2:maxsize),err_max_cc(4,:),'-|');
xlabel('$M$'); xlim([2^5,2^13]); legend('$c=2$','$c=4$','$c=6$','$c=8$')
colororder(["#009E73";"#BCEA1C";"#FF00FF";"#808000"]);
title({'Analytic weights','obtained by Algorithm 5.5'});
subplot(1,2,2); loglog(2.^(2:maxsize),err_max_ls(1,:),'-*',2.^(2:maxsize),err_max_ls(2,:),'-o',2.^(2:maxsize),err_max_ls(3,:),'-^',2.^(2:maxsize),err_max_ls(4,:),'-|');
xlabel('$M$'); xlim([2^5,2^13]); legend('$c=2$','$c=4$','$c=6$','$c=8$')
title({'Least squares weights','obtained by means of (5.7)'});
sgtitle({'Figure 5.5: Maximum approximation error (5.27) of the fast $\mathrm{sinc}$ transform', 'in Algorithm 5.10 for certain bandwidths $M=2^s$, $s=5,\dots,13$,', 'shown for several choices of $n=cM$ with $c\in\{2,4,6,8\}$ when using', 'random nodes $a_k, b_\ell \in \big[-\frac 12,\frac 12\big]$ and random coefficients $c_k\in [0,1]$', 'with $d=1$ and $K=L=M$, compared for different weights $w_j$, $j=0,\dots,n$,', 'obtained by Algorithm 5.5 or the least squares approach (5.7)', 'with $z_j$ in (5.8) each.'});

%% Generate tables for tikz

if (save_results == 1)
% Least squares errors
fprintf(fileID,'Maximum Error LSQR\n\n');
for k = 1:length(Nu)
    fprintf(fileID,['n=',num2str(Nu(k)),'M\n\n']);
    matrix = [2.^(2:maxsize)',err_max_ls(k,:).'];
    fprintf(fileID,format,transpose(matrix));
    fprintf(fileID,'\n\n');
end%for
fprintf(fileID,'-------------------------------------------------\n\n');
% Clenshaw-Curtis errors
fprintf(fileID,'Maximum Error Clenshaw-Curtis\n\n');
for k = 1:length(Nu)
    fprintf(fileID,['n=',num2str(Nu(k)),'M\n\n']);
    matrix = [2.^(2:maxsize)',err_max_cc(k,:).'];
    fprintf(fileID,format,transpose(matrix));
    fprintf(fileID,'\n\n');
end%for
fclose(fileID);
end%if