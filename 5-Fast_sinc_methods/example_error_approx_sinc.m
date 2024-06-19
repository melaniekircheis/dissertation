% Code file for Figures 5.2 and 5.3

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('..\nfft\matlab\nfft')  % NFFT software
error('NFFT software needs to be installed on your system first!') % path above needs to be adjusted and this line to be deleted

%% Setup

% Set switches
flag_legendre = 1; % Switch for comparison to Gauss-Legendre quadrature
flag_equi = 1; % Switch for comparison to equispaced quadrature
flag_mc = 1; % Switch for comparison to Monte-Carlo quadrature
flag_2d = 1; % Switch for 2-dimensional tests
flag_3d = 1; % Switch for 3-dimensional tests
save_results = 0; % Switch flag for saving results to txt-file

% Set grid points for the comparison
if flag_3d == 1
    R = 10^2;
elseif flag_2d == 1
    R = 10^4; 
else
    R = 10^6;
end%if 
x = (-1:2/R:1-1/R)';

% Set parameters
maxN = 15; % Maximum size of nonharmonic bandwidth
maxnu = 10; % Maximum size for scaling parameter

%% Clenshaw-Curtis quadrature

% Error constant
const_cc = 48/35*(exp(2.^(1:maxN)'*(3/4*pi-(1:maxnu)*log(2)))+exp(-2.^(1:maxN)'*(3/4*pi+(1:maxnu)*log(2))))/2; %% = 48/35*2.^(-n).*cosh(3/4*pi*N);
if flag_2d == 1; const_cc_2d = 2*const_cc; end%if
if flag_3d == 1; const_cc_3d = 3*const_cc; end%if

% Approximation error
err_max_cc = zeros(maxN,maxnu);
if flag_2d == 1; err_max_cc_2d = zeros(maxN,maxnu); end%if
if flag_3d == 1; err_max_cc_3d = zeros(maxN,maxnu); end%if
for j = 1:maxN
    % Set evaluations
    N = 2^j;                 % Degree of sinc kernel
    ev_x = sin(N*pi*x)./(N*pi*x); ev_x(x==0)=1; 
    % Compute approximation
    for nu = 1:maxnu
        n = nu*N;
        k = (0:n)';
        eps = ones(n+1,1); eps([1,end])=1/2;
        alpha = (2./(1-4*(0:n/2).^2)).';
        % Fast compuation of the weights using FFT
        beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
        y = [beta;flip(beta(2:end-1))];
        w_cc = ifft(y);
        w_cc = eps.*real(w_cc(1:n+1));
        % Fast evaluation by means of an adjoint NFFT
        z_cc = 1/2*(cos(k*pi/n));
        plan = nfft(1,R,n+1); 
        plan.x = -z_cc*2*N/R;
        plan.f = w_cc; 
        nfft_adjoint(plan); 
        approx_x_cc = plan.fhat; 
        % Compute approximation error
        err_max_cc(j,nu) = max(abs(ev_x-approx_x_cc));
        % Two-dimensional computations if requested
        if flag_2d == 1
            err_cc_2d = zeros(R,1);
            for r = 1:R
                err_cc_2d(r) = max(abs(ev_x(r)*ev_x-approx_x_cc(r)*approx_x_cc));
            end%for
            err_max_cc_2d(j,nu) = max(err_cc_2d);
        end%if
        % Three-dimensional computations if requested
        if flag_3d == 1
            err_cc_3d = zeros(R,1);
            for r = 1:R
                h = zeros(R,1);
                for p = 1:R
                    h(p) = max(abs(ev_x(p)*ev_x(r)*ev_x-approx_x_cc(p)*approx_x_cc(r)*approx_x_cc));
                end%for
                err_cc_3d(r) = max(h);
            end%for
            err_max_cc_3d(j,nu) = max(err_cc_3d);
        end%if
    end%for
end%for

% Visualization of error constant (with respect to N)
figure(1); subplot(1,2,1); loglog(2.^(1:maxN),const_cc(:,1),'--|',2.^(1:maxN),const_cc(:,2),'--*',2.^(1:maxN),const_cc(:,3),'--^',2.^(1:maxN),const_cc(:,4),'--o',2.^(1:maxN),const_cc(:,5),'--x')
legend('$\nu=1$','$\nu=2$','$\nu=3$','$\nu=4$','$\nu=5$','Location','northwest'); xlabel('$M$');
title('Error constant (5.9)'); xlim([2,32]); xticks([2,4,8,16,32]); ylim([1e-6,1e18]);
set(gca, 'ColorOrder', [14, 65, 61;0, 158, 115;129, 87, 29;188, 234, 28;0, 114, 178]./255, 'NextPlot', 'replacechildren'); % "#0E413D";"#009E73";"#81571D";"#BCEA1C";"#0072B2";

% Visualization of error and constant (with respect to n)
subplot(1,2,2); semilogy(1:maxnu,const_cc(2,:),'--o',1:maxnu,const_cc(4,:),'--o',1:maxnu,const_cc(6,:),'--o',1:maxnu,const_cc(8,:),'--o',1:maxnu,const_cc(10,:),'--o',1:maxnu,err_max_cc(2,:),'-o',1:maxnu,err_max_cc(4,:),'-o',1:maxnu,err_max_cc(6,:),'-o',1:maxnu,err_max_cc(8,:),'-o',1:maxnu,err_max_cc(10,:),'-o')
legend('','','','','','$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Maximum error (5.47)',' and error constant (5.8)'}); xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]);
set(gca, 'ColorOrder', [204, 121, 167;129, 28, 160;213, 94, 0;141, 138, 138;0, 0, 0]./255, 'NextPlot', 'replacechildren'); % "#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"
sgtitle({'Figure 5.2: Error constant (5.9) (dashed) and maximum error (5.47) ','(solid) of the approximation of $\mathrm{sinc}(M\pi x)$, $x\in[-1,1]$, with $R=10^6$ ','for certain bandwidths $M=2^s$, $s\in\{1,\dots,15\}$, using the ','Chebyshev nodes (5.8) and the weigths computed by Algorithm 5.5, ','where $n=cM$ with $c\in\{1,\dots,10\}$.'});

if flag_2d == 1
figure(2); semilogy(1:maxnu,const_cc_2d(2,:),'--o',1:maxnu,const_cc_2d(4,:),'--o',1:maxnu,const_cc_2d(6,:),'--o',1:maxnu,const_cc_2d(8,:),'--o',1:maxnu,const_cc_2d(10,:),'--o',1:maxnu,err_max_cc_2d(2,:),'-o',1:maxnu,err_max_cc_2d(4,:),'-o',1:maxnu,err_max_cc_2d(6,:),'-o',1:maxnu,err_max_cc_2d(8,:),'-o',1:maxnu,err_max_cc_2d(10,:),'-o')
legend('','','','','','$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Error constant (5.9) (dashed) and maximum error (5.47) (solid)','using Clenshaw-Curtis quadrature -- two-dimensional setting'});
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
end%if
if flag_3d == 1
figure(3); semilogy(1:maxnu,const_cc_3d(2,:),'--o',1:maxnu,const_cc_3d(4,:),'--o',1:maxnu,const_cc_3d(6,:),'--o',1:maxnu,const_cc_3d(8,:),'--o',1:maxnu,const_cc_3d(10,:),'--o',1:maxnu,err_max_cc_3d(2,:),'-o',1:maxnu,err_max_cc_3d(4,:),'-o',1:maxnu,err_max_cc_3d(6,:),'-o',1:maxnu,err_max_cc_3d(8,:),'-o',1:maxnu,err_max_cc_3d(10,:),'-o')
legend('','','','','','$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Error constant (5.9) (dashed) and maximum error (5.47) (solid)','using Clenshaw-Curtis quadrature -- three-dimensional setting'});
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
end%if

%% Gauss-Legendre quadrature

% Error constant
if (flag_legendre == 1)
const_leg = (pi./(2*(1:maxnu)*exp(2))).^(2*(1:maxnu).*(2.^(1:maxN)'));

% Approximation error
err_max_leg = zeros(maxN,maxnu);
if flag_2d == 1; err_max_leg_2d = zeros(maxN,maxnu); end%if
if flag_3d == 1; err_max_leg_3d = zeros(maxN,maxnu); end%if
for j = 1:maxN
    % Set evaluations
    N = 2^j;                 % Degree of sinc kernel
    ev_x = sin(N*pi*x)./(N*pi*x); ev_x(x==0)=1; 
    % Compute approximation
    for nu = 1:maxnu
        n = nu*N;
        [z_leg,w_leg] = legpts(n+1,[-1/2,1/2]);
        % Fast evaluation by means of an adjoint NFFT
        plan = nfft(1,R,n+1); 
        plan.x = -z_leg*2*N/R;
        plan.f = w_leg'; 
        nfft_adjoint(plan); 
        approx_x_leg = plan.fhat; 
        % Compute approximation error
        err_max_leg(j,nu) = max(abs(ev_x-approx_x_leg));
        % Two-dimensional computations if requested
        if flag_2d == 1
            err_leg_2d = zeros(R,1);
            for r = 1:R
                err_leg_2d(r) = max(abs(ev_x(r)*ev_x-approx_x_leg(r)*approx_x_leg));
            end%for
            err_max_leg_2d(j,nu) = max(err_leg_2d);
        end%if
        % Three-dimensional computations if requested
        if flag_3d == 1
            err_leg_3d = zeros(R,1);
            for r = 1:R
                h = zeros(R,1);
                for p = 1:R
                    h(p) = max(abs(ev_x(p)*ev_x(r)*ev_x-approx_x_leg(p)*approx_x_leg(r)*approx_x_leg));
                end%for
                err_leg_3d(r) = max(h);
            end%for
            err_max_leg_3d(j,nu) = max(err_leg_3d);
        end%if
    end%for
end%for

% Visualization of error and constant (with respect to n)
if flag_2d == 1
figure(4); semilogy(1:maxnu,err_max_leg_2d(2,:),'-o',1:maxnu,err_max_leg_2d(4,:),'-o',1:maxnu,err_max_leg_2d(6,:),'-o',1:maxnu,err_max_leg_2d(8,:),'-o',1:maxnu,err_max_leg_2d(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Error constant (dashed) and maximum error (5.47) (solid)','using Gauss-Legendre quadrature -- two-dimensional setting'});
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
end%if
if flag_3d == 1
figure(5); semilogy(1:maxnu,err_max_leg_3d(2,:),'-o',1:maxnu,err_max_leg_3d(4,:),'-o',1:maxnu,err_max_leg_3d(6,:),'-o',1:maxnu,err_max_leg_3d(8,:),'-o',1:maxnu,err_max_leg_3d(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Error constant (dashed) and maximum error (5.47) (solid)','using Gauss-Legendre quadrature -- three-dimensional setting'});
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
end%if
end%if

%% Comparison of both approaches in terms of computation times

if (flag_legendre == 1)
t_cc = zeros(maxN,1); t_gl = zeros(maxN,1);
for j = 1:maxN
    % Set parameters
    N = 2^j;                 % Degree of sinc kernel
    
    % Computation of Clenshaw-Curtis points and weights
    n = 4*N;
    tic; k = (0:n)';
    eps = ones(n+1,1); eps([1,end])=1/2;
    alpha = (2./(1-4*(0:n/2).^2)).';
    % Fast compuation of weights using FFT
    beta = zeros(n+1,1); beta(mod(k,2)==0) = alpha;
    y = [beta;flip(beta(2:end-1))];
    w_cc = ifft(y);
    w_cc = eps.*real(w_cc(1:n+1)); t_cc(j) = toc;
    
    % Comutation of Gauss-Legendre points and weights
    n = 2*N;
    tic; [z_leg,w_leg] = legpts(n+1,[-1/2,1/2]); t_gl(j) = toc;
end%for

% Visualization of computation times
figure(6); subplot(1,2,1); semilogy(1:maxnu,const_leg(2,:),'--o',1:maxnu,const_leg(4,:),'--o',1:maxnu,const_leg(6,:),'--o',1:maxnu,const_leg(8,:),'--o',1:maxnu,const_leg(10,:),'--o',1:maxnu,err_max_leg(2,:),'-o',1:maxnu,err_max_leg(4,:),'-o',1:maxnu,err_max_leg(6,:),'-o',1:maxnu,err_max_leg(8,:),'-o',1:maxnu,err_max_leg(10,:),'-o')
legend('','','','','','$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]);
set(gca, 'ColorOrder', [204, 121, 167;129, 28, 160;213, 94, 0;141, 138, 138;0, 0, 0]./255, 'NextPlot', 'replacechildren'); % "#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"
title({'Maximum error (5.47) using the', 'weights mentioned in Remark 5.6'});
subplot(1,2,2); loglog(2.^(1:maxN),t_cc,'-o',2.^(1:maxN),t_gl,'-x')
legend('Clenshaw-Curtis with $n=4M$','Gauss-Legendre with $n=2M$','Location','northwest'); xlabel('$M$');
title({'Time needed to compute the', 'approximation (in seconds)'})
set(gca, 'ColorOrder', [188, 234, 28;0, 158, 115]./255, 'NextPlot', 'replacechildren'); % "#0E413D";"#009E73";"#81571D";"#BCEA1C";"#0072B2";
sgtitle({'Figure 5.3: Comparison of the approximation of $\mathrm{sinc}(M\pi x)$, $x\in[-1,1]$,', 'for certain bandwidths $M=2^s$, $s\in\{1,\dots,15\}$,', 'using the Gauss--Legendre weights mentioned in Remark 5.6', 'and the Clenshaw--Curtis weights computed by Algorithm 5.5.'});
end%if

%% Equispaced quadrature

% Approximation error
if (flag_equi == 1)
err_max_eq = zeros(maxN,maxnu);
if flag_2d == 1; err_max_eq_2d = zeros(maxN,maxnu); end%if
if flag_3d == 1; err_max_eq_3d = zeros(maxN,maxnu); end%if
for j = 1:maxN
    % Set evaluations
    N = 2^j;                 % Degree of sinc kernel
    ev_x = sin(N*pi*x)./(N*pi*x); ev_x(x==0)=1; 
    % Compute approximation
    for nu = 1:maxnu
        n = nu*N;
        z_eq = (-0.5:1/n:0.5)'; w_eq = 1/(n+1)*ones(n+1,1);
        % Fast evaluation by means of an adjoint NFFT
        plan = nfft(1,R,n+1); 
        plan.x = -z_eq*2*N/R;
        plan.f = w_eq; 
        nfft_adjoint(plan); 
        approx_x_eq = plan.fhat; 
        % Compute approximation error
        err_max_eq(j,nu) = max(abs(ev_x-approx_x_eq));
        % Two-dimensional computations if requested
        if flag_2d == 1
            err_eq_2d = zeros(R,1);
            for r = 1:R
                err_eq_2d(r) = max(abs(ev_x(r)*ev_x-approx_x_eq(r)*approx_x_eq));
            end%for
            err_max_eq_2d(j,nu) = max(err_eq_2d);
        end%if
        % Three-dimensional computations if requested
        if flag_3d == 1
            err_eq_3d = zeros(R,1);
            for r = 1:R
                h = zeros(R,1);
                for p = 1:R
                    h(p) = max(abs(ev_x(p)*ev_x(r)*ev_x-approx_x_eq(p)*approx_x_eq(r)*approx_x_eq));
                end%for
                err_eq_3d(r) = max(h);
            end%for
            err_max_eq_3d(j,nu) = max(err_eq_3d);
        end%if
    end%for
end%for

% Visualization of error (with respect to n)
figure(7); semilogy(1:maxnu,err_max_eq(2,:),'-o',1:maxnu,err_max_eq(4,:),'-o',1:maxnu,err_max_eq(6,:),'-o',1:maxnu,err_max_eq(8,:),'-o',1:maxnu,err_max_eq(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
title({'Maximum error (5.47) using equispaced quadrature'});
end%if
if flag_2d == 1
figure(8); semilogy(1:maxnu,err_max_eq_2d(2,:),'-o',1:maxnu,err_max_eq_2d(4,:),'-o',1:maxnu,err_max_eq_2d(6,:),'-o',1:maxnu,err_max_eq_2d(8,:),'-o',1:maxnu,err_max_eq_2d(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Maximum error (5.47) using equispaced quadrature -- two-dimensional setting'});
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
end%if
if flag_3d == 1
figure(9); semilogy(1:maxnu,err_max_eq_3d(2,:),'-o',1:maxnu,err_max_eq_3d(4,:),'-o',1:maxnu,err_max_eq_3d(6,:),'-o',1:maxnu,err_max_eq_3d(8,:),'-o',1:maxnu,err_max_eq_3d(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
title({'Maximum error (5.47) using equispaced quadrature -- three-dimensional setting'});
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
end%if

%% Monte Carlo quadrature

% Approximation error
if (flag_mc == 1)
err_max_mc = zeros(maxN,maxnu);
if flag_2d == 1; err_max_mc_2d = zeros(maxN,maxnu); end%if
if flag_3d == 1; err_max_mc_3d = zeros(maxN,maxnu); end%if
for j = 1:maxN
    % Set evaluations
    N = 2^j;                 % Degree of sinc kernel
    ev_x = sin(N*pi*x)./(N*pi*x); ev_x(x==0)=1; 
    % Compute approximation
    for nu = 1:maxnu
        n = nu*N;
        z_mc = rand(n+1,1)-0.5; w_mc = 1/(n+1)*ones(n+1,1);
        % Fast evaluation by means of an adjoint NFFT
        plan = nfft(1,R,n+1); 
        plan.x = -z_mc*2*N/R;
        plan.f = w_mc; 
        nfft_adjoint(plan); 
        approx_x_mc = real(plan.fhat); 
        % Compute approximation error
        err_max_mc(j,nu) = max(abs(ev_x-approx_x_mc));
        % Two-dimensional computations if requested
        if flag_2d == 1
            err_mc_2d = zeros(R,1);
            for r = 1:R
                err_mc_2d(r) = max(abs(ev_x(r)*ev_x-approx_x_mc(r)*approx_x_mc));
            end%for
            err_max_mc_2d(j,nu) = max(err_mc_2d);
        end%if
        % Three-dimensional computations if requested
        if flag_3d == 1
            err_mc_3d = zeros(R,1);
            for r = 1:R
                h = zeros(R,1);
                for p = 1:R
                    h(p) = max(abs(ev_x(p)*ev_x(r)*ev_x-approx_x_mc(p)*approx_x_mc(r)*approx_x_mc));
                end%for
                err_mc_3d(r) = max(h);
            end%for
            err_max_mc_3d(j,nu) = max(err_mc_3d);
        end%if
    end%for
end%for

% Visualization of error (with respect to n)
figure(10); semilogy(1:maxnu,err_max_mc(2,:),'-o',1:maxnu,err_max_mc(4,:),'-o',1:maxnu,err_max_mc(6,:),'-o',1:maxnu,err_max_mc(8,:),'-o',1:maxnu,err_max_mc(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
title({'Maximum error (5.47) using a Monte Carlo approach'});
end%if
if flag_2d == 1
figure(11); semilogy(1:maxnu,err_max_mc_2d(2,:),'-o',1:maxnu,err_max_mc_2d(4,:),'-o',1:maxnu,err_max_mc_2d(6,:),'-o',1:maxnu,err_max_mc_2d(8,:),'-o',1:maxnu,err_max_mc_2d(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
title({'Maximum error (5.47) using a Monte Carlo approach -- two-dimensional setting'});
end%if
if flag_3d == 1
figure(12); semilogy(1:maxnu,err_max_mc_3d(2,:),'-o',1:maxnu,err_max_mc_3d(4,:),'-o',1:maxnu,err_max_mc_3d(6,:),'-o',1:maxnu,err_max_mc_3d(8,:),'-o',1:maxnu,err_max_mc_3d(10,:),'-o')
legend('$M=4$','$M=16$','$M=64$','$M=256$','$M=1024$','Location','northeast');
xlabel('$c$'); xlim([1,10]); ylim([3e-16,4e0]); colororder(["#CC79A7";"#811CA0";"#D55E00";"#8D8A8A";"#000000"]);
title({'Maximum error (5.47) using a Monte Carlo approach -- three-dimensional setting'});
end%if

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('error_approx_sinc.txt','w');
format = '%d %1.4e \n';
% Clenshaw-Curtis error constant with respect to N
fprintf(fileID,'\n\n Clenshaw-Curtis error constant with respect to N');
for i = 1:maxnu
    fprintf(fileID,['\n\n nu=',num2str(i),'\n\n']);
    matrix = [2.^(1:maxN)',(const_cc(:,i))];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n-------------------------------------------');
% Clenshaw-Curtis error and error constant with respect to nu
fprintf(fileID,'\n\n Clenshaw-Curtis error and error constant with respect to nu');
for i = 1:maxN
    fprintf(fileID,['\n\n N=',num2str(2^i)]);
    % Error constant
    fprintf(fileID,'\n\n Error constant \n');
    matrix = [(1:maxnu)',(const_cc(i,:)).'];
    fprintf(fileID,format,matrix.');
    % Error
    fprintf(fileID,'\n\n Error \n');
    matrix = [(1:maxnu)',(err_max_cc(i,:)).'];
    fprintf(fileID,format,matrix.');
    fprintf(fileID,'\n-------------------------------------------');
end%for
if (flag_legendre == 1)
% Gauss-Legendre error and error constant with respect to nu
fprintf(fileID,'\n-------------------------------------------');
fprintf(fileID,'\n\n Gauss-Legendre error and error constant with respect to nu');
for i = 1:maxN
    fprintf(fileID,['\n\n N=',num2str(2^i)]);
    % Error constant
    fprintf(fileID,'\n\n Error constant \n');
    matrix = [(1:maxnu)',(const_leg(i,:)).'];
    fprintf(fileID,format,matrix.');
    % Error
    fprintf(fileID,'\n\n Error \n');
    matrix = [(1:maxnu)',(err_max_leg(i,:)).'];
    fprintf(fileID,format,matrix.');
    fprintf(fileID,'\n-------------------------------------------');
end%for
% Computation times for quadraure points and weights
fprintf(fileID,'\n-------------------------------------------');
fprintf(fileID,'\n\n Computation times for quadraure points and weights');
fprintf(fileID,'\n\n Clenshaw-Curtis with $n=4N$ \n');
matrix = [2.^(1:maxN)',t_cc];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Gauss-Legendre with $n=2N$ \n');
matrix = [2.^(1:maxN)',t_gl];
fprintf(fileID,format,matrix.');
end%if
fclose(fileID);
end%if