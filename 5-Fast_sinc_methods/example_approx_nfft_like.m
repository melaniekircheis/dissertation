% Code file for Figure 5.8

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('../4-Regularized_Shannon_sampling_formulas'); % cardinal_bspline.m

%% Setup

% Switch flag for saving results to txt-file
save_results = 0; 

% Switch flag for choice of regularized sinc function ('Gauss','bspline','sinh','cKB')
switch_func = 'sinh';

% Set parameters
lambda = 1; % oversampling parameter
m = 5; % truncation parameter
Mmax = 50; % maximum bandwidth
err_nonequi = zeros(Mmax,2); % Init error vector

%% Main loop
for j = 1:Mmax
M = 20*j;
k = -M/2:M/2-1;
L = (1+lambda)*M;

%% Initialization of the bandlimited function and its Fourier transform

% f = @(x) sqrt(M)*my_sinc(M*pi,x);  % sinc
% fhat = @(v) 1/sqrt(M)*rectangularPulse(-M/2,M/2,v); % characteristic function
f = @(x) (my_sinc(M*pi/2,x)).^2;  % sinc^2
fhat = @(v) (abs(v)<=M/2).*2/M.*(1-abs(2*v/M)); % triangular pulse
% f = @(x) sqrt(4*M/5)*(my_sinc(M*pi,x)+my_sinc(M*pi,(x-1))/2);  % usual test function
% fhat = @(v) (abs(v)<=M/2).*(sqrt(4/(5*M)))+(abs(v)<=M/2).*(sqrt(4/(5*M)).*1/2*exp(-2*pi*1i*v)); % corresponding Fourier transform
% n = M/(1+lambda); f = @(x) (n+M)/2.*my_sinc((n+M)/2*pi,x).*my_sinc((M-n)/2*pi,x);  % psi_lin
% fhat = @(v) (abs(v)<=n/2) + (n/2<abs(v) & abs(v)<=M/2).*(1-(2*abs(v)-n)/(M-n));  % psihat_lin
% n = M/(1+lambda); f = @(x) (abs(x)<=eps).*((M+n)/2) + (abs(x)>eps).*fillmissing((12./((M-n)^3*pi^4*x.^4).*(cos(n*pi*x)-cos(M*pi*x))-6./((M-n)^2*pi^3*x.^3).*(sin(n*pi*x)+sin(M*pi*x))),'constant',(M+n)/2);  % psi_cub
% fhat = @(v) (abs(v)<=n/2) + (n/2<abs(v) & abs(v)<=M/2).*(48/(3*(M-n)^3)*(abs(v)-M/2).^2.*(abs(v)-(3*n-M)/4));  % psihat_cub
% n = M/(1+lambda); f = @(x) (abs(x)<=eps).*((M-n)/pi*sin((n*pi)/(M-n))+((M-n)*(2*pi*cos((n*pi)/(M-n))+4*sin((M*pi)/(M-n))-5*sin((n*pi)/(M-n))-sin(((-2*M+n)*pi)/(M-n))))/(8*pi)) + (abs(x)>eps).*fillmissing(-1./(x.^2*(M-n)^2-1).*(n/2*my_sinc(n*pi,x)+M/2*my_sinc(M*pi,x)),'constant',(M-n)/pi*sin((n*pi)/(M-n))+((M-n)*(2*pi*cos((n*pi)/(M-n))+4*sin((M*pi)/(M-n))-5*sin((n*pi)/(M-n))-sin(((-2*M+n)*pi)/(M-n))))/(8*pi));  % psi_cos
% fhat = @(v) (abs(v)<=n/2) + (n/2<abs(v) & abs(v)<=M/2).*(1/2.*(1+cos(pi*(2*abs(v)-n)/(M-n))));  % psihat_cos

%% Precomputation of regularized sinc function

switch switch_func
    case 'Gauss' % Gaussian regularized sinc function
        sigma = sqrt(m./(pi*L*(L-M))); % Set variance of Gaussian
        phi = @(x) (abs(x)<=m/L).*exp(-x.^2/(2*sigma^2));
        psi = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*exp(-x.^2/(2*sigma^2));
        phihat = @(v) sqrt(2*pi)*sigma*exp(-2*pi^2*sigma^2*v.^2);
        psihat = @(v) (erf(sqrt(2)*pi*sigma*(v(:)+L/2))-erf(sqrt(2)*pi*sigma*(v(:)-L/2)))/(2*L);
    case 'bspline' % B-spline regularized sinc function
        s = ceil((m+1)/2);
        phi = @(x) (abs(x)<=m/L).*(cardinal_bspline(s*L*x/m,2*s)./cardinal_bspline(0,2*s));
        psi = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*(cardinal_bspline(s*L*x/m,2*s)./cardinal_bspline(0,2*s));
        phihat = @(u) m./(s*L*cardinal_bspline(0,2*s))*my_sinc(pi*m/(s*L),u).^(2*s);
        psihat = @(v) arrayfun(@(j) 1/L*integral(phihat,v(j)-L/2,v(j)+L/2),(1:length(v)).');
    case 'sinh' % sinh-type regularized sinc function
        beta = m*pi*lambda./(1+lambda); % Set parameter of sinh
        phi = @(x) (abs(x)<=m/L).*(sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta));
        psi = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*(sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta));
        rad = @(u) sqrt(u.^2/L^2-(lambda)^2/((2+2*lambda)^2));
        phihat = @(u) beta./(2*L*sinh(beta)).*fillmissing(besselj(1,2*pi*m*rad(u))./rad(u),'constant',2*pi*m*1/2);
        psihat = @(v) arrayfun(@(j) 1/L*integral(phihat,v(j)-L/2,v(j)+L/2),(1:length(v)).');
    case 'cKB' % cKB regularized sinc function
        beta = m*pi*lambda./(1+lambda); % Set parameter of cKB
        phi = @(x) (abs(x)<=m/L).*((besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1));
        psi = @(x) (abs(x)<=m/L).*my_sinc(L*pi,x).*((besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1));
        subst = @(u) 2*pi*m/(beta*L)*u;
        phihat = @(u) 2*m./((besseli(0,beta)-1)*L).*(sinh(beta*sqrt(1-subst(u).^2))./(beta*sqrt(1-subst(u).^2))-my_sinc(beta,subst(u)));
        psihat = @(v) arrayfun(@(j) 1/L.*integral(phihat,v(j)-L/2,v(j)+L/2),(1:length(v)).');
end%switch

%% Reconstruction of equispaced samples

% Compute original samples
Omega = L;
l = -Omega/2:Omega/2-1;
nu = f(l'/L);

% Reconstruction
theta_hat = fhat(k')./psihat(k);
theta = 1/Omega*exp(2*pi*1i*l'/L*k)*theta_hat*M/length(k); % FFT

% Error computation
err_equi = norm(theta-nu,'Inf');

%% Reconstruction of nonequispaced samples

% Init of nonequispaced points
N = M/2;
% x = (-1/2:1/N:1/2-1/N)'*(1-2*m/L);
% x = sort((rand(N,1)-1/2));
x = sort(cos((0:N-1)'*pi/N))*(1/2-m/L);
% x = (-1/2:2/N:1/2-2/N)'; x = sort([x;x+1/(2*N)])*(1/2-m/L);
fx = f(x); %exp(2*pi*1i*x*k)*fhat(k');

% Multiplication with sparse matrix
Psi = psi(x-l/L);
f_tilde = Psi*theta;
f_tilde_nu = Psi*nu; err_nu = norm(fx-f_tilde_nu,'Inf'); % Comparison to original samples

% Comparison to NFFT
f_tilde_ndft = exp(2*pi*1i*x*k)*fhat(k'); % NDFT
h = (fhat(k)./phihat(k))'; % Multiplication with diagonal matrix D
h = 1/Omega*exp(2*pi*1i*l'/L*k)*h*M/length(k); % FFT
B = phi(mod((x-l/L)+0.5,1)-0.5);
f_tilde_nfft = B*h; % Multiplication with sparse matrix B

% Error computation
err_nonequi(j,:) = [norm(fx-f_tilde,'Inf'), norm(fx-f_tilde_nfft,'Inf')];
end%for

%% Visualization

% Visualization of the samples of the Fourier transform
figure(1); plot(k,fhat(k),'.'); xlabel('$k$'); title('$\hat f(k)$');

% Visualiazation of the pointwise reconstruction of the function f
figure(2); plot(x,fx,'bsquare',x,real(f_tilde_nfft),'x',x,real(f_tilde),'o',x,f_tilde_nu,'g.');
legend('$\mathbf{f}$','$\mathbf{BFD \hat f}$','$\mathbf{\Psi F D_{\hat\psi} \hat f}$','$\mathbf{\Psi \nu}$');
xlabel('$x_j$'); xticks([-1/2+m/L,0,1/2-m/L]); xticklabels({'$-\frac  12+\frac mL$','$0$','$\frac 12-\frac mL$'});
colororder(["#5E2316";"#B22D10"]); title('Reconstruction of $f(x_j)$');
figure(3); semilogy(x,abs(fx-real(f_tilde_nfft)),'x',x,abs(fx-real(f_tilde)),'o',x,abs(fx-f_tilde_nu),'g.');
legend('$\mathbf{BFD \hat f}$','$\mathbf{\Psi F D_{\hat\psi} \hat f}$','$\mathbf{\Psi \nu}$');
xlabel('$x_j$'); xticks([-1/2+m/L,0,1/2-m/L]); xticklabels({'$-\frac 12+\frac mL$','$0$','$\frac 12-\frac mL$'});
colororder(["#5E2316";"#B22D10"]); title('Error of reconstruction of $f(x_j)$');

% Visualization of maximum error
figure(4); semilogy(20:20:20*Mmax,err_nonequi(:,2),'-x',20:20:20*Mmax,err_nonequi(:,1),'-o'); 
legend('$\mathbf{BFD}$','$\mathbf{\Psi FD_{\hat\psi}}$');
xlabel('$M$');
title({'Figure 5.8: Maximum approximation error (5.51) of the NFFT-like procedure for','bandlimited functions in Algorithm 5.14 and the classical NFFT in Algorithm 2.2','using the $\sinh$-type window function (4.62) computed for the','function $f(x) = \mathrm{sinc}^2 \big(\frac{M}{2}\pi x\big)$ using several bandwidth',' parameters $M \in \{20,40,\dots,1000\}$ and the scaled Chebyshev nodes (5.52) ','with $N = \frac{M}{2}$, $m=5$, $M_\sigma=L=(1+\lambda)M$, as well as $\lambda= 1$ and $d=1$.'});
colororder(["#5E2316";"#B22D10"])

%% Generate tables for tikz

if save_results==1
fileID = fopen('approx_nfft_like_approach.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'NFFT\n\n');
matrix = [(20:20:1000)',err_nonequi(:,2)];
fprintf(fileID,format,matrix');
fprintf(fileID,'\n-------------------\n\n');
fprintf(fileID,'\nnew approach\n\n');
matrix = [(20:20:1000)',err_nonequi(:,1)];
fprintf(fileID,format,matrix');
fclose(fileID);
end%if

%% Nested functions

function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end
