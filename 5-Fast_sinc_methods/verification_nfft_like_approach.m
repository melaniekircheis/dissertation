% Code file for Figure 5.1

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

addpath('../4-Regularized_Shannon_sampling_formulas'); % cardinal_bspline.m

%% Setup

% Switch flag for saving results to txt-file
save_results = 0; 
if save_results == 1
    fileID = fopen('verification_nfft_like_approach.txt','w');
    format = '%d %1.4e \n';
end%if

% Switch flag for choice of regularized sinc function ('Gauss','bspline','sinh','cKB')
switch_func = 'sinh';

% Set parameters
M = 20; % bandwidth
lambda = 1; % oversampling parameter
L = (1+lambda)*M; % oversampled bandwidth
m = 5; % truncation parameter

%% Init of the problem

% Set non-integer samples in frequency domain
v = -M/2-m:1/32:M/2+m;

% Set samples in spatial domain
l = -L/2:L/2-1;
R = 10^5;
for j = 1:2
switch j
    case 1
        y = (-1/2:1/R:1/2-1/R)';
        text = '$x\in[-\frac 12, \frac 12)$';
    case 2
        y = (-1/2+m/L:(1-2*m/L)/R:1/2-m/L-(1-2*m/L)/R)';
        text = '$x\in[-\frac 12+\frac mL, \frac 12-\frac mL)$';
end%switch

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

%% Comparison of the quality of the two approaches

% Compute approximation using the NFFT approach
approx_nfft = phi(mod((y-l/L)+0.5,1)-0.5)*(exp(2*pi*1i*l'/L*v)./(L*phihat(v)));
err_nfft = max(abs(exp(2*pi*1i*y*v) - approx_nfft));

% Compute approximation using the new NFFT-like approach for bandlimited functions
approx_bandlim = psi(y-l/L)*(exp(2*pi*1i*l'/L*v)./(L*psihat(v).'));
err_bandlim = max(abs(exp(2*pi*1i*y*v) - approx_bandlim));

%% Visualization

% Visualization of the function psihat
figure(1), plot(v,L*psihat(v)); title('$L \hat\psi(v)$'); xlim([-M/2+m/L,M/2-m/L]);
xlabel('$v$'); xticks([-M/2,0,M/2]); xticklabels({'$-\frac M2$','$0$','$\frac M2$'});

% Visualization of the approximation
figure(2); subplot(1,2,j); semilogy(v,err_nfft,'-x',v,err_bandlim,'-o'); xlim([min(v),max(v)]);
xlabel('$v$'); xticks([-M/2,0,M/2]); xticklabels({'$-\frac M2$','$0$','$\frac M2$'});
legend('$\mathbf{BFD}$','$\mathbf{\Psi F D_{\hat\psi}}$','Location','north'); title(text);
sgtitle({'Figure 5.1: Approximation error (5.45) for $P=10^5$ computed',' for (5.46) with $S=32$ using the $\sinh$-type window function (4.62)','as well as the parameters $M=20$, $\lambda=1$, $L=(1+\lambda)M$, and $m=5$',' in the one-dimensional setting $d=1$.'});
colororder(["#5E2316";"#B22D10"]);

%% Generate tables for tikz

if save_results == 1
fprintf(fileID,[text,'\n\n']);
    fprintf(fileID,'NFFT\n\n');
    matrix = [v.',err_nfft.'];
    fprintf(fileID,format,matrix');
    fprintf(fileID,'\nnew approach\n\n');
    matrix = [v.',err_bandlim.'];
    fprintf(fileID,format,matrix');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%if

end%for

if save_results == 1, fclose(fileID); end%if

%% Nested functions

function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end