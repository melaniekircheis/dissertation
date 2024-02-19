% Code file for Figures 4.4, 4.5, 4.7 and 4.10

%% Setup

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Switch flag for saving results to txt-file
save_results = 0;

%% Initialization of parameters

N = 80; % bandwidth
m = 5; % truncation parameter
lambda = 1; % oversampling parameter
L = (lambda+1)*N; % oversampling
S = 1e3; % number of points

%% Plot window functions in spatial domain

% Evaluation of sinc function
R = 2*m/L; % range of points
x = (-R:2*R/S:R-2*R/S)'; % set evaluation points
f_sinc = my_sinc(L*pi,x); % sinc

% Evaluation of rectangular regularized sinc function
phi_rect = abs(x)<=m/L; % rectangular window
psi_rect = f_sinc.*phi_rect; % rectangular regularized sinc function

% Evaluation of Gaussian regularized sinc function
sigma = 1/N*sqrt(m./(pi*(1+lambda)*lambda)); % set variance of Gaussian
phi_Gauss = exp(-x.^2/(2*sigma^2)); % Gaussian window
psi_Gauss = f_sinc.*phi_Gauss; % Gaussian regularized sinc function

% Evaluation of B-spline regularized sinc function
s = ceil((m+1)/2); % set shape parameter
phi_B = cardinal_bspline(s*L*x/m,2*s)./cardinal_bspline(0,2*s); % modified B-Spline window
psi_B = f_sinc.*phi_B; % B-spline regularized sinc function

% Evaluation of sinh-type regularized sinc function
beta = m*pi*lambda./(1+lambda); % set shape parameter
phi_sinh = sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta); phi_sinh = real(phi_sinh); % sinh-type window
psi_sinh = f_sinc.*phi_sinh; % sinh-type regularized sinc function

% Evaluation of the continuous Kaiser-Bessel regularized sinc function
beta = m*pi*lambda./(1+lambda); % set shape parameter
phi_cKB = (besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1); phi_cKB = real(phi_cKB); % cKB window
psi_cKB = f_sinc.*phi_cKB; % cKB regularized sinc function

% Comparison of the regularized sinc functions
figure(1); plot(x,f_sinc,'-k',x,psi_rect,x,psi_Gauss,x,psi_B,x,psi_sinh,x,psi_cKB); 
legend('$\mathrm{sinc}(L \pi x)$','$\psi_\mathrm{rect}(x)$','$\psi_\mathrm{Gauss}(x)$','$\psi_\mathrm{B}(x)$','$\psi_\mathrm{sinh}(x)$','$\psi_\mathrm{cKB}(x)$');
title('Comparison of the regularized $\mathrm{sinc}$ functions $\psi$');
xlabel('$t$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])

%% Plot the corresponding functions in frequency domain

% Evaluation of the characteristic function
k = (-L:2*L/S:L)'; % Set evaluation points
chi = zeros(length(k),1); chi(abs(k)<N/2) = 1;

% Evaluation of Gaussian window function
sigma = 1/N*sqrt(m./(pi*(1+lambda)*lambda)); % set variance of Gaussian
psihat_Gauss = (erf(sqrt(2)*pi*sigma*(k+L/2))-erf(sqrt(2)*pi*sigma*(k-L/2)))/(2);

% Evaluation of modified B-spline window function
s = ceil((m+1)/2); % set shape parameter
psihat_B = zeros(length(k),1);
func = @(u) my_sinc(pi*m/(s*L),u).^(2*s);
for i = 1:length(k)
    psihat_B(i) = integral(func,k(i)-L/2,k(i)+L/2);
end%for
psihat_B = m*psihat_B./(s*L*cardinal_bspline(0,2*s)); 

% Evaluation of sinh-type window function
beta = m*pi*lambda./(1+lambda); % set shape parameter
psihat_sinh = zeros(length(k),1);
rad = @(u) sqrt(u.^2/L^2-(lambda)^2/((2+2*lambda)^2));
func = @(u) besselj(1,2*pi*m*rad(u))./rad(u);
for i = 1:length(k)
    psihat_sinh(i) = integral(func,k(i)-L/2,k(i)+L/2);
end%for
psihat_sinh = beta.*psihat_sinh./(2*L*sinh(beta));

% Evaluation of the continuous Kaiser-Bessel window function
beta = m*pi*lambda./(1+lambda); % set shape parameter
psihat_cKB = zeros(length(k),1);
subst = @(u) 2*pi*m/(beta*L)*u;
func = @(u) sinh(beta*sqrt(1-subst(u).^2))./(beta*sqrt(1-subst(u).^2))-my_sinc(beta,subst(u));
for i = 1:length(k)
    psihat_cKB(i) = integral(func,k(i)-L/2,k(i)+L/2);
end%for
psihat_cKB = 2*m./((besseli(0,beta)-1)*L).*psihat_cKB;

% Comparison of the functions in frequency domain
figure(2); plot(k,chi,'k--',k,psihat_Gauss,k,psihat_B,k,psihat_sinh,k,psihat_cKB); 
legend('$\frac 1L\chi_{[-\frac L2,\frac L2]}(v)$','$\hat\psi_{\mathrm{Gauss}}(v)$','$\hat\psi_{\mathrm{B}}(v)$','$\hat\psi_{\mathrm{sinh}}(v)$','$\hat\psi_{\mathrm{cKB}}(v)$'); 
title('Comparison of the function $\hat\psi$ in frequency domain');
xlim([-L,L]); xticks([-L,-L/2,0,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$0$','$\frac{L}{2}$','$L$'}); 
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'});

%% Generate the paper plots

% Visualization for the Gaussian window function
figure(3); subplot(1,2,1); p = plot(x,f_sinc,'-k',x,phi_Gauss,'--',x,psi_Gauss); p(2).Color = [0.75 0.75 0.75]; p(3).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\varphi_{\mathrm{Gauss}}(x)$','$\psi_\mathrm{Gauss}(x)$');
title('$\psi_{\mathrm{Gauss}}$ in (4.89)');
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])
subplot(1,2,2); plot(k,chi,'k--',k,psihat_Gauss)
legend('$\frac 1L\,\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{Gauss}}(v)$','Location','south'); 
title('$\hat\psi_{\mathrm{Gauss}}$ in (4.90)');
xlabel('$v$'); xlim([-L,L]); xticks([-L,-L/2,-N/2,0,N/2,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M}{2}$','$\frac{L}{2}$','$L$'});
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'});
sgtitle('Figure 4.4: The Gaussian regularized $\mathrm{sinc}$ function $\psi_{\mathrm{Gauss}}$ as well as its Fourier transform $\hat\psi_{\mathrm{Gauss}}$ with $m=5$ and $\alpha=\frac{1}{L}\sqrt{\frac{10}{\pi}}$.');

% Visualization for the modified B-spline window function
figure(4); subplot(1,2,1); p = plot(x,f_sinc,'-k',x,phi_B,'--',x,psi_B); p(2).Color = [0.75 0.75 0.75]; p(3).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\varphi_{\mathrm{B}}(x)$','$\psi_{\mathrm{B}}(x)$');
title('$\psi_{\mathrm{B}}$ in (4.96)');
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])
subplot(1,2,2); plot(k,chi,'k--',k,psihat_B)
legend('$\frac 1L\,\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{B}}(v)$','Location','south'); 
title('$\hat\psi_{\mathrm{B}}$ in (4.97)');
xlabel('$v$'); xlim([-L,L]); xticks([-L,-L/2,-N/2,0,N/2,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M}{2}$','$\frac{L}{2}$','$L$'});
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'}); 
sgtitle('Figure 4.5: The B--spline regularized $\mathrm{sinc}$ function $\psi_{\mathrm{B}}$ as well as its Fourier transform $\hat\psi_{\mathrm{B}}$ with $m=5$ and $s=3$.');

% Visualization for the sinh-type window function
figure(5); subplot(1,2,1); p = plot(x,f_sinc,'-k',x,phi_sinh,'--',x,psi_sinh); p(2).Color = [0.75 0.75 0.75]; p(3).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\varphi_{\sinh}(x)$','$\psi_\mathrm{sinh}(x)$');
title('$\psi_{\mathrm{sinh}}$ in (4.105)');
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])
subplot(1,2,2); plot(k,chi,'k--',k,psihat_sinh)
legend('$\frac 1L\,\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{sinh}}(v)$','Location','south'); 
title('$\hat\psi_{\mathrm{sinh}}$ in (4.106)');
xlabel('$v$'); xlim([-L,L]); xticks([-L,-L/2,-N/2,0,N/2,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M}{2}$','$\frac{L}{2}$','$L$'});
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'});
sgtitle('Figure 4.7: The sinh-type regularized $\mathrm{sinc}$ function $\psi_{\sinh}$ as well as its Fourier transform $\hat\psi_{\sinh}$ with $m=5$ and $\beta=\frac{5}{2}\pi$.');

% Visualization for the cKB window function
figure(6); subplot(1,2,1); p = plot(x,f_sinc,'-k',x,phi_cKB,'--',x,psi_cKB); p(2).Color = [0.75 0.75 0.75]; p(3).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\varphi_{\mathrm{cKB}}(x)$','$\psi_{\mathrm{cKB}}(x)$');
title('$\psi_{\mathrm{cKB}}$ in (4.113)');
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{2m}{L}$','$-\frac{m}{L}$','0','$\frac{m}{L}$','$\frac{2m}{L}$'})
yticks([0 0.5 1])
subplot(1,2,2); plot(k,chi,'k--',k,psihat_cKB)
legend('$\frac 1L\,\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{cKB}}(v)$','Location','south'); 
title('$\hat\psi_{\mathrm{cKB}}$ in (4.114)');
xlabel('$v$'); xlim([-L,L]); xticks([-L,-L/2,-N/2,0,N/2,L/2,L]); xticklabels({'$-L$','$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M}{2}$','$\frac{L}{2}$','$L$'});
ylim([0,1]); yticks([0,1/(2),1]); yticklabels({'$0$','$\frac{1}{2L}$','$\frac 1L$'}); 
sgtitle('Figure 4.10: The continuous Kaiser--Bessel regularized $\mathrm{sinc}$ function $\psi_{\mathrm{cKB}}$ as well as its Fourier transform $\hat\psi_{\mathrm{cKB}}$ with $m=5$ and $\beta=\frac{5}{2}\pi$.');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('time_windows.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,['xmin=',num2str(-L),', xmax=',num2str(L),', xtick = {',num2str(-L),',',num2str(-L/2),',',num2str(-N/2),',0,',num2str(N/2),',',num2str(L/2),',',num2str(L),'}, \n\n']);
fprintf(fileID,'chi\n\n');
matrix = [k,chi];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Gaussian window function\n\n');
matrix = [k,psihat_Gauss];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n modified B-spline window function\n\n');
matrix = [k,psihat_B];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n sinh-type window function\n\n');
matrix = [k,psihat_sinh];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Continuous Kaiser-Bessel window function\n\n');
matrix = [k,psihat_cKB];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n -----------------------------------------------\n\n');
fprintf(fileID,['xmin=',num2str(-10*R),', xmax=',num2str(10*R),', xtick = {',num2str(-10*R),',',num2str(-10*R/2),',0,',num2str(10*R/2),',',num2str(10*R),'}, \n\n']);
fprintf(fileID,'sinc\n\n');
matrix = [10*x,f_sinc];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Gaussian window function\n\n');
matrix = [10*x,phi_Gauss];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Gaussian regularized sinc function\n\n');
matrix = [10*x,psi_Gauss];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Modified B-Spline window function\n\n');
matrix = [10*x,phi_B];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n B-spline regularized sinc function\n\n');
matrix = [10*x,psi_B];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n sinh-type window function\n\n');
matrix = [10*x,phi_sinh];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n sinh-type regularized sinc function\n\n');
matrix = [10*x,psi_sinh];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Continuous Kaiser-Bessel window function\n\n');
matrix = [10*x,phi_cKB];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n cKB regularized sinc function\n\n');
matrix = [10*x,psi_cKB];
fprintf(fileID,format,matrix.');
fclose(fileID);
end%if

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end%function