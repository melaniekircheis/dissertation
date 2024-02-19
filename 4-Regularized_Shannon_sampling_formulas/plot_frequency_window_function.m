% Code file for Figures 4.1, 4.2 and 4.3

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

%% Plot frequency window functions

% Evaluation of the characteristic function
k = (-3*N/2:3*N/S:3*N/2)'; % Set evaluation points
chi = zeros(length(k),1); chi(abs(k)<=N/2) = 1;

% Evaluation of linear frequency window function
ind = N/2<=abs(k) & abs(k)<=L/2;
psihat_lin = chi; psihat_lin(ind) = 1-(2*abs(k(ind))-N)/(L-N);

% Evaluation of cubic frequency window function
psihat_cub = chi; psihat_cub(ind) = 48/(3*(L-N)^3)*(abs(k(ind))-L/2).^2.*(abs(k(ind))-(3*N-L)/4);

% Evaluation of raised cosine frequency window function
psihat_cos = chi; psihat_cos(ind) = 1/2.*(1+cos(pi*(2*abs(k(ind))-N)/(L-N)));

% Evaluation of convolutional frequency window function
psihat_conv = zeros(length(k),10);
k_ind2 = k(ind)-sign(k(ind))*N/2;
for n = 1:10
psihat_conv(:,n) = chi; psihat_conv(ind,n) = 1./cardinal_bspline(0,n+1).*cardinal_bspline((2*n)/(L-N)*k_ind2,n+1);
end%for

% Comparison of the frequency window function
figure(1); plot(k,chi,'k--',k,psihat_lin,k,psihat_cub,k,psihat_cos)
legend('$\mathbf 1_{[-\frac N2,\frac N2]}(v)$','$\hat\psi_{\mathrm{lin}}(v)$','$\hat\psi_{\mathrm{cub}}(v)$','$\hat\psi_{\cos}(v)$','Location','south');  
title('Comparison of the frequency window functions $\hat\psi$'); set(gca,'fontsize',14);
xlabel('$v$'); xlim([-5*N/4,5*N/4]); xticks([-L/2,-N/2,0,N/2,L/2]); xticklabels({'$-\frac{L}{2}$','$-\frac{N}{2}$','$0$','$\frac{N}{2}$','$\frac{L}{2}$'});
yticks([0,1/(2),1]); 

%% Plot the corresponding functions in time domain

% Evaluation of sinc function
R = 2*m/L; % range of points
x = (-R:2*R/S:R-2*R/S)'; % set evaluation points
f_sinc = my_sinc(L*pi,x); % sinc

% Evaluation of linear window function
psi_lin = (N+L)/2.*my_sinc((N+L)/2*pi,x).*my_sinc((L-N)/2*pi,x);

% Evaluation of cubic window function
psi_cub = 12./((L-N)^3*pi^4*x.^4).*(cos(N*pi*x)-cos(L*pi*x))-6./((L-N)^2*pi^3*x.^3).*(sin(N*pi*x)+sin(L*pi*x));
psi_cub(abs(x)<3e-6) = (L+N)/2;

% Evaluation of raised cosine window function
ind = abs(abs(x)-1/(L-N))<eps;
psi_cos = -1./(x.^2*(L-N)^2-1).*(N/2*my_sinc(N*pi,x)+L/2*my_sinc(L*pi,x));
psi_cos(ind) = (L-N)/pi*sin((N*pi)/(L-N))+((L-N)*(2*pi*cos((N*pi)/(L-N))+4*sin((L*pi)/(L-N))-5*sin((N*pi)/(L-N))-sin(((-2*L+N)*pi)/(L-N))))/(8*pi);

% Evaluation of convolutional window function
psi_conv = zeros(length(x),10);
for n = 1:10
psi_conv(:,n) = (N+L)/2.*my_sinc((N+L)/2*pi,x).*(my_sinc((L-N)/(2*n)*pi,x)).^n;
end%for
psi_conv_inf = (N+L)/2.*my_sinc((N+L)/2*pi,x);

% Comparison of the functions in time domain
figure(2); p = plot(x,f_sinc,'-k',x,psi_lin/L,x,psi_cub/L,x,psi_cos/L); p(2).Color = [0.8500 0.3250 0.0980];
legend('$\mathrm{sinc}(L \pi x)$','$\psi_{\mathrm{lin}}(x)$','$\psi_{\mathrm{cub}}(x)$','$\psi_{\cos}(x)$');
title('Comparison of the functions $\psi$ in time domain'); set(gca,'fontsize',14);
xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{10}{L}$','$-\frac{5}{L}$','0','$\frac{5}{L}$','$\frac{10}{L}$'})
yticks([0 0.5 1])

%% Generate the paper plots

% Visualization for the linear frequency window function
figure(3); subplot(1,2,1); plot(k,chi,'k--',k,psihat_lin)
legend('$\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{lin}}(v)$','Location','south');  
title('$\hat\psi_{\mathrm{lin}}$ in (4.44)'); set(gca,'fontsize',11);
xlabel('$v$'); xlim([-90,90]); xticks([-L/2,-N/2,0,N/2,L/2]); xticklabels({'$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M{2}$','$\frac{L}{2}$'});
yticks([0,1/(2),1]); 
subplot(1,2,2), plot(x,f_sinc,'-k',x,psi_lin/L); 
legend('$\mathrm{sinc}(L \pi x)$','$\frac 1L\psi_{\mathrm{lin}}(x)$');
title('$\frac 1L\psi_{\mathrm{lin}}$ in (4.45)'); set(gca,'fontsize',11);
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{10}{L}$','$-\frac{5}{L}$','0','$\frac{5}{L}$','$\frac{10}{L}$'})
yticks([0 0.5 1])
sgtitle('Figure 4.1: The frequency window function (4.44) and its scaled inverse Fourier transform (4.45).');

% Visualization for the cubic frequency window function
figure(4); subplot(1,2,1); plot(k,chi,'k--',k,psihat_cub,k,psihat_cos,'b:')
legend('$\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{cub}}(v)$','$\hat\psi_{\cos}(v)$','Location','south');  
title('$\hat\psi_{\mathrm{cub}}$ and $\hat\psi_{\cos}$'); set(gca,'fontsize',11);
xlabel('$v$'); xlim([-90,90]); xticks([-L/2,-N/2,0,N/2,L/2]); xticklabels({'$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M}{2}$','$\frac{L}{2}$'});
yticks([0,1/(2),1]); 
subplot(1,2,2), plot(x,f_sinc,'-k',x,psi_cub/L,x,psi_cos/L,'b:'); 
legend('$\mathrm{sinc}(L \pi x)$','$\frac 1L\psi_{\mathrm{cub}}(x)$','$\frac 1L\psi_{\cos}(x)$');
title('$\frac 1L\psi_{\mathrm{cub}}$ and $\frac 1L\psi_{\cos}$'); set(gca,'fontsize',11);
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{10}{L}$','$-\frac{5}{L}$','0','$\frac{5}{L}$','$\frac{10}{L}$'})
yticks([0 0.5 1])
sgtitle('Figure 4.2: The frequency window functions (4.47) and (4.50), and their scaled inverse Fourier transforms.');

% Visualization for the convolutional frequency window functions
figure(5); subplot(1,2,1); plot(k,chi,'k--',k,psihat_conv)
legend('$\chi_{[-\frac M2,\frac M2]}(v)$','$\hat\psi_{\mathrm{conv,1}}(v)$','$\hat\psi_{\mathrm{conv,2}}(v)$','$\hat\psi_{\mathrm{conv,3}}(v)$','$\hat\psi_{\mathrm{conv,4}}(v)$','$\hat\psi_{\mathrm{conv,5}}(v)$','$\hat\psi_{\mathrm{conv,6}}(v)$','$\hat\psi_{\mathrm{conv,7}}(v)$','$\hat\psi_{\mathrm{conv,8}}(v)$','$\hat\psi_{\mathrm{conv,9}}(v)$','$\hat\psi_{\mathrm{conv,10}}(v)$','Location','south');  
title('$\hat\psi_{\mathrm{conv,n}}$ in (4.53)'); set(gca,'fontsize',11);
xlabel('$v$'); xlim([-90,90]); xticks([-L/2,-N/2,0,N/2,L/2]); xticklabels({'$-\frac{L}{2}$','$-\frac{M}{2}$','$0$','$\frac{M}{2}$','$\frac{L}{2}$'});
yticks([0,1/(2),1]); 
subplot(1,2,2), plot(x,f_sinc,'-k',x,psi_conv/L,x,psi_conv_inf/L); 
legend('$\mathrm{sinc}(L \pi x)$','$\frac 1L\psi_{\mathrm{conv,1}}(x)$','$\frac 1L\psi_{\mathrm{conv,2}}(x)$','$\frac 1L\psi_{\mathrm{conv,3}}(x)$','$\frac 1L\psi_{\mathrm{conv,4}}(x)$','$\frac 1L\psi_{\mathrm{conv,5}}(x)$','$\frac 1L\psi_{\mathrm{conv,6}}(x)$','$\frac 1L\psi_{\mathrm{conv,7}}(x)$','$\frac 1L\psi_{\mathrm{conv,8}}(x)$','$\frac 1L\psi_{\mathrm{conv,9}}(x)$','$\frac 1L\psi_{\mathrm{conv,10}}(x)$','$\frac{L+M}{2L} \mathrm{sinc}\big(\frac{L+M}{2} \pi x\big)$');
title('$\frac 1L\psi_{\mathrm{conv,n}}$ in (4.55)'); set(gca,'fontsize',11);
xlabel('$x$'); xlim([-R,R]); xticks([-R,-R/2,0,R/2,R]); xticklabels({'-$\frac{10}{L}$','$-\frac{5}{L}$','0','$\frac{5}{L}$','$\frac{10}{L}$'})
yticks([0 0.5 1])
sgtitle('Figure 4.3: The convolutional frequency window functions~\mbox{$\hat\psi_{\mathrm{conv},n}$} in (4.53) with \mbox{$\rho_n(v) =\frac{2 n}{L - M}\,B_n\big(\frac{2n}{L - M} v \big)$} for~\mbox{$n\in\{1,\dots,6\}$}, and their scaled inverse Fourier transforms (4.55).');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('frequency_windows.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,['xmin=',num2str(-3*N/2),', xmax=',num2str(3*N/2),', xtick = {',num2str(-L/2),',',num2str(-N/2),',0,',num2str(N/2),',',num2str(L/2),'}, \n\n']);
fprintf(fileID,'chi\n\n');
matrix = [k,chi];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Linear frequency window function\n\n');
matrix = [k,psihat_lin];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Cubic frequency window function\n\n');
matrix = [k,psihat_cub];
fprintf(fileID,format,matrix.');
fprintf(fileID,'\n\n Raised cosine frequency window function\n\n');
matrix = [k,psihat_cos];
fprintf(fileID,format,matrix.');
for n=1:10
fprintf(fileID,['\n\n Convolutional frequency window function with n=',num2str(n),'\n\n']);
matrix = [k,psihat_conv(:,n)];
fprintf(fileID,format,matrix.');
end%for
% 
fprintf(fileID,'\n\n -----------------------------------------------\n\n');
fprintf(fileID,['xmin=',num2str(-10*R),', xmax=',num2str(10*R),', xtick = {',num2str(-10*R),',',num2str(-10*R/2),',0,',num2str(10*R/2),',',num2str(10*R),'}, \n\n']);
fprintf(fileID,'\n\n sinc\n\n');
matrix = [10*x,f_sinc];
fprintf(fileID,format,matrix.');
fprintf(fileID,'Linear frequency window function\n\n');
matrix = [10*x,psi_lin/L];
fprintf(fileID,format,matrix.');
fprintf(fileID,'Cubic frequency window function\n\n');
matrix = [10*x,psi_cub/L];
fprintf(fileID,format,matrix.');
fprintf(fileID,'Raised cosine frequency window function\n\n');
matrix = [10*x,psi_cos/L];
fprintf(fileID,format,matrix.');
for n=1:10
fprintf(fileID,['\n\n Convolutional frequency window function with n=',num2str(n),'\n\n']);
matrix = [10*x,psi_conv(:,n)/L];
fprintf(fileID,format,matrix.');
end%for
fclose(fileID);
end%if

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end%function