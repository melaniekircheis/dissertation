% Code file for Figures 4.8 and 4.9

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Set parameters
m = 5; lambda = 0.25; beta = pi*m*lambda/(1+lambda);

% Read known zeros j_{1,n} of the Bessel function J1, see [Watson 80]
j1 = [3.8317,7.0156,10.1735,13.3237,16.4706,19.6159];

% Compute zeros of the integrand
xi = sqrt((j1./beta).^2+1);

%% Initialization of the Bessel function

func = @(w) besselj(1,beta*sqrt(w.^2-1))./sqrt(w.^2-1);
func_subst = @(x) besselj(1,beta*sinh(x));
w1_max = (2+lambda)/lambda; amax = acosh(w1_max);

%% Visualization of the integrands

wj = linspace(-xi(end),xi(end),1000); yj = func(wj);
figure(1); plot(wj,yj); legend('$\frac{J_1(\beta\sqrt{w^2-1})}{\sqrt{w^2-1}}$');
xlabel('$w$'); xlim([-xi(end),xi(end)]); xticks([-flip(xi(1:end-1)),0,xi(1:end-1)]); xticklabels({'$-\xi_5$','','$-\xi_3$','','$-\xi_1$','$0$','$\xi_1$','','$\xi_3$','','$\xi_5$'});
ylim([min(yj)*2,max(yj)*1.1]); yticks([0,besseli(1,beta)]); yticklabels({'$0$','$I_1(\beta)$'})
title('Figure 4.8: The integrand $\frac{J_1(\beta\sqrt{w^2-1})}{\sqrt{w^2-1}}$.');

xj = linspace(0,acosh(w1_max),1000); fj = func_subst(xj);
figure(2); plot(xj,fj); legend('$J_1(\beta \,\sinh x)$'); 
xlabel('$x$'); xlim([0,amax]); xticks([0,amax]); xticklabels({0,'$\mathrm{arcosh}(w_1(\frac{M}{2}))$'});
yticks(0);
title('Figure 4.9: The integrand $J_1(\beta\sinh x)$ on the interval $\big[0,\mathrm{arcosh}\big(w_1\big(\frac M2\big)\big)\big]$.');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('assertion_sinh.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'first integrand\n\n');
matrix = [wj',yj'];
fprintf(fileID,format,matrix.');
fprintf(fileID,'----------------------------------------\n\n');
fprintf(fileID,'second integrand\n\n');
matrix = [xj',fj'];
fprintf(fileID,format,matrix.');
fclose(fileID);
end%if
