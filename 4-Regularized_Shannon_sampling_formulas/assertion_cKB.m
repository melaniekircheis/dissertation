% Code file for Figures 4.11 and 4.12

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Set parameters
m = 2; lambda = 1; beta = pi*m*lambda/(1+lambda);

%% Initialization of the Bessel function

func = @(w) sinh(beta*sqrt(1-w.^2))./(beta*sqrt(1-w.^2))-my_sinc(beta,w);

% Approximate the zeros 
tmax = sqrt((20./beta).^2+1);
tj = linspace(-tmax,tmax,1e5); zj = func(tj);
ind = find(abs(zj)<3e-5);
xi = unique(round(tj(ind)*100)/100);

%% Visualization 

wj = linspace(-xi(end),xi(end),1000); yj = func(wj);
figure(1); plot(wj,yj); legend('$\frac{\sinh\big(\beta \sqrt{1 - w^2}\,\big)}{\beta \sqrt{1 - w^2}} - \mathrm{sinc}(\beta w)$');
xlabel('$w$'); xlim([-xi(end),xi(end)]); xticks(sort([-xi,0])); xticklabels({'$-\xi_5$','','$-\xi_3$','','$-\xi_1$','$0$','$\xi_1$','','$\xi_3$','','$\xi_5$'});
yticks(0); yticklabels({'$0$'})
title('Figure 4.11: The integrand $\left(\frac{\sinh\big(\beta \sqrt{1 - w^2}\,\big)}{\beta \sqrt{1 - w^2}} - \mathrm{sinc}(\beta w)\right)$.');

x = linspace(0,30,200);
h = besseli(0,x)-struvem(0,x)+2/pi*sinint(x)-1;
figure(2); plot(x,h,x,zeros(200,1),'k-.',x,9/20*ones(200,1),'k-.')
legend('$I_0(\beta) - {\textbf L}_0(\beta) - 1 + \frac{2}{\pi}\,\mathrm{Si}(\beta)$');
xlabel('$\beta$'); xlim([0,30]); xticks([0,10,20,30]);
ylim([-0.05,0.5]); yticks([0,0.45]); yticklabels({'$0$','$\frac{9}{20}$'});
title('Figure 4.12: Visualization of $\big[ I_0(\beta) - {\textbf L}_0(\beta) - 1 + \frac{2}{\pi}\,\mathrm{Si}(\beta) \big] \in \left[0,\,\frac{9}{20}\right)$ for suitable $\beta>0$.');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('assertion_cKB.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'first assertion\n\n');
matrix = [wj',yj'];
fprintf(fileID,format,matrix.');
fprintf(fileID,'----------------------------------------\n\n');
fprintf(fileID,'second assertion\n\n');
matrix = [x',h'];
fprintf(fileID,format,matrix.');
fprintf(fileID,'-------------------------\n');
matrix = [[0;30],zeros(2,1)];
fprintf(fileID,format,matrix.');
fprintf(fileID,'-------------------------\n');
matrix = [[0;30],9/20*ones(2,1)];
fprintf(fileID,format,matrix.');
fclose(fileID);
end%if

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end

function f=struvem(v,x,n)
% Drew (2024). struvem.m (https://www.mathworks.com/matlabcentral/fileexchange/27218-struvem-m), MATLAB Central File Exchange. Retrieved February 19, 2024.
% 
% Calculates the Modified Struve Function L_v(x) and n is the length of
% the series calculation (n=100 if unspecified)
%
% from: Abramowitz and Stegun: Handbook of Mathematical Functions
% 		http://www.math.sfu.ca/~cbm/aands/page_498.htm
% 
if nargin<3
n=100;
end
k=0:n;
x=x(:)';
k=k(:);
xx=repmat(x,length(k),1);
kk=repmat(k,1,length(x));
TOP=1;
BOT=gamma(kk+1.5).*gamma(kk+v+1.5);
RIGHT=(xx./2).^(2.*kk+v+1);
FULL=TOP./BOT.*RIGHT;
f=sum(FULL);
end%function