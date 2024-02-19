% Code file for Figure 4.14

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Switch flag for plotting part (a) or (b) of Figure 4.14
flag_experiment = 'a';

% Set bandwidth
a = 7;
N = 2^a;

% Set test function
f = @(x) sqrt(N)*my_sinc(N*pi,x);
% f = @(x) sqrt(3*N/4)*(my_sinc(N*pi/2,x)).^2;
% f = @(x) sqrt(4*N/5)*(my_sinc(N*pi,x)+my_sinc(N*pi,(x-1))/2);

% Set parameters to study
switch flag_experiment
    case 'a'
        n = 10:10:500; % truncation parameter
    case 'b'
        n = 2.^(0:15); % truncation parameter
end%switch
lambda = [0.5;1;2]; % oversampling parameter
sigma = 1+lambda; % auxiliary parameter

% Initialization of error vectors
const = zeros(length(n),length(lambda)); 
err = zeros(length(n),length(lambda)); 

%% Error constant

% Computation of error constant for the linear frequency window function
for k = 1:length(lambda) 
    const(:,k) = real(sqrt(2*sigma(k)*N/3).*2*(1+lambda(k))./(pi^2*lambda(k)).*(n-sigma(k)*N).^(-3/2));
end%for

%% Reconstruction error

% Initialization of a fine grid for evaluation of reconstruction error
S = 1e5;
s = (-S:S)';
t = s/S;

% Initialization of vectors
Rm = zeros(length(t),1);
err_tilde = zeros(length(t),1);

% Loop for computation of the error
% Set function evaluations
ft = f(t);

for i2 = 1:length(n) 
    % Set truncation parameter
    T = n(i2);
    j = (-T:T)'; % Corresponding index set
        
    for k = 1:length(sigma)
        % Set oversampling
        L = sigma(k)*N;

        % Set function evaluations
        fj = f(j/L);

        % Setup
        for i3 = 1:length(t)
            x = t(i3,:)-j/L;
            psi = (N+L)/2.*my_sinc((N+L)/2*pi,x).*my_sinc((L-N)/2*pi,x);

            % Evaluation of regularized WKS sums
            Rm(i3) = psi.'/L*fj;
        end%for

        % Computation of reconstruction errors
        err(i2,k) = norm(Rm-ft,inf);
    end%for
    
    fprintf(['T=',num2str(n(i2)),' done %s\n'], datestr(datetime('now')))
end%for

%% Visualization 

switch flag_experiment
    case 'a'
        figure(1); semilogy(n,const(:,1),'--',n,err(:,1),'-s',[192,192],[1e+03,1e-14],'-.',n,const(:,2),'--',n,err(:,2),'-diamond',[256,256],[1e+03,1e-14],'-.',n,const(:,3),'--',n,err(:,3),'-*',[384,384],[1e+03,1e-14],'-.'); 
        xlabel('$T$'); xticks([0 100 200 300 400 500]); legend('',['$\lambda=$ ',num2str(lambda(1))],'','',['$\lambda=$ ',num2str(lambda(2))],'','',['$\lambda=$ ',num2str(lambda(3))],''); title(['$N=$ ',num2str(N)])
        title({'Figure 4.14 (a): Maximum approximation error (4.149) (solid) and error constant (4.40) (dashed)', 'using the linear frequency window $\psi_{\mathrm{lin}}$ from (4.45) in (4.38) for the function', ' $f(x) = \sqrt{M} \,\mathrm{sinc}(M \pi x)$ with $M=128$, $\lambda\in\{0.5,1,2\}$ and $T\in\{10,20,\dots,500\}$.'})
        colororder(["#FF007F";"#FF007F";"#FF007F";"#D95319";"#D95319";"#D95319";"#008080";"#008080";"#008080"])
        ylim([5e-13,2e+2])
    case 'b'
        figure(1); loglog(n,const(:,1),'--',n,err(:,1),'-s',[192,192],[1e+03,1e-14],'-.',n,const(:,2),'--',n,err(:,2),'-diamond',[256,256],[1e+03,1e-14],'-.',n,const(:,3),'--',n,err(:,3),'-*',[384,384],[1e+03,1e-14],'-.'); 
        xlabel('$T$'); legend('',['$\lambda=$ ',num2str(lambda(1))],'','',['$\lambda=$ ',num2str(lambda(2))],'','',['$\lambda=$ ',num2str(lambda(3))],''); title(['$N=$ ',num2str(N)])
        title({'Figure 4.14 (b): Maximum approximation error (4.149) (solid) and error constant (4.40) (dashed)', 'using the linear frequency window $\psi_{\mathrm{lin}}$ from (4.45) in (4.38) for the function', ' $f(x) = \sqrt{M} \,\mathrm{sinc}(M \pi x)$ with $M=128$, $\lambda\in\{0.5,1,2\}$ and $T = 2^c$, $c \in \{0,\dots,15\}$.'})
        colororder(["#FF007F";"#FF007F";"#FF007F";"#D95319";"#D95319";"#D95319";"#008080";"#008080";"#008080"])
        ylim([5e-13,2e+2])
end%switch

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('error_psihat_lin.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,['\n\n Reconstruction error for different lambda with N= ',num2str(N),'\n']);
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Error constant
fprintf(fileID,'Error constant');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',const(:,k)];
    fprintf(fileID,format,matrix.');
end%for
% Reconstruction error 
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Error');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',err(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fclose(fileID);
end%if
fprintf('\n Finished %s\n', datestr(datetime('now')))

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end