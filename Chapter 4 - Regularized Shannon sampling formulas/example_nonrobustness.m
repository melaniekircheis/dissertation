% Code file for Figure 4.13

%% Setup

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Set bandwidth
a = 7;
N = 2^a;

% Set parameters to study
% n = 1:10; % truncation parameter
n = 2.^(0:15); % truncation parameter
lambda = [0;0.5;1;2]; % oversampling parameter
sigma = 1+lambda; % auxiliary parameter

% Set maximum pertubation parameter
epsilon = 1e-3;

% Initialization of vectors
lower = zeros(length(n),1); 
upper = zeros(length(n),1); 
err = zeros(length(n),length(lambda)); 

%% Error constant

for i1 = 1:length(n)
    lower(i1,:) = 2/pi*log(n(i1)) + 5/4;
    upper(i1,:) = 2/pi*log(n(i1)) + 5/4 + 1/(2*n(i1));
end%if 

%% Reconstruction error

% Initialize fine grid for evaluation of reconstruction error
S = 1e5;
s = (-S:S)';
t = s/S;
err_term = zeros(length(t),1);

% Loop for computation of the error
f = @my_sinc;
for i1 = 1:length(n) 
    % Set truncation parameter
    T = n(i1);
    k = (-T:T)'; % Corresponding index set

    % Set function evaluations
    eps_k = sign(prod(f(pi,1/2-k),2));

    for i2 = 1:length(sigma)
        % Set oversampling
        L = sigma(i2)*N;

        % Setup
        for i3 = 1:length(t)
            x = t(i3,:)-k/L;
            phi = f(L*pi,x); % Sinc 
            phi = prod(phi,2);

            % Evaluation of regularized WKS sums
            err_term(i3) = phi.'*eps_k;
       end%for

        % Computation of reconstruction errors
        err(i1,i2) = norm(err_term,inf);
    end%for
    
    fprintf(['T=',num2str(T),' done %s\n'], datestr(datetime('now')))
end%for

%% Visualization 

figure(1); 
for d=1:3
    subplot(1,3,d), semilogx(n,epsilon.*(lower(:,1)).^d,'--',n,epsilon.*(upper(:,1)).^d,'--',n,epsilon.*(err(:,1)).^d,'-o',n,epsilon.*(err(:,2)).^d,'-square',n,epsilon.*(err(:,3)).^d,'-diamond',n,epsilon.*(err(:,4)).^d,'-*'); 
    xlabel('$T$'); legend('','',['$\lambda=$ ',num2str(lambda(1))],['$\lambda=$ ',num2str(lambda(2))],['$\lambda=$ ',num2str(lambda(3))],['$\lambda=$ ',num2str(lambda(4))],'location','east'); title(['$N=$ ',num2str(N)])
    title(['$d=$ ',num2str(d)]);
    colororder(["#000000";"#000000";"#1f38ab";"#FF007F";"#D95319";"#008080"])
end%for
sgtitle({'Figure 4.13: The maximum approximation error (4.148) as well as its lower and upper bounds (4.29) and (4.27),',' respectively, for $T = 2^c$, $c\in\{0,\dots, 15\}$, and $L= M(1+\lambda)$ with $\lambda \in \{0,0.5,1,2\}$, where $M=128$ and $\varepsilon=10^{-3}$ are chosen.'});

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('non_robustness.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,['\n\n Robustness error for different lambda with N= ',num2str(N),', epsilon= ',num2str(epsilon),'\n']);
fprintf(fileID,'\n---------------------------------------------------------------');
fprintf(fileID,'\n---------------------------------------------------------------\n');
% Loop for several dimensions
for d = 1:3
    fprintf(fileID,['d=',num2str(d)]);
    fprintf(fileID,'\n---------------------------------------------------------------\n\n');
    % Error constants
    fprintf(fileID,'Lower error constant\n');
    matrix = [n.',(lower(:,1)).^d];
    fprintf(fileID,format,matrix.');
    fprintf(fileID,'\n\nUpper error constant\n');
    matrix = [n.',(upper(:,1)).^d];
    fprintf(fileID,format,matrix.');
    % Error 
    fprintf(fileID,'\n-------------------------\n\n');
    fprintf(fileID,'Error');
    for k = 1:length(lambda)
        fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
        matrix = [n.',(err(:,k)).^d];
        fprintf(fileID,format,matrix.');
    end%for
    fprintf(fileID,'\n---------------------------------------------------------------');
    fprintf(fileID,'\n---------------------------------------------------------------\n');
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