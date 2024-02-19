% Code file for Figures 4.15, 4.24 and 4.25

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Switch flag for the rational approximation of Natterer (very time consuming !!)
rat_flag = 0;

% Set bandwidth
a = 8;
N = 2^a;

% Set test function
% f = @(x) sqrt(N)*my_sinc(N*pi,x);
% f = @(x) sqrt(3*N/4)*(my_sinc(N*pi/2,x)).^2;
f = @(x) sqrt(4*N/5)*(my_sinc(N*pi,x)+my_sinc(N*pi,(x-1))/2);

% Set parameters to study
n = 2:10; % truncation parameter
lambda = [1/2;1;2]; % oversampling parameter
sigma = 1+lambda; % auxiliary parameter

% Initialization of error vectors
err_shannon = zeros(length(n),length(lambda)); 
err_lin = zeros(length(n),length(lambda)); 
err_cub = zeros(length(n),length(lambda)); 
err_cos = zeros(length(n),length(lambda)); 
err_conv2 = zeros(length(n),length(lambda)); 
err_conv3 = zeros(length(n),length(lambda)); 
err_conv4 = zeros(length(n),length(lambda)); 
err_conv_inf = zeros(length(n),length(lambda)); 
err_rat = zeros(length(n),length(lambda));
err_rect = zeros(length(n),length(lambda)); 
err_Gauss = zeros(length(n),length(lambda)); 
err_B = zeros(length(n),length(lambda)); 
err_sinh = zeros(length(n),length(lambda)); 
err_cKB = zeros(length(n),length(lambda)); 

% Initialization of vectors for error constants
const_shannon = zeros(length(n),length(lambda)); 
const_lin = zeros(length(n),length(lambda)); 
const_cub = zeros(length(n),length(lambda)); 
const_rect = zeros(length(n),length(lambda)); 
const_Gauss = zeros(length(n),length(lambda)); 
const_B = zeros(length(n),length(lambda)); 
const_sinh = zeros(length(n),length(lambda)); 
const_cKB = zeros(length(n),length(lambda)); 

% Initialization of vectors for computation time
time_shannon = zeros(length(n),length(lambda)); 
time_lin = zeros(length(n),length(lambda)); 
time_cub = zeros(length(n),length(lambda)); 
time_cos = zeros(length(n),length(lambda)); 
time_conv2 = zeros(length(n),length(lambda)); 
time_conv3 = zeros(length(n),length(lambda)); 
time_conv4 = zeros(length(n),length(lambda)); 
time_conv_inf = zeros(length(n),length(lambda)); 
time_rect = zeros(length(n),length(lambda)); 
time_Gauss = zeros(length(n),length(lambda)); 
time_B = zeros(length(n),length(lambda)); 
time_sinh = zeros(length(n),length(lambda)); 
time_cKB = zeros(length(n),length(lambda)); 
if (rat_flag == 1), time_rat = zeros(length(n),length(lambda)); end%if

%% Error constants

for k = 1:length(lambda)     
    % Computation of error constants for frequency regularizations
    const_shannon(:,k) = sqrt(2*(1+lambda(k))*N/pi^2).*(n).^(-1/2); % classical Shannon sampling sums
    const_lin(:,k) = sqrt(2*(1+lambda(k))*N/3).*2*(1+lambda(k))./(pi^2*lambda(k)).*(n).^(-3/2); % linear frequency window function
    const_cub(:,k) = sqrt(2*(1+lambda(k))*N/5).*14*(1+lambda(k))^2./(pi^3*lambda(k)^2).*(n).^(-5/2); % cubic frequency window function

    % Computation of error constants for time regularizations
    const_rect(:,k) = sqrt(2*sigma(k)*N)./pi.*sqrt(1./n+1./n.^2);
    % Gaussian window function
    const_Gauss(:,k) = (sqrt(2*n*sigma(k)*N)+sqrt(2*(sigma(k)*N)*lambda(k)*(1+n)))./(pi*n*sqrt(lambda(k))).*exp(-n*pi*lambda(k)/(2*sigma(k)));
    % B-Spline window function
    s = ceil((n+1)/2); % set shape parameter
    const_B(:,k) = (3*sqrt(s.*N)./(sqrt(2)*pi*(2*s-1))).*exp(-n.*(log(pi*n.*(lambda(k)))-log(2*s*sigma(k))));
    % sinh-type window function
    beta = pi*n.*(lambda(k))./(sigma(k)); % set shape parameter
    const_sinh(:,k) = sqrt(N).*exp(-beta);
    % Continuous Kaiser-Bessel window function
    const_cKB(:,k) = sqrt(N).*(7*pi*n*lambda(k).*(1+lambda(k)+4*n*lambda(k)))./(4*(sigma(k))^2).*exp(-beta);
end%for

%% Reconstruction error

% Initialization of a fine grid for evaluation of reconstruction error
S = 1e5;
s = (-S:S)';
t = s/S;

% Initialization of vectors for reconstructions
Rm_shannon = zeros(length(t),1);
Rm_lin = zeros(length(t),1);
Rm_cub = zeros(length(t),1);
Rm_cos = zeros(length(t),1);
Rm_conv2 = zeros(length(t),1);
Rm_conv3 = zeros(length(t),1);
Rm_conv4 = zeros(length(t),1);
Rm_conv_inf = zeros(length(t),1);
Rm_rect = zeros(length(t),1);
Rm_Gauss = zeros(length(t),1);
Rm_B = zeros(length(t),1);
Rm_sinh = zeros(length(t),1);
Rm_cKB = zeros(length(t),1);
if (rat_flag == 1), Rm_rat = zeros(length(t),1); end%if

% Loop for computation of the error
% Set function evaluations for comparison
ft = f(t);

for i2 = 1:length(n) 
    % Set truncation parameter
    m = n(i2);

    for k = 1:length(sigma)
        % Set oversampling
        L = sigma(k)*N;
        T = m+L;
        j = (-T:T)'; % Corresponding index set

        % Set equispaced function evaluations
        fj = f(j/L);

        % Setup
        for i3 = 1:length(t)
            x = t(i3,:)-j.'/L; % evaluation points
            phi = my_sinc(L*pi,x); % sinc function
            ind_delta = (abs(x)-m/L<=eps); % characteristic function
        
            % Evaluation of classical Shannon sampling sums
            tic; 
            Rm_shannon(i3) = phi*fj; 
            time_shannon(i2,k) = time_shannon(i2,k) + toc;

            % Evaluation of linear frequency regularization
            tic; 
            psi_lin = (N+L)/2.*my_sinc((N+L)/2*pi,x).*my_sinc((L-N)/2*pi,x);
            Rm_lin(i3) = psi_lin/L*fj; 
            time_lin(i2,k) = time_lin(i2,k) + toc;

            % Evaluation of cubic frequency regularization
            tic; 
            psi_cub = 12./((L-N)^3*pi^4*x.^4).*(cos(N*pi*x)-cos(L*pi*x))-6./((L-N)^2*pi^3*x.^3).*(sin(N*pi*x)+sin(L*pi*x));
            psi_cub(abs(x)<3e-6) = (L+N)/2;
            Rm_cub(i3) = psi_cub/L*fj; 
            time_cub(i2,k) = time_cub(i2,k) + toc;

            % Evaluation of raised cosine regularization
            ind = abs(abs(x)-1/(L-N))<eps;
            psi_cos = -1./(x.^2*(L-N)^2-1).*(N/2*my_sinc(N*pi,x)+L/2*my_sinc(L*pi,x));
            psi_cos(ind) = (L-N)/pi*sin((N*pi)/(L-N))+((L-N)*(2*pi*cos((N*pi)/(L-N))+4*sin((L*pi)/(L-N))-5*sin((N*pi)/(L-N))-sin(((-2*L+N)*pi)/(L-N))))/(8*pi);
            Rm_cos(i3) = psi_cos/L*fj; 
            time_cos(i2,k) = time_cos(i2,k) + toc;

            % Evaluation of convolution based regularizations
            tic; psi_conv2 = (N+L)/2.*my_sinc((N+L)/2*pi,x).*(my_sinc((L-N)/(2*2)*pi,x).^2);
            Rm_conv2(i3) = psi_conv2/L*fj; 
            time_conv2(i2,k) = time_conv2(i2,k) + toc;
            tic; psi_conv3 = (N+L)/2.*my_sinc((N+L)/2*pi,x).*(my_sinc((L-N)/(2*3)*pi,x).^3);
            Rm_conv3(i3) = psi_conv3/L*fj; 
            time_conv3(i2,k) = time_conv3(i2,k) + toc;
            tic; psi_conv4 = (N+L)/2.*my_sinc((N+L)/2*pi,x).*(my_sinc((L-N)/(2*4)*pi,x).^4);
            Rm_conv4(i3) = psi_conv4/L*fj; 
            time_conv4(i2,k) = time_conv4(i2,k) + toc;
            % Evaluation of convolutional window function
            tic; psi_conv_inf = (N+L)/2.*my_sinc((N+L)/2*pi,x);
            Rm_conv_inf(i3) = psi_conv_inf/L*fj; 
            time_conv_inf(i2,k) = time_conv_inf(i2,k) + toc;

            % Evaluation of rectangular regularization
            tic;
            psi_rect = phi; psi_rect(~ind_delta) = 0;
            Rm_rect(i3) = psi_rect*fj; 
            time_rect(i2,k) = time_rect(i2,k) + toc;

            % Evaluation of Gaussian regularization
            tic; mu = 1/N*sqrt(m./(pi*lambda(k)*(sigma(k)))); % Set variance of Gaussian
            psi_Gauss = phi.*exp(-x.^2/(2*mu.^2)); psi_Gauss(~ind_delta) = 0;
            Rm_Gauss(i3) = psi_Gauss*fj; 
            time_Gauss(i2,k) = time_Gauss(i2,k) + toc;

            % Evaluation of B-spline regularization
            tic; s = ceil((m+1)/2); 
            psi_B = phi./cardinal_bspline(0,2*s).*cardinal_bspline(s*L*x/m,2*s); psi_B(~ind_delta) = 0;
            Rm_B(i3) = psi_B*fj; 
            time_B(i2,k) = time_B(i2,k) + toc;

            % Evaluation of sinh-type regularization
            tic; beta = m*pi*lambda(k)./(1+lambda(k));
            if beta==0
                psi_sinh = phi.*sqrt(1-(L*x/m).^2);
            else
                psi_sinh = phi.*sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta); 
            end
            psi_sinh(~ind_delta) = 0;       
            Rm_sinh(i3) = psi_sinh*fj; 
            time_sinh(i2,k) = time_sinh(i2,k) + toc;

            % Evaluation of continuous Kaiser-Bessel regularization
            tic; beta = m*pi*lambda(k)./(1+lambda(k));
            if beta==0
                psi_cKB = phi.*(1-(L*x/m).^2);
            else
                psi_cKB = phi.*(besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1);
            end
            psi_cKB(~ind_delta) = 0;
            Rm_cKB(i3) = psi_cKB*fj; 
            time_cKB(i2,k) = time_cKB(i2,k) + toc;

            if (rat_flag == 1)
            % Evaluation of the rational regularization by Natterer
            tic; 
%             psi_rat = phi.*natterer((L-N)/2*x); % different test function
            psi_rat = (N+L)/(2*L)*my_sinc((N+L)*pi/2,x).*natterer((L-N)/4*x);
            Rm_rat(i3) = psi_rat*fj; time_rat(i2,k) = time_rat(i2,k) + toc;
            end%if

        end%for

        % Computation of reconstruction errors
        err_shannon(i2,k) = norm(Rm_shannon-ft,inf); 
        err_lin(i2,k) = norm(Rm_lin-ft,inf); 
        err_cub(i2,k) = norm(Rm_cub-ft,inf); 
        err_cos(i2,k) = norm(Rm_cos-ft,inf); 
        err_conv2(i2,k) = norm(Rm_conv2-ft,inf); 
        err_conv3(i2,k) = norm(Rm_conv3-ft,inf); 
        err_conv4(i2,k) = norm(Rm_conv4-ft,inf); 
        err_conv_inf(i2,k) = norm(Rm_conv_inf-ft,inf); 
        if (rat_flag == 1), err_rat(i2,k) = norm(Rm_rat-ft,inf); end%if
        err_rect(i2,k) = norm(Rm_rect-ft,inf);
        err_Gauss(i2,k) = norm(Rm_Gauss-ft,inf);
        err_B(i2,k) = norm(Rm_B-ft,inf);
        err_sinh(i2,k) = norm(Rm_sinh-ft,inf);
        err_cKB(i2,k) = norm(Rm_cKB-ft,inf);
    end%for

    fprintf(['m=',num2str(n(i2)),' done %s\n'], datestr(datetime('now')))
end%for

%% Visualization 

T = sigma*N+(2:10);
for ind_lambda = 1:length(lambda)
figure(1); subplot(1,3,ind_lambda); semilogy(n,err_shannon(:,ind_lambda),'-|',n,const_shannon(:,ind_lambda),'--',n,err_lin(:,ind_lambda),'-o',n,const_lin(:,ind_lambda),'--',n,err_rat(:,ind_lambda),'-*',n,err_cub(:,ind_lambda),'-square',n,const_cub(:,ind_lambda),'--',n,err_cos(:,ind_lambda),'-diamond',n,err_conv2(:,ind_lambda),'-hexagram',n,err_conv3(:,ind_lambda),'-x',n,err_conv4(:,ind_lambda),'-^',n,err_conv_inf(:,ind_lambda),':.'); 
legend('$\mathrm{sinc}(L\pi\cdot)$','','$\psi_{\mathrm{lin}}$','','$\psi_{\mathrm{rat}}$','$\psi_{\mathrm{cub}}$','','$\psi_{\cos}$','$\psi_{\mathrm{conv,2}}$','$\psi_{\mathrm{conv,3}}$','$\psi_{\mathrm{conv,4}}$','$\psi_{\mathrm{conv,\infty}}$'); 
xlabel('$m$'); title(['$\lambda=$ ',num2str(lambda(ind_lambda))])
colororder(["#808080";"#808080";"#EDB120";"#EDB120";"#f7f54d";"#4DBEEE";"#4DBEEE";"#115f7d";"#D95319";"#FE5BAC";"#7E2F8E";"#000000"])
end%for
sgtitle({'Figure 4.15: Maximum approximation error (solid) and error constants (dashed) using ','classical Shannon sampling sums compared to regularizations (4.36) with several frequency window functions for the ','function (4.151), where $M=256$, $T=L+m$ with $m\in\{2, 3, \ldots, 10\}$, and $\lambda\in\{0.5,1,2\}$.'});

for ind_lambda = 1:length(lambda)
figure(2); subplot(1,3,ind_lambda); semilogy(n,err_rect(:,ind_lambda),'-|',n,const_rect(:,ind_lambda),'--',n,err_Gauss(:,ind_lambda),'-o',n,const_Gauss(:,ind_lambda),'--',n,err_B(:,ind_lambda),'-*',n,const_B(:,ind_lambda),'--',n,err_sinh(:,ind_lambda),'-diamond',n,const_sinh(:,ind_lambda),'--',n,err_cKB(:,ind_lambda),'-square',n,const_cKB(:,ind_lambda),'--'); 
legend('$\varphi_{\mathrm{const}}$','','$\varphi_{\mathrm{Gauss}}$','','$\varphi_{\mathrm{B}}$','','$\varphi_{\sinh}$','','$\varphi_{\mathrm{cKB}}$',''); 
xlabel('$m$'); title(['$\lambda=$ ',num2str(lambda(ind_lambda))])
colororder(["#808080";"#808080";"#801919";"#801919";"#10520e";"#10520e";"#77AC30";"#77AC30";"#00008B";"#00008B"])
end%for
sgtitle({'Figure 4.24: Maximum approximation error (4.152) and error constant (4.79) using $\varphi\in\{\varphi_{\mathrm{const}},\varphi_{\mathrm{Gauss}},\varphi_{\mathrm{B}},\varphi_{\sinh},\varphi_{\mathrm{cKB}}\}$','for the function (4.151) with $M=256$, $d=1$, as well as $m\in\{2, 3, \ldots, 10\}$, and $\lambda\in\{0.5,1,2\}$.'});

for ind_lambda = 1:length(lambda)
figure(3); subplot(1,3,ind_lambda); semilogy(n,err_lin(:,ind_lambda),'-o',n,const_lin(:,ind_lambda),'--',n,err_conv2(:,ind_lambda),'-hexagram',n,const_cub(:,ind_lambda),'--',n,err_conv3(:,ind_lambda),'-x',n,err_conv4(:,ind_lambda),'-^',n,err_sinh(:,ind_lambda),'-diamond',n,const_sinh(:,ind_lambda),'--',n,err_cKB(:,ind_lambda),'-square',n,const_cKB(:,ind_lambda),'--'); 
legend('$\psi_{\mathrm{lin}}$','','$\psi_{\mathrm{conv,2}}$','','$\psi_{\mathrm{conv,3}}$','$\psi_{\mathrm{conv,4}}$','$\psi_{\sinh}$','','$\psi_{\mathrm{cKB}}$',''); 
xlabel('$m$'); title(['$\lambda=$ ',num2str(lambda(ind_lambda))])
colororder(["#EDB120";"#EDB120";"#D95319";"#D95319";"#FE5BAC";"#7E2F8E";"#77AC30";"#77AC30";"#00008B";"#00008B"])
end%for
sgtitle({'Figure 4.25: Maximum approximation error (solid) and error constants (dashed) using ','regularizations (4.36) in frequency domain compared to regularizations (4.163) in spatial domain with selected window functions each, ',' for the function (4.151) where $M=256$, $m\in\{2, 3, \ldots, 10\}$, and $\lambda\in\{0.5,1,2\}$.'});

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('comparison_windows.txt','w');
fprintf(fileID,['\n\n Test for different lambda with N=',num2str(N)]);
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
format = '%d %1.4e \n';

% Classical Shannon sampling sums
fprintf(fileID,'Error for classical Shannon sampling sums');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_shannon(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for classical Shannon sampling sums');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_shannon(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Linear frequency window function
fprintf(fileID,'Error for linear frequency window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_lin(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for linear frequency window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_lin(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Convolution based window functions
fprintf(fileID,'Error for convolution based frequency window function with m=2');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_conv2(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error for convolution based frequency window function with m=3');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_conv3(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Error for convolution based frequency window function with m=4');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_conv4(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error for convolution based frequency window function with m=inf');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_conv_inf(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Cubic frequency window function
fprintf(fileID,'Error for cubic frequency window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_cub(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for cubic frequency window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_cub(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Raised cosine frequency window function
fprintf(fileID,'Error for cosine frequency window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_cos(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Rectangular window function
fprintf(fileID,'Error for rectangular window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_rect(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for rectangular window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_rect(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Gaussian window function
fprintf(fileID,'Error for Gaussian window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_Gauss(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for Gaussian window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_Gauss(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% B-spline window function
fprintf(fileID,'Error for B-spline window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_B(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for B-spline window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_B(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% sinh-type window function
fprintf(fileID,'Error for sinh-type window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_sinh(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for sinh-type window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_sinh(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Continuous Kaiser-Bessel window function
fprintf(fileID,'Error for continuous Kaiser-Bessel window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_cKB(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for continuous Kaiser-Bessel window function');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',const_cKB(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
if (rat_flag == 1)
% Rational approximation by Natterer
fprintf(fileID,'Rational approximation by Natterer');
for ind_lambda = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(ind_lambda)),'\n\n']);
    matrix = [n.',err_rat(:,ind_lambda)];
    fprintf(fileID,format,matrix.');
end%for
end%if
fclose(fileID);
end%if

%% Function definitions

% Definition of the sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end%function


% Definition of the rational approximation by Natterer
function f = natterer(x)
f = zeros(size(x));

for j = 1:length(x)
    y = abs(x(j));
    a = zeros(1,4);
    b = zeros(1,4);
     
    if (y>=0 && y<5)
        m = 2.5;
        a = [0.584728 -0.251411 -0.010858 0.007110];
        b = [1 0.040010 0.022789 0.000517];
    elseif (y>=5 && y<10)
        m = 7.5;
        a = [-0.071488 0.038915 0.015432 -0.004603];
        b = [1 0.113894 0.034231 0.001948];
    elseif (y>=10 && y<15)
        m = 12.5;
        a = [0.002007 -0.020457 -0.000241 0.001651];
        b = [1 0.181669 0.057192 0.007547];
    elseif (y>=15 && y<20)
        m = 17.5;
        a = [0.008740 0.003174 -0.002948 -0.000050];
        b = [1 0.118630 0.095243 0];
    elseif (y>=20 && y<25)
        m = 22.5;
        a = [-0.002973 0.003391 0.000693 -0.000354];
        b = [1 0.030828 0.048077 -0.002685];
    elseif (y>=25 && y<30)
        m = 27.5;
        a = [-0.001414 -0.001887 0.000306 0.000207];
        b = [1 0.241197 0.078265 0.018285];
    elseif (y>=30 && y<35)
        m = 32.5;
        a = [0.001169 -0.000767 -0.000311 0.000117];
        b = [1 -0.071068 0.049563 -0.012647];
    elseif (y>=35 && y<=40)
        m = 37.5;
        a = [0.000257 0.000685 -0.000060 -0.000066];
        b = [1 0.154374 0.066485 0.010572];
    elseif y>40 
        m = 0;
        b = ones(1,4);
%         error('Only arguments in [0,40] allowed.'); 
    end%if
    
    f(j) = sum(a.*(y-m).^(0:3)) / sum(b.*(y-m).^(0:3));
end%for
end%function