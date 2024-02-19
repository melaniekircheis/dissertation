% Code file for Figures 4.16 -- 4.23

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0; 

% Switch flag for choice of spatial window function ('Gauss','bspline','sinh','cKB')
switch_func = 'Gauss';

% Set dimension
d = 1;

% Set truncation and oversampling parameters
n = 2:10;
lambda = [0;0.5;1;2];
sigma = 1+lambda;

% Initialization of error vectors
const = zeros(length(n),length(lambda)); 
err = zeros(length(n),length(lambda)); 
const_pert = zeros(length(n),length(lambda)); 
err_pert = zeros(length(n),length(lambda)); 
cond_true = zeros(length(n),length(lambda));

% Set bandwidth
a = 6/d;
N = 2^a;
tau = 0.5;
delta = tau*N;

% Set test function
% f = @(x) N^(d/2)*prod(my_sinc(N*pi,x),2);
f = @(x) (3*N/4)^(d/2)*prod((my_sinc(N*pi/2,x)).^2,2);
% f = @(x) (4*N/5)^(d/2)*prod((my_sinc(N*pi,x)+my_sinc(N*pi,(x-1))/2),2);

% Set maximum pertubation parameter
epsilon = 1e-3;
maxiter = 10;

%% Error constants

switch switch_func
    case 'Gauss'
        % Computation of error constant for Gaussian regularized sinc
        for k = 1:length(lambda)
            const(:,k) = (d*sqrt(2*n*sigma(k)*N^d)+sqrt(2*d*(sigma(k)*N)^d*lambda(k)*(1+n)))./(pi*n*sqrt(lambda(k))).*exp(-n*pi*lambda(k)/(2*sigma(k)));
            const_pert(:,k) = epsilon*(2+sqrt((2*n.*(sigma(k)))./(lambda(k)))).^d;
        end
    case 'bspline'
        % Computation of error constant for B-Spline regularized sinc
        for k = 1:length(lambda)
            b = ceil((n+1)/2);
            const(:,k) = (3*d*sqrt(b*N^d)./(sqrt(2)*pi*(2*b-1))).*exp(-n.*(log(pi*n*lambda(k))-log(2*b*sigma(k))));
            const_pert(:,k) = epsilon*(2+3/2.*sqrt(n)).^d;
            cond_true(:,k) = (lambda(k)>(n+2)./(pi*n-(n+2)));
        end
    case 'sinh'
        % Computation of error constant for sinh-type regularized sinc
        for k = 1:length(lambda)
            beta = pi*n.*(lambda(k))./(sigma(k));
            const(:,k) = (2^d-1)*N^(d/2).*exp(-beta);
            const_pert(:,k) = epsilon*(2+sqrt((2+2*lambda(k))./(lambda(k))).*sqrt(n)./(1-exp(-2*beta))).^d;
        end%for
    case 'cKB'
        % Computation of error constant for cKB regularized sinc
        for k = 1:length(lambda)
            beta = pi*n.*(lambda(k))./(sigma(k));
            const(:,k) = (2^d-1)*N^(d/2).*(7*pi*n*lambda(k).*(1+lambda(k)+4*n*lambda(k)))./(4*(sigma(k))^2).*exp(-beta);
            const_pert(:,k) = epsilon*(2+sqrt((2*n.*(sigma(k)))./(lambda(k)))).^d;
            cond_true(:,k) = (lambda(k)>1./(n-1));
        end%for
end%switch

%% Reconstruction error

% Initialize fine grid for evaluation of reconstruction error
S = 2^(12/d);
s = (-S:S)';
if d==2
    [S1,S2] = ndgrid(s);
    s = [S1(:) S2(:)];
elseif d==3
    [S1,S2,S3] = ndgrid(s);
    s = [S1(:) S2(:) S3(:)];
end%if
t = s/S;

% Initialization of vectors
Rm = zeros(length(t),1);
err_tilde = zeros(length(t),1);

% Loop for computation of the error
% Set function evaluations
ft = f(t);

for i2 = 1:length(n) 
    % Set truncation parameter
    m = n(i2);
        
    for k = 1:length(sigma)
        % Set oversampling
        L = sigma(k)*N;
        j = (-m-L:m+L)'; % Corresponding index set
        if d==2
            [J1,J2] = ndgrid(j);
            j = [J1(:) J2(:)];
        elseif d==3
            [J1,J2,J3] = ndgrid(j);
            j = [J1(:) J2(:) J3(:)];
        end%if

        % Set function evaluations
        fj = f(j/L);
        eps_j = (rand(length(j),maxiter)*2-1)*epsilon;
        fj_tilde = fj+eps_j;

        % Setup
        for i3 = 1:length(t)
            x = t(i3,:)-j/L;
            phi = my_sinc(L*pi,x); % Sinc 
            ind_delta = (abs(x)-m/L<=eps); % Characteristic function

            switch switch_func
                case 'Gauss'
                    % Evaluation of Gaussian regularized sinc function
                    mu = sqrt((m)./(pi*L*(L-2*delta))); % Set variance of Gaussian
                    if lambda(k)==0
                        psi = phi;
                    else
                        psi = phi.*exp(-x.^2/(2*mu.^2)); 
                    end
                    psi(~ind_delta) = 0;
                case 'bspline'
                    % Evaluation of B-Spline regularized sinc function
                    b = ceil((m+1)/2);
                    psi = phi./cardinal_bspline(0,2*b).*cardinal_bspline(b*L*x/m,2*b); psi(~ind_delta) = 0;
                case 'sinh'
                    % Evaluation of modified sinh-type regularized sinc function
                    beta = m*pi*lambda(k)./(1+lambda(k));
                    if beta==0
                        psi = phi.*sqrt(1-(L*x/m).^2);
                    else
                        psi = phi.*sinh(beta*sqrt(1-(L*x/m).^2))/sinh(beta); 
                    end
                    psi(~ind_delta) = 0;
                    cond_true(i2,k) = (lambda(k)>=log(2)/(pi*m-log(2)));
                case 'cKB'
                    % Evaluation of continuous Kaiser-Bessel regularized sinc function
                    beta = m*pi*lambda(k)./(1+lambda(k));
                    if beta==0
                        psi = phi.*(1-(L*x/m).^2);
                    else
                        psi = phi.*(besseli(0,beta*sqrt(1-(L*x/m).^2))-1)/(besseli(0,beta)-1);
                    end
                    psi(~ind_delta) = 0;
                case 'KB'
                    % Evaluation of Kaiser-Bessel regularized sinc function
                    beta = m*pi*lambda(k)./(1+lambda(k));
                    psi = phi.*besseli(0,beta*sqrt(1-(L*x/m).^2))/besseli(0,beta); psi(~ind_delta) = 0;
            end%switch
            psi = prod(psi,2);

            % Evaluation of regularized WKS sums
            Rm(i3) = psi.'*fj;
            Rm_tilde = psi.'*fj_tilde;
            err_tilde(i3) = norm(Rm(i3)-Rm_tilde,inf);
        end%for

        % Computation of reconstruction errors
        err(i2,k) = norm(Rm-ft,inf);
        err_pert(i2,k) = max(err_tilde);
    end%for
    
    fprintf(['m=',num2str(n(i2)),' done %s\n'], datestr(datetime('now')))
end%for

%% Visualization 

% Visualization of reconstruction error
ind_lambda = (const>1e-14); 
figure(1); semilogy(n,err(:,1),'-o',n(ind_lambda(:,2)),const(ind_lambda(:,2),2),'--',n,err(:,2),'-square',n(ind_lambda(:,3)),const(ind_lambda(:,3),3),'--',n,err(:,3),'-^',n(ind_lambda(:,4)),const(ind_lambda(:,4),4),'--',n,err(:,4),'-*'); 
xlabel('$m$'); legend('$e_{m,0}(f)$','$E_{m,0.5}$','$e_{m,0.5}(f)$','$E_{m,1}$','$e_{m,1}(f)$','$E_{m,2}$','$e_{m,2}(f)$','Location','northeastoutside'); 
colororder(["#111e99";"#FF007F";"#FF007F";"#D95319";"#D95319";"#008080";"#008080"])
switch switch_func
    case 'Gauss'
        title({'Figure 4.16: Maximum approximation error (4.152) and error constant (4.153)',' using  the continuous         Gaussian window function $\varphi_{\mathrm{Gauss}}$ in (4.60) with $\alpha = \frac{1}{M}\,\sqrt{\frac{m}{\pi(1+\lambda)\lambda}}$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
    case 'bspline'
        title({'Figure 4.18: Maximum approximation error (4.152) and error constant (4.157)',' using  the modified $\mathrm B$--spline window function $\varphi_{\mathrm{B}}$ in (4.61) with $s = \lceil\frac{m+1}{2}\rceil$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
    case 'sinh'
        title({'Figure 4.20: Maximum approximation error (4.152) and error constant (4.159)',' using  the $\sinh$-type window function $\varphi_{\mathrm{sinh}}$ in (4.62) with $\beta = \frac{\pi m \lambda}{1+\lambda}$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
    case 'cKB'
        title({'Figure 4.22: Maximum approximation error (4.152) and error constant (4.161)',' using  the continuous Kaiser--Bessel window function $\varphi_{\mathrm{cKB}}$ in (4.63) with $\beta = \frac{\pi m \lambda}{1+\lambda}$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
end%switch

% Visualization of robustness error
ind_lambda = (const_pert>1e-14); 
figure(2); semilogy(n,err_pert(:,1),'-o',n(ind_lambda(:,2)),const_pert(ind_lambda(:,2),2),'--',n,err_pert(:,2),'-square',n(ind_lambda(:,3)),const_pert(ind_lambda(:,3),3),'--',n,err_pert(:,3),'-^',n(ind_lambda(:,4)),const_pert(ind_lambda(:,4),4),'--',n,err_pert(:,4),'-*'); 
xlabel('$m$'); legend('$\tilde e_{m,0}(f)$','$\tilde E_{m,0.5}$','$\tilde e_{m,0.5}(f)$','$\tilde E_{m,1}$','$\tilde e_{m,1}(f)$','$\tilde E_{m,2}$','$\tilde e_{m,2}(f)$','Location','northeastoutside'); 
colororder(["#111e99";"#FF007F";"#FF007F";"#D95319";"#D95319";"#008080";"#008080"])
switch switch_func
    case 'Gauss'
        title({'Figure 4.17: Maximum perturbation error (4.155) and error constant (4.156)',' using the Gaussian window function $\varphi_{\mathrm{Gauss}}$ in (4.60) with $\alpha = \frac{1}{M}\,\sqrt{\frac{m}{\pi(1+\lambda)\lambda}}$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $\varepsilon=10^{-3}$, $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
    case 'bspline'
        title({'Figure 4.19: Maximum perturbation error (4.155) and error constant (4.158)',' using the modified $\mathrm B$--spline window function $\varphi_{\mathrm{B}}$ in (4.61) with $s = \lceil\frac{m+1}{2}\rceil$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $\varepsilon=10^{-3}$, $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
    case 'sinh'
        title({'Figure 4.21: Maximum perturbation error (4.155) and error constant (4.160)',' using the $\sinh$-type window function $\varphi_{\mathrm{sinh}}$ in (4.62) with $\beta = \frac{\pi m \lambda}{1+\lambda}$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $\varepsilon=10^{-3}$, $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
    case 'cKB'
        title({'Figure 4.23: Maximum perturbation error (4.155) and error constant (4.162)',' using the continuous Kaiser--Bessel window function $\varphi_{\mathrm{cKB}}$ in (4.63) with $\beta = \frac{\pi m \lambda}{1+\lambda}$ ','for the function $f(x) = \big(\frac{3M}{4}\big)^{d/2} \,\mathrm{sinc}^2\big(\frac M2 \pi x\big)$ with fixed $M=2^{6/d}$ ','as well as $\varepsilon=10^{-3}$, $m\in\{2,3,\ldots,10\}$, and $\lambda\in\{0,0.5,1,2\}$.'});
end%switch

%% Generate tables for tikz

if (save_results == 1)
switch switch_func
    case 'Gauss'       
        fileID = fopen(['error_Gauss_',num2str(d),'d.txt'],'w');
    case 'sinh'
        fileID = fopen(['error_sinh_',num2str(d),'d.txt'],'w');
    case 'bspline'
        fileID = fopen(['error_bspline_',num2str(d),'d.txt'],'w');
    case 'cKB'
        fileID = fopen(['error_cKB_',num2str(d),'d.txt'],'w');
end%switch

format = '%d %1.4e \n';

% % Test for different lambda
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,['\n\n Reconstruction error for different lambda with N= ',num2str(N),', tau= ',num2str(tau),'\n']);
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Error constant
fprintf(fileID,'Error constant');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',const(:,k)];
    fprintf(fileID,format,matrix.');
end%for
% Error 
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Error');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',err(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,['\n\n Pertubation error for different lambda with N= ',num2str(N),', tau= ',num2str(tau),'\n']);
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Error constant
fprintf(fileID,'Error constant');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',const_pert(:,k)];
    fprintf(fileID,format,matrix.');
end%for
% Error 
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
fprintf(fileID,'Error');
for k = 1:length(lambda)
    fprintf(fileID,['\n\n lambda=',num2str(lambda(k)),'\n\n']);
    matrix = [n.',err_pert(:,k)];
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