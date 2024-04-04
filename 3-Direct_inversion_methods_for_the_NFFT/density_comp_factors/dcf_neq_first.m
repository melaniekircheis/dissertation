function [w,t,opt] = dcf_neq_first(x,M,M_sigma,kernel,mode)
%function [w,t,opt] = dcf_neq_first(x,M,M_sigma,kernel,mode)
% Function for computing the Weighted density Compensation Factors, cf. [Sedarat et al.], [Rosenfeld]
% 
% INPUT:
% x - Vector of nonequispaced nodes in spatial domain
% M - Number of Fourier coefficients
% 
% OPTIONAL INPUT:
% kernel - Specification of used convolution kernel: 'sinc' (default),'Dirichlet','nfft'
% mode - Computation mode: 'lsqr','pcg' (default),'direct'         
% 
% OUTPUT:
% w - Computed DCF
% t - Computation time
% 
% Author: Melanie Kircheis
 
 
tic
 
% Input check - column vector needed
if size(x,1) < size(x,2)
    x = x.';
end%if

% Input check - validation of kernel input
if nargin<4 || isempty(kernel)
    kernel = 'sinc';
elseif not(ismember(kernel, {'sinc', 'Dirichlet','nfft'}))
    error('Invalid convolution kernel. Possible options are: ''sinc'', ''Dirichlet'' or ''nfft''.')
end%if

% Set equispaced points in spatial domain
if size(x,2) == 1
    k = (-M/2:M/2-1)';
    l = (-M_sigma/2:M_sigma/2-1)';
elseif size(x,2) == 2
    k = (-M/2:M/2-1)';
    [K1,K2] = ndgrid(k,k);
    k = [K1(:) K2(:)];
    l = (-M_sigma/2:M_sigma/2-1)';
    [L1,L2] = ndgrid(l,l);
    l = [L1(:) L2(:)];
end%if
 
% Set default for computation mode
if nargin < 5
    mode = 'pcg';
elseif not(ismember(mode, {'direct', 'lsqr', 'pcg'}))
    error('Invalid computation mode. Possible options are: ''direct'', ''lsqr'' or ''pcg''.')
end%if
 
% Computation of the weights
C = ones(); % Initialize convolution matrix
if strcmp(kernel,'sinc')
    for j=1:size(x,2)
        C = C .* (M*my_sinc(M*pi*(x(:,j)-l(:,j)'/M_sigma)));
    end%for
    const = M_sigma^size(x,2);
elseif strcmp(kernel,'Dirichlet')
    for j=1:size(x,2)
        C = C .* (M*my_sinc((M-1/2)*pi*(x(:,j)-l(:,j)'/M_sigma))./my_sinc(pi/2*(x(:,j)-l(:,j)'/M_sigma)));
    end%for
    const = M_sigma^size(x,2);
elseif strcmp(kernel,'nfft')
    C = exp(2*pi*1i*x*k');
    const = 1;
end%if
CC = C*C';
A = abs(CC).^2;  % Set up system matrix
b = const*diag(CC);  % Set up right side
% Set weights as solution to the system of equations
if strcmp(mode,'lsqr')
    [w,opt.flag,opt.relres,opt.iter]  = lsqr(A,b,1e-6,200);  
elseif strcmp(mode,'pcg')
    [w,opt.flag,opt.relres,opt.iter]  = pcg(A,b,1e-6,200);  
else
    w = A\b;
end%if
 
t = toc;
end%function