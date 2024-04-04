function [w,t] = dcf_pinv(x,M,M_sigma,kernel)
%function [w,t] = dcf_pinv(x,M,M_sigma,kernel)
% Function for computing optimal Density Compensation Factors, cf. [Sedarat et al.]
% 
% INPUT:
% x - Vector of nonequispaced nodes in spatial domain
% M - Number of Fourier coefficients
% kernel - specification of used convolution kernel ('sinc','Dirichlet','nfft')
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
elseif not(ismember(kernel, {'sinc','Dirichlet','nfft'}))
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

% Computation of the weights
% Compute enumerator
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
[U,~,~] = svd(C);  % Compute SVD
r = rank(C);
H = zeros(size(U));
H(1:r,1:r) = eye(r);
enum = diag(U*H*U');  % Compute C*inv(C'*C)*C' using SVD
% Compute denominator
denom = sum(C.^2,2);  % Compute CC' directly
% Set weights as the fraction
w = const*enum./denom;  

t = toc;
end%function