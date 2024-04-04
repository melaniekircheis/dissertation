function [fcheck,t] = nfft_adjoint_opt_2d(f,B,m,M,window,n)
% Performs an 2d adjoint NFFT with given sparse matrix B

tic

% Initialization of grid
k1 = -M(1)/2:M(1)/2-1;
k2 = -M(2)/2:M(2)/2-1;
[K1,K2] = ndgrid(k1,k2);
k = [K1(:) K2(:)];

% Multiplication with given matrix
fcheck = B'*f;

% Perform an FFT
fcheck = reshape(fcheck,m(1),m(2));
fcheck = fft2(fcheck);
fcheck = fcheck .* (transpose(exp(pi*1i*(0:m(1)-1)))*exp(pi*1i*(0:m(2)-1)));
fcheck = fftshift(fcheck);

% Multiplication with diagonal matrix
phihat_1 = phi_hat(k(:,1),M(1),m(1),window,n);
phihat_2 = phi_hat(k(:,2),M(2),m(2),window,n);
fcheck = fcheck(m(1)/2-M(1)/2+1:m(1)/2+M(1)/2,m(2)/2-M(2)/2+1:m(2)/2+M(2)/2);
fcheck = 1/prod(m)*1./(phihat_1.*phihat_2) .* reshape(fcheck,prod(M),1);

% Computation time
t = toc;

end%function