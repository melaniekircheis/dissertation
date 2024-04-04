function [B,t,info] = optB_overdet_2d(N,M,m,x,window,n,mode)
% Optimization of the entries of B
% via   Min ||A'BFD-I||_F
% 
% INPUT:
% N ... Number of nodes
% M ... Numbers of Fourier coefficients
% m ... Oversampled M
% x ... Nodes
% window ... Specification of window function - 'BSpline', 'Gaussian', 'Sinc', 'Bessel', 'Kaiser', 'Dirichlet'
% n ... Cut-off parameter
% mode ... Specification of the method to compute the entries - 'naive', 'normal_equations', 'lsqr', 'pcg'
% 
% OUTPUT:
% B ... Optimized matrix
% t ... Run-time
% cond_H ... Condition number of used matrices (Note that this will cause
% the most amount of time)
% 
% Author: Melanie Kircheis


tic
m = m(:);

% Default method for computation
if nargin==6
    mode = 'normal_equations';
end

% Suppress certain warnings that matrix may be singular
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')

% Log condition number of used matrices if desired
% Note that this will cause the most amount of time
if nargout==3
    info_flag = 1;
    info.cond_H = zeros(prod(m),1);
else
    info_flag = 0;
end

% Initialize matrix B
B = spalloc(N,prod(m),(2*n+1)^2*N);

% Determine necessary index sets
k1 = (-M(1)/2:M(1)/2-1)';
[K1,K2] = ndgrid(k1,k1);
k = [K1(:) K2(:)];
l1 = -m(1)/2:m(1)/2-1;
[L1,L2] = ndgrid(l1,l1);
l = [L1(:) L2(:)];

% Precomputations for index set
[mx,ind_sort] = sortrows(floor(x.*m'),[2,1]);
aux = zeros(prod(m),2);
j = 1;
h = 1;
while j<=N && h<=prod(m)
    if mx(j,:)==l(h,:)
        if aux(h,1)==0
            aux(h,:) = j;
        else
            aux(h,2) = j;
        end%if
        j = j+1;
    else
        h = h+1;
    end%if
end%while

% Loop counter
h = 1;

% Position of loop counter on grid l
row = 1;
col = 1;

while h <= prod(m)
        
    % Determine all l in window
    rows = mod(row-n:row+n-1,m(1));
    rows(rows==0) = m(1);
    cols = mod(col-n:col+n-1,m(1));
    cols(cols==0) = m(1);
    ind_l = (rows-1)*m(1)+cols';
    ind_l = ind_l(:);
    
    % Determine indices of non-zero entries
    ind = [];
    for j = 1:length(ind_l)
        if not(aux(ind_l(j),1)==0)
            ind = [ind, aux(ind_l(j),1) : aux(ind_l(j),2)];
        end%if
    end%for
    
    % Get right positions
    ind_x = ind_sort(ind);
    
    if not(isempty(ind_x))
    
    switch mode
        case 'naive'
        % Direct computation of bl -- solve Hl*bl=v with backslash
            % Generate matrix Hl
            Hl = exp(-2*pi*1i*(k*x(ind_x,:).'));

            % Computation of Kronecker factor of diagonal entries of D
            phihat = phi_hat(k1,M(1),m(1),window,n);

            % Computation of Kronecker factors of columns of F
            ind1 = ceil(h/m(1));
            ind2 = mod(h,m(1));
            if ind2 == 0, ind2 = m(1); end
            f1 = exp(-2i*pi*k1*(l1(ind1)'./m(1)));
            f2 = exp(-2i*pi*k1*(l1(ind2)'./m(1)));

            % Computation of the right side
            v = kron(f1.*phihat,f2.*phihat);

            % Computation of bl
            b = Hl\v;

        case 'normal_equations'
        % Direct computation of bl via normal equations -- use of backslash
            switch window
                case 'Dirichlet'
                    v = zeros(length(ind_x),2);
                    z = x(ind_x,:)-l(h,:)./m';
                    o = z==0;
                    v(o) = M(1);
                    v(~o(:,1),1) = sin((M(1)-1)*pi*z(~o(:,1),1))./sin(pi*z(~o(:,1),1)) + exp(-M(1)*pi*1i*z(~o(:,1),1));
                    v(~o(:,2),2) = sin((M(2)-1)*pi*z(~o(:,2),2))./sin(pi*z(~o(:,2),2)) + exp(-M(1)*pi*1i*z(~o(:,2),2));
                    v = v(:,1).*v(:,2);
                otherwise
                    % Computation of Kronecker factor of diagonal entries of D
                    phihat = phi_hat(k1,M(1),m(1),window,n);
                    d = kron(phihat,phihat);

                    % Computation of the right side
                    v = exp(2i*pi*(x(ind_x,:)-l(h,:)./m')*k')*d;
            end%switch

            % Determine differences of nodes
            y = zeros(length(ind_x)*(length(ind_x)-1)/2,2);
            for j = 1:length(ind_x)-1
                y(length(ind_x)*(j-1)-j*(j-1)/2+1:length(ind_x)*j-j/2*(j+1),:) = x(ind_x(j+1:end),:) - x(ind_x(j),:);
            end%for
            
            % Determine matrix of normal equations via Dirichlet kernel
            vec1 = zeros(size(y,1),1);
            vec2 = zeros(size(y,1),1);
            o1 = y(:,1)==0;  % Index set
            o2 = y(:,2)==0;  % Index set
            vec1(o1) = M(1);
            vec1(~o1) = (sin((M(1)-1)*pi*y(~o1,1))./sin(pi*y(~o1,1)) + exp(-M(1)*pi*1i*y(~o1,1)));
            vec2(o2) = M(2);
            vec2(~o2) = (sin((M(2)-1)*pi*y(~o2,2))./sin(pi*y(~o2,2)) + exp(-M(2)*pi*1i*y(~o2,2)));
            HlHl = zeros(length(ind_x));
            % Compute only lower triangle
            ind = tril(true(length(ind_x)),-1);
            HlHl(ind) = vec1 .* vec2;
            % Determine remaining entries (matrix Hermitesch, diagonal known)
            HlHl = HlHl + HlHl' + prod(M)*eye(length(ind_x));

            % Computation of bl
            if info_flag; info.cond_H(h) = cond(HlHl); end%if
            b = HlHl\v;

        case 'lsqr'
        % Iterative computation of bl via LSQR with NFFT handle -- applied to Hl*bl=v
            % Computation of Kronecker factor of diagonal entries of D
            phihat = phi_hat(k1,M(1),m(1),window,n);

            % Computation of Kronecker factors of columns of F
            ind1 = ceil(h/m(1));
            ind2 = mod(h,m(1));
            if ind2 == 0, ind2 = m(1); end
            f1 = exp(-2i*pi*k1*(l1(ind1)'./m(1)));
            f2 = exp(-2i*pi*k1*(l1(ind2)'./m(1)));

            % Computation of the right side
            v = kron(f1.*phihat,f2.*phihat);

            plan = nfft(2,M,length(ind_x));
            plan.x = -x(ind_x,:);  % NFFT defined with an additional minus
            [b,~] = lsqr(@(x,transp_flag) nfft_adj_mult(plan,x,transp_flag),v,eps,100);

        case 'pcg'
        % Iterative computation of bl via PCG with NFFT handle -- applied to ???
            % Computation of Kronecker factor of diagonal entries of D
            phihat = phi_hat(k1,M(1),m(1),window,n);

            % Computation of Kronecker factors of columns of F
            ind1 = ceil(h/m(1));
            ind2 = mod(h,m(1));
            if ind2 == 0, ind2 = m(1); end
            f1 = exp(-2i*pi*k1*(l1(ind1)'./m(1)));
            f2 = exp(-2i*pi*k1*(l1(ind2)'./m(1)));

            % Computation of the right side
            v = kron(f1.*phihat,f2.*phihat);

            plan = nfft(2,M,length(ind_x));
            plan.x = -x(ind_x,:);  % NFFT defined with an additional minus
            plan.fhat = v;
            nfft_trafo(plan);
            right = plan.f;
            [b,~] = pcg(@(x) fastsum(plan,x),right,eps,100);
    end%switch
    
    if info_flag
        Hl = exp(-2*pi*1i*(k*x(ind_x,:).'));
        Hl_pinv=pinv(Hl);
        T1=Hl_pinv*Hl; T2=Hl*Hl_pinv; 
        info.t1(h)=norm(T1-eye(length(ind)),'fro')^2;
        info.t2(h)=norm(T2-eye(prod(M)),'fro')^2;
        
        % Computation of Kronecker factor of diagonal entries of D
        phihat = phi_hat(k1,M(1),m(1),window,n);

        % Computation of Kronecker factors of columns of F
        ind1 = ceil(h/m(1));
        ind2 = mod(h,m(1));
        if ind2 == 0, ind2 = m(1); end
        f1 = exp(-2i*pi*k1*(l1(ind1)'./m(1)));
        f2 = exp(-2i*pi*k1*(l1(ind2)'./m(1)));

        % Computation of the right side
        v = kron(f1.*phihat,f2.*phihat);

    	info.t3(h) = norm(Hl*b-v,2)^2;
        info.d = kron(phihat,phihat);
    end%if
        
    % Compose the whole matrix
    B(:,h) = sparse(ind_x,1,b,N,1);
    end%if
    
    % Increase loop counter
    h = h+1;
    
    % Adjust position of index h
    if col < m(1)
        col = col + 1;
    elseif row < m(1)
        col = 1;
        row = row + 1;
    end
end%while

% Computation time
t = toc;

end%function