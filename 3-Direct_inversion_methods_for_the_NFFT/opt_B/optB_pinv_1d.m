function [S,t,d,info] = optB_pinv_1d(N,M,m,x,window,n,r,mode)
% Optimization of the entries of pinv(B')
% via   Min ||B'S-1/m*I||_F
% 
% INPUT:
% N ... Number of nodes
% M ... Number of Fourier coefficients
% m ... Oversampled M
% x ... Nodes
% window ... Specification of window function - 'BSpline', 'Gaussian', 'Sinc', 'Bessel', 'Kaiser', 'Dirichlet'
% n ... Cut-off parameter
% mode ... Specification of the method to compute the entries - 'naive', 'normal_equations', 'lsqr', 'pcg'
% 
% OUTPUT:
% S ... Optimized matrix
% t ... Run-time
% d ... Optimized scaling factors in D (optional)
% info. Additional information on matrix of linear system of equations to solve (optional)
% 
% Author: Melanie Kircheis


tic

% Default method for computation
if ( (nargin<7) || (isempty(r)) )
    r = m/2;
end%if
if ( (nargin<8) || (isempty(mode)) )
    mode = 'normal_equations';
end%if

% Set switch for additional information
if nargout<4
    info_flag = 0;
else 
    info_flag = 1;
end%if

% Set switch for additional optimization of D
if nargout<3
    diag_flag = 0;
else 
    diag_flag = 1;
end%if

% Sort if necessary
if not(issorted(x))
    [x,ind_sort] = sort(x);
    sort_flag = 1;
else 
    sort_flag = 0;
end%if

% Precomputations for index set
mx = floor(m*x);
aux = zeros(m,2);
j = 1;
l = -m/2;
while j<=N && l<m/2
    if mx(j)==l
        if aux(l+m/2+1,1)==0
            aux(l+m/2+1,:) = j;
        else
            aux(l+m/2+1,2) = j;
        end%if
        j = j+1;
    else
        l = l+1;
    end%if
end%while

% Generate adjoint of matrix B
B_adj = sparse(prod(m),N);
l = -m/2:m/2-1;
for j = 1:N
  y = x(j)-l./m';
  y = mod(y+0.5,1)-0.5;
  phi1 = phi(y,M(1),m(1),window,n);
  B_adj = B_adj + sparse(1:prod(m),j,phi1,prod(m),N);
end

% Init
S = spalloc(N,m,(2*n+1)*N);

% Computation of the optimal entries of sparse matrix
for l = -m/2:m/2-1
  
  % Determine indices of non-zero entries
  ind_l = mod((l-n:l+n-1)+m/2,m)+1;
  ind = [];
  for j = 1:length(ind_l)
      if not(aux(ind_l(j),1)==0)
          ind = [ind, aux(ind_l(j),1) : aux(ind_l(j),2)];
      end%if
  end%for
  
  % Choose only points in an r-vinicity of l
  ind2 = unique(mod((l-r:l+r)+m/2,m)+1);
  
  % Set auxiliary matrix and right side
  Gl = B_adj(ind2,ind);
  v = zeros(m,1); v(l+m/2+1) = 1/m; v = v(ind2);

  % Computation of the nonzero elements
  if not(isempty(ind))
  switch mode
      case 'naive'
          % Direct computation of sl
          s = Gl\v;
          % Compute some additional information if desired
          if info_flag
              info.rank(l+m/2+1) = rank(full(Gl));
          end%if
      case 'lsqr'
          % Iterative computation of sl via LSQR
          [s,~] = lsqr(Gl,v);
          % Compute some additional information if desired
          if info_flag
              info.rank(l+m/2+1) = rank(full(Gl));
          end%if
      case 'normal_equations'
          % Direct computation of sl via normal equations
          GlGl = Gl'*Gl;
          s = (GlGl\(Gl'*v));
          % Compute some additional information if desired
          if info_flag
              info.cond(l+m/2+1) = cond(GlGl);
              info.rank(l+m/2+1) = rank(full(GlGl));
          end%if
      case 'pcg'
          % Iterative computation of sl via PCG
          GlGl = Gl'*Gl; 
          [s,~] = pcg(GlGl,Gl'*v);
          % Compute some additional information if desired
          if info_flag
              if not(min(eigs(GlGl,length(ind))>0)), warning('PCG is used for a matrix that is not positive definite!'); end
              info.cond(l+m/2+1) = condest(GlGl);
              info.rank(l+m/2+1) = rank(full(GlGl));
          end%if
  end%switch
  
  % Get correct positions for nonzeros
  if sort_flag
      ind_x = ind_sort(ind);
  else 
      ind_x = ind;
  end%if

  % Compose the whole matrix
  S(:,l+m/2+1) = sparse(ind_x,1,s,N,1);

  end%if
end%for


% Optional computation of optimized scaling factors in D
if diag_flag
    % Generate matrix A
    A = zeros(N,M);
    j = 1:N;
    for k = -M/2 : M/2-1
        A(:,k + M/2 + 1) = exp(2*pi*1i*k*x(j));
    end

    % Generate matrix F
    F = zeros(m,M);
    l = -m/2 : m/2-1;
    for k = -M/2 : M/2-1
        F(:,k + M/2 + 1) = exp(2*pi*1i*k*l/m);
    end

    % Computation of optimal matrix D
    d = zeros(M,1); 
    % Precompute necessary matrix
    H = A'*S;
    for s = -M/2:M/2-1
        % Computation of corresponding vector
        h = H*(F(:,s+M/2+1));
        % Computation of the enumerator
        enum = h(s+M/2+1);
        % Computation of the denominator
        denom = h'*h;
        % Computation of the optimal scaling factors
        d(s+M/2+1) = enum./denom;
    end%for
end%if

t = toc;

end%function