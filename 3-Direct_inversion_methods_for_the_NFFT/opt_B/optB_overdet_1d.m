function [B,t,info] = optB_overdet_1d(N,M,m,x,window,n,mode)
% Optimization of the entries of B
% via   Min ||A'B-D^{-1}F'||_F
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
% B ... Optimized matrix
% t ... Run-time
% info. Additional information (optional)
% 
% Author: Melanie Kircheis


tic

% Default method for computation
if nargin==6
    mode = 'normal_equations';
end%if

% Suppress certain warnings that matrix may be singular
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')

% Log some additional information if desired
if nargout==3
    info_flag = 1;
else
    info_flag = 0;
end

if not(issorted(x))
    [x,ind_sort] = sort(x);
    sort_flag = 1;
else 
    sort_flag = 0;
end%if
x = x(:)';

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

B = spalloc(N,m,(2*n+1)*N);
k = (-M/2:M/2-1)';

d = phi_hat(k,M,m,window,n);

for l = -m/2:m/2-1
  
  % Determine indices of non-zero entries
  ind_l = mod((l-n:l+n-1)+m/2,m)+1;
  ind = [];
  for j = 1:length(ind_l)
      if not(aux(ind_l(j),1)==0)
          ind = [ind, aux(ind_l(j),1) : aux(ind_l(j),2)];
      end%if
  end%for
  
  if not(isempty(ind))
%         if rank(exp(-2*pi*1i*k*x(ind)))<length(ind) && ( strcmp(mode,'normal_equations') || strcmp(mode,'pcg') )
%             mode = 'naive';
%         end

  switch mode
      case 'naive'
          % Direct computation of bl
          fl = exp(-2i*pi*k*l/m);
          v = fl.*d;
          Hl = zeros(M, length(ind));
          for j = 1:length(ind)
            Hl(:,j) = exp(-2*pi*1i*k'*x(ind(j)));
          end
          b = Hl\v;
          
      case 'normal_equations'
          % Direct computation of bl via normal equations
          switch window 
              case 'Dirichlet'
                  v = zeros(length(ind),1);
                  z = x(ind)-l/m;
                  o = z==0;
                  v(o) = M;
                  v(~o) = sin((M-1)*pi*z(~o))./sin(pi*z(~o)) + exp(-M*pi*1i*z(~o));
              otherwise
                  v = exp(2i*pi*(x(ind)-l/m)'*k')*d;
          end
          y = x(ind)'-x(ind);
          HlHl = sin((M-1)*pi*y)./sin(pi*y) + exp(-M*pi*1i*y);
          HlHl(1:length(ind)+1:end) = M;
          b = HlHl\v;
          
      case 'lsqr'
          % Iterative computation of bl via LSQR with NFFT handle
          fl = exp(-2i*pi*k*l/m);
          plan = nfft(1,M,length(ind));
          plan.x = -x(ind)';  % NFFT defined with an additional minus
          [b,~] = lsqr(@(x,transp_flag) nfft_adj_mult(plan,x,transp_flag),(fl.*d),1e-10,100);
          
      case 'pcg'
          % Iterative computation of bl via PCG using NFFT handle
          fl = exp(-2i*pi*k*l/m);
          plan = nfft(1,M,length(ind));
          plan.x = -x(ind)';  % NFFT defined with an additional minus
          plan.fhat = fl.*d;
          nfft_trafo(plan);
          right = plan.f;
          [b,~] = pcg(@(x) my_fastsum(plan,x),right);
  end%switch
  
  % Get correct positions for nonzeros
  if sort_flag
      ind_x = ind_sort(ind);
  else 
      ind_x = ind;
  end%if

  % Compose the whole matrix
  B(:,l+m/2+1) = sparse(ind_x,1,b,N,1);

  if info_flag
    Hl = exp(-2*pi*1i*k*x(ind));
    Hl_pinv=pinv(Hl);
    T1=Hl_pinv*Hl; T2=Hl*Hl_pinv; 
    info.t1(l+m/2+1)=norm(T1-eye(length(ind)),'fro')^2;
    info.t2(l+m/2+1)=norm(T2-eye(prod(M)),'fro')^2;
    
    fl = exp(-2i*pi*k*l/m);
    v = fl.*d;
    info.t3(l+m/2+1) = norm(Hl*b-v,2)^2;
    info.d = d;
    info.condH(l+m/2+1) = cond(Hl);
    info.rankH(l+m/2+1) = rank(Hl)-length(ind);
    
    [Q,~]=qr(Hl);
    r=Q'*v; c = rank(Hl);
    [U,S,~]=svd(Hl);
    s=U'*v;
    info.t4(l+m/2+1) = norm(Hl*b-v,2)-norm(r(c+1:end),2)<5e-15;
    
  end%if
  end%if
end%for
t = toc;

end%function