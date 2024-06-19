classdef fastsinc < handle   
% This class provides a Matlab interface for the sinc-transform, i.e., the fast evaluation of sinc sums.
%
% Usage example:
% plan = fastsinc(d,N);
% plan.a = a;
% plan.c = c;
% fastsinc_precompute_least_squares(plan) %or fastsinc_precompute_analytic(plan);
% plan.b = b;
% fastsinc_trafo(plan);
% result = plan.result;

    properties
        d       % Dimension
        N       % Nonharmonic bandwith (natural number)
        w       % Quadrature weights
        z       % Quadrature points (optional input)
        c       % Coefficients for sinc-trafo
        a       % Points a_k in sinc-trafo
        b       % Points b_l in sinc-trafo
        result  % Computed sinc transform
        time    % Computation times
        precomputations_done = false;   % Flag if precomputations were done
        flag_a_equi                     % Flag if the points a_k are equispaced
        flag_b_equi                     % Flag if the points b_l are equispaced
        flag_a_tensorizable             % Flag if the points a_k and the coefficients c_k can be tensorized
        flag_c_factors_equal            % Flag if the factors of the coefficients c_k are equal (if tensorizable)
        flag_b_tensorizable             % Flag if the points b_l can be tensorized
    end%properties
    
    properties(Hidden=true,SetAccess='private',GetAccess='private')
        N_star                          % Scaled bandwidth (only needed for NNFFT)
        n                               % Number of quadrature points (optional input)
        y                               % Evaluation points for least squares precomputations (optional input)
        P                               % Number of evalution points for least squares precomputations (optional input)
        flag_y_equi                     % Flag if the evaluation points y for LSQR are equispaced
        tol                             % Tolerance for LSQR (optional input)
        maxit                           % Maximum iteration number for LSQR (optional input)
        a_is_set = false;               % Flag if the points a_k are set
        c_is_set = false;               % Flag if the coefficients c_k are set
        b_is_set = false;               % Flag if the points b_l are set
    end%properties
    
    
    methods
    % Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.a(h,a)
            h.a_is_set = false;
            if( isempty(a) || ~isnumeric(a) || ~isreal(a) || all(sum(abs(a)>1/2)) )
                error('The points a_k have to lie in the interval [-1/2,1/2)^d.');
            % Save the points a_k as a column vector
            elseif size(a,1) > size(a,2)
                % Check whether the tensor decomposition is given
                if size(a,2) == 1
                    h.flag_a_tensorizable = true;
                    h.flag_a_equi = all(diff(a,2)<eps) && (h.N==length(a)); % Check equispacing
                elseif size(a,2) == h.d
                    h.flag_a_tensorizable = false;
                else
                    error('The points a_k have to be given in their d-dimensional or tensorized form.');
                end%if
                h.a = a;
            elseif size(a,1) < size(a,2)
                % Check whether the tensor decomposition is given
                if size(a,1) == 1
                    h.flag_a_tensorizable = true;
                    h.flag_a_equi = all(diff(a,2)<eps) && (h.N==length(a)); % Check equispacing
                elseif size(a,1) == h.d
                    h.flag_a_tensorizable = false;
                else
                    error('The points a_k have to be given in their d-dimensional or tensorized form.');
                end%if
                h.a = a.';
            else
                warning('Cannot determine whether the points a_k are given as a column vector. This may cause unexpected behavior.');
                h.a = a;
            end%if
                h.a_is_set = true;
        end%function
        
        function set.c(h,c)
            h.c_is_set = false;
            if ~h.a_is_set
                error('The coefficients c_k can only be set after the points a_k were given.');
            % Save the coefficients c_k as a column vector
            elseif size(c,1) == length(h.a)              
                % Check the tensor decomposition
                if size(c,2) == 1
                    h.flag_c_factors_equal = h.flag_a_tensorizable;
                elseif size(c,2) == h.d
                    if h.flag_a_tensorizable == true
                        h.flag_c_factors_equal = false;
                    else
                        error('The coefficients c_k have to be 1-dimensional.');
                    end%if
                else
                    error('The coefficients c_k have to be given as a 1-dimensional vector or in their d-dimensional tensorized form.');
                end%if
                h.c = c;
            elseif size(c,2) == length(h.a)              
                % Check the tensor decomposition
                if size(c,1) == 1
                    h.flag_c_factors_equal = h.flag_a_tensorizable;
                elseif size(c,1) == h.d
                    if h.flag_a_tensorizable == true
                        h.flag_c_factors_equal = false;
                    else
                        error('The coefficients c_k have to be 1-dimensional.');
                    end%if
                else
                    error('The coefficients c_k have to be given as a 1-dimensional vector or in their d-dimensional tensorized form.');
                end%if
                h.c = c.';
            else
                error('The number of coefficients c_k has to equal the number of points a_k.');
            end %if
            h.c_is_set = true;
        end %function
        
        function set.b(h,b)
            h.b_is_set = false;
            if( isempty(b) || ~isnumeric(b) || ~isreal(b) || all(sum(abs(b)>1/2)) )
                error('The points b_l have to lie in the interval [-1/2,1/2)^d.');
            % Save the points b_l as a column vector
            elseif size(b,1) > size(b,2)
                % Check whether the tensor decomposition is given
                if size(b,2) == 1
                    h.flag_b_tensorizable = true;
                    h.flag_b_equi = all(diff(b,2)<eps) && (h.N==length(b)); % Check equispacing
                elseif size(b,2) == h.d
                    h.flag_b_tensorizable = false;
                else
                    error('The points b_l have to be given in their d-dimensional or tensorized form.');
                end%if
                h.b = b;
            elseif size(b,1) < size(b,2)
                % Check whether the tensor decomposition is given
                if size(b,1) == 1
                    h.flag_b_tensorizable = true;
                    h.flag_b_equi = all(diff(b,2)<eps) && (h.N==length(b)); % Check equispacing
                elseif size(b,1) == h.d
                    h.flag_b_tensorizable = false;
                else
                    error('The points b_l have to be given in their d-dimensional or tensorized form.');
                end%if
                h.b = b.';
            else
                warning('Cannot determine whether the points b_l are given as a column vector. This may cause unexpected behavior.');
                h.b = b;
            end%if
                h.b_is_set = true;
        end%function

        
    % Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = fastsinc(d,N,varargin)
            % Input check for nonharmonic bandwidth
            if( isempty(d) || length(d)~=1 || ~isnumeric(d) || ~isreal(d) || mod(d,1)~=0 || d<=0 )
                error('The dimension d has to be a natural number.');
            else
                h.d = d;
            end%if
            if( isempty(N) || length(N)~=1 || ~isnumeric(N) || ~isreal(N) || mod(N,1)~=0 || N<=0 )
                error('The bandwidth N has to be a natural number.');
            else
                h.N = N;
                h.N_star = ceil((1+2*8/h.N)*h.N);
            end%if
            
            % Parse for optional input
            params = inputParser;
            params.KeepUnmatched = true;
            
            % Set the number of quadrature points and weights
            addParameter(params,'n',4*h.N, @(x) assert(mod(x,1)==0 && x>=2*h.N,'n has to be natural numbers bigger than 2N.'))
            parse(params,varargin{:});
            h.n = params.Results.n;   
            
            % Set additional parameters for least squares precomputation
            addParameter(params,'P',2.5*h.n, @(x) assert(mod(x,2)==0 && x>0,'P has to be an even natural number.'))
            addParameter(params,'tol',1e-11, @(x) assert(x>=0 && x<1 && not(ischar(x)),'tol has to be a real number smaller than 1.'))
            addParameter(params,'maxit',1000, @(x) assert( mod(x,1)==0 && x>0,'maxit must be natural number.'))
            parse(params,varargin{:});
            h.P = params.Results.P;
            h.tol = params.Results.tol;
            h.maxit = params.Results.maxit;
            
            % Set quadrature points
            addParameter(params,'z',1/2*(cos((0:h.n).'*pi/h.n)), @(x) assert(sum(abs(x)>1/2)==0 && length(x)==h.n+1 && (size(x,1)==1 || size(x,2)==1),'The quadrature points z in [-1/2,1/2] have to be a 1-dimensional vector of length n+1.'))
            parse(params,varargin{:});
            h.z = params.Results.z(:);
            
            % Set evaluation points for least squares precomputation
            addParameter(params,'y',[], @(x) assert(sum(abs(x)>1)==0 && length(x)==h.P+1 && (size(x,1)==1 || size(x,2)==1),'The valuation points y in [-1,1] have to be a 1-dimensional vector of length P+1.'))
            parse(params,varargin{:});
            h.y = params.Results.y;
            if isempty(h.y)
                h.y = (-1:2/h.P:1)';
                h.flag_y_equi = true;
            elseif all(diff(y,2)<eps)
                h.flag_y_equi = true;
            else 
                h.flag_y_equi = false;
            end%if
            
            % Read out flags for possible speedup
            addParameter(params,'flag_a_equi',false, @(x) assert(islogical(x),'Flags can only be logicals.'))
            addParameter(params,'flag_b_equi',false, @(x) assert(islogical(x),'Flags can only be logicals.'))
            addParameter(params,'flag_a_tensorizable',false, @(x) assert(islogical(x),'Flags can only be logicals.'))
            addParameter(params,'flag_c_factors_equal',false, @(x) assert(islogical(x),'Flags can only be logicals.'))
            addParameter(params,'flag_b_tensorizable',false, @(x) assert(islogical(x),'Flags can only be logicals.'))
            parse(params,varargin{:});
            h.flag_a_equi = params.Results.flag_a_equi;
            h.flag_b_equi = params.Results.flag_b_equi;
            h.flag_a_tensorizable = params.Results.flag_a_tensorizable;
            h.flag_c_factors_equal = params.Results.flag_c_factors_equal;
            h.flag_b_tensorizable = params.Results.flag_b_tensorizable;
        end%function
        
        
    % User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function for least squares precomputations
        function fastsinc_precompute_least_squares(h)
            % Evaluate the sinc function at given points y
            ev = sin(h.N*pi*h.y)./(h.N*pi*h.y); ev(h.y==0)=1; 
            % Iterative computation of weights using LSQR with function handles
            switch h.flag_y_equi
                case 1 % equispaced points y
                    % NFFT
                    tic; plan_nfft = nfft(1,h.P+2,h.n+1); const = 2*h.N/h.P; plan_nfft.x = -h.z*const;% NFFT init
                    [h.w,~] = lsqr(@(x,flag) nfft_adj_mult(plan_nfft,x,flag,1),ev,h.tol,h.maxit); 
                    h.time.t_precompute = toc;
                case 0 % nonequispaced points y
                    % NNFFT
                    tic; plan_nnfft = nnfft(1,h.n+2,h.P+1,2*h.N_star); % create plan of class type nfft
                    plan_nnfft.x = h.y/2; plan_nnfft.v = h.N/h.N_star*[0;h.z]; % set nodes in plan
                    plan_adjoint = nnfft(1,h.P+2,h.n+1,2*h.N_star); % create plan of class type nfft
                    plan_adjoint.x = -h.z; plan_adjoint.v = h.N/h.N_star*[0;h.y]/2; % set nodes in plan
                    nnfft_precompute_psi(plan_nnfft); nnfft_precompute_psi(plan_adjoint); % precomputations
                    [h.w,~] = lsqr(@(x,flag) nnfft_mult(plan_nnfft,plan_adjoint,x,flag,1,true),ev,h.tol,h.maxit);
                    h.time.t_precompute = toc;
            end%switch
            
            % Tensorize quadrature weights and points in fully multidimensional setting
            if ( h.d > 1 && ~h.flag_a_tensorizable )
                if ~h.a_is_set
                    error('The points a_k need to be set before doing the precomputation.');
                end%if
                w_1d = h.w; z_1d = h.z;
                for j = 1:h.d-1
                    h.w = kron(h.w,w_1d);
                    h.z = [kron(z_1d,ones((h.n+1)^j,1)), kron(ones(h.n+1,1),h.z)];
                end%for
            end%if

            h.precomputations_done = true;
        end%function
        
        % Function for analytic precomputations
        function fastsinc_precompute_analytic(h)
            % Precompute Clenshaw-Curtis weights
            tic; k = (0:h.n)';
            eps = ones(h.n+1,1); eps([1,end])=1/2;
            alpha = (2./(1-4*(0:h.n/2).^2)).';
            % Fast compuation of weights using FFT
            beta = zeros(h.n+1,1); beta(mod(k,2)==0) = alpha;
            h.w = [beta;flip(beta(2:end-1))];
            h.w = ifft(h.w);
            h.w = eps.*real(h.w(1:h.n+1)); h.time.t_precompute = toc;
            % Overwrite the quadrature points
            h.z = 1/2*(cos(k*pi/h.n));
            
            % Tensorize quadrature weights and points in fully multidimensional setting
            if ( h.d > 1 && ~h.flag_a_tensorizable )
                if ~h.a_is_set
                    error('The points a_k need to be set before doing the precomputation.');
                end%if
                w_1d = h.w; z_1d = h.z;
                for j = 1:h.d-1
                    h.w = kron(h.w,w_1d);
                    h.z = [kron(z_1d,ones((h.n+1)^j,1)), kron(ones(h.n+1,1),h.z)];
                end%for
            end%if

            h.precomputations_done = true;
        end%function
        

        % Fast sinc transform
        function result = fastsinc_trafo(h)
            % Check if all the necessary input was given
            if(~h.precomputations_done)
                error('Before doing a fast sinc transform you have to do precomputations.');
            elseif(~h.c_is_set)
                error('Before doing a fast sinc transform you have to set the coefficients c_k.');
            elseif(~h.a_is_set)
                error('Before doing a fast sinc transform you have to set the points a_k.');
            elseif(~h.b_is_set)
                error('Before doing a fast sinc transform you have to set the points b_l.');
            else
                % Fast evaluation of sinc sum by means of NNFFT and NFFT
                % step 1
                tic;
                switch h.flag_a_equi
                    case 1 % points a equispaced
                        % NFFT
                        plan_nfft = nfft(1,h.N,h.n+1);
                        plan_nfft.x = -h.z;
                        plan_nfft.fhat = h.c;
                        nfft_trafo(plan_nfft)
                        g = plan_nfft.f;
                    case 0 % points a nonequispaced
                        % NNFFT
                        plan_nnfft = nnfft(1,length(h.a),h.n+1,h.N_star); % create plan of class type nnfft
                        plan_nnfft.x = h.z; plan_nnfft.v = -h.N/h.N_star*h.a; % set nodes in plan
                        nnfft_precompute_psi(plan_nnfft); % precomputations
                        g = nnfft_mult([],plan_nnfft,h.c,'transp',0,true);
                end%switch

                % step 2 - Multiplication by precomputed weights
                g = g.*h.w;
                
                % step 3
            switch h.flag_b_equi
                case 1 % points b equispaced
                    % NFFT
                    plan_nfft = nfft(1,h.N,h.n+1);
                    plan_nfft.x = -h.z;
                    plan_nfft.f = g;
                    nfft_adjoint(plan_nfft)
                    result = plan_nfft.fhat;
                case 0 % points b nonequispaced
                    % NNFFT
                    plan_nnfft = nnfft(1,h.n+2,length(h.b),h.N_star); 
                    plan_nnfft.x = h.b; plan_nnfft.v = h.N/h.N_star*[0;h.z]; % set nodes in plan
                    nnfft_precompute_psi(plan_nnfft); % precomputations
                    result = nnfft_mult(plan_nnfft,[],g,'notransp',[0,1],true);
            end%switch
            h.result = result;
            h.time.t_trafo = toc;
            end%if
        end%function
        

    % Inner methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function handle for adjoint NFFT
        function x = nfft_adj_mult(plan,x,transp_flag,truncation,w)
          if nargin<4 || isempty(truncation)
              truncation = 0;
          end%if

          if strcmp(transp_flag,'transp')
            plan.fhat = [zeros(truncation,1);x];
            nfft_trafo(plan);
            if nargin < 5 
                x = plan.f;
            else
                x = w.*plan.f;
            end%if
          else
            if nargin < 5 
                plan.f = x;
            else
                plan.f = w.*x;
            end%if
            nfft_adjoint(plan);
            x = plan.fhat;
            if truncation>=1
                x = x(truncation+1:end);
            end%if
          end
        end%function
        
        
        % Function handle for NNFFT
        function f = nnfft_mult(plan_nnfft,plan_adjoint,fhat,transp_flag,truncation,split_flag)
              if nargin==3, transp_flag=[]; elseif nargin<3, error('Wrong number of input arguments.'); end
              if nargin<5 || isempty(truncation)
                truncation = 0;
              end%if
              if nargin<6 || isempty(split_flag)
                split_flag = false;
              end%if

              % NNFFT
              if not(strcmp(transp_flag,'transp'))
                plan_nnfft.fhat = [zeros(truncation,1);fhat];
                nnfft_trafo(plan_nnfft);
                f = plan_nnfft.f;
              else
                plan_adjoint.fhat = [zeros(truncation,1);fhat];
                nnfft_trafo(plan_adjoint);
                f = plan_adjoint.f;
                if split_flag
                    f(1) = conj(f(1));
                end%if
              end
        end%function
    end%methods
end%classdef