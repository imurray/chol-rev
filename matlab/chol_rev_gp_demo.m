function chol_rev_gp_demo()

% Gaussian Process (GP) demo. Iain Murray, January 2016.

% The purpose here is to look at strategies for computing some derivatives, not
% to produce a full GP implementation. There are plenty of them. See
% http://www.gaussianprocess.org/gpml/ for more information on GPs.

% GP setup:
D = 20; % Dimension of input features
N = 2001; % Number of datapoints
ell = 1 + 0.1*randn(D, 1); % kernel lengthscales
sigma_n = 0.01; % fixed noise variance
xx = randn(D, N); % random input positions
% labels drawn from the prior:
yy = chol(plus_diag(kernel_fwd(xx, ell), sigma_n), 'lower')*randn(N, 1);

% Time three equivalent gradient computations, see the comments in each function
% below for more information on the different strategies.
tic; [F, d_ell] = f_ell_grad(ell, sigma_n, xx, yy); time_new_method = toc
tic; [F, d_ell2] = f_ell_grad_old_way(ell, sigma_n, xx, yy, 1); time_old_method = toc
tic; [F, d_ell3] = f_ell_grad_old_way(ell, sigma_n, xx, yy, 0); time_older_method = toc

fprintf('You can try to make "old methods" times better by compiling solve_chol from GPML.\n');
fprintf('"new method" times will be bad if you don''t compile the chol_rev mex file.\n');
fprintf('Having chol_rev.mex doesn''t make things *much* faster than the old way.\n');
fprintf('(In fact on some machines, the old way is faster.)\n');
fprintf('The new code is arguably neater and easier to derive though.\n');
fprintf('Pushing derivatives through reparameterizations defined using chol\n');
fprintf('will *require* chol_rev, which is the real motivation.\n');


function [F, ell_bar] = f_ell_grad(ell, sigma_n, xx, yy)
% Compute F, the -ve log probability of labels yy under a Gaussian determined by
% the other arguments. Also return derivatives ell_bar = dF/d_ell, found by
% back-propagating through the computation involving the Cholesky decomposition.
N = size(xx, 2);
[K, obj] = kernel_fwd(xx, ell);
L = chol(plus_diag(K, sigma_n), 'lower');
ww = L\yy;
F = 0.5*ww'*ww + sum(log(diag(L)));
% Backpropagation starts here:
L_bar = -tril((L'\ww)*ww');
L_bar = plus_diag(L_bar, 1./diag(L));
K_bar = chol_rev(L, L_bar);
K_bar = 0.5*(K_bar + K_bar'); % CARE: kernel code assumes full K matters
ell_bar = kernel_rev(obj, K_bar);


function [F, ell_bar] = f_ell_grad_old_way(ell, sigma_n, xx, yy, backprop)
% Compute F, the -ve log probability of yy under Gaussian determined by the
% other arguments and derivatives ell_bar = dF/d_ell. Traditionally the
% derivatives are performed on the expression involving matrix inverses, and
% then all the resulting expressions are evaluated using Cholesky decomposition
% solving where possible. We still end up inverting matrices though.
N = size(xx, 2);
[K, obj] = kernel_fwd(xx, ell);
L = chol(plus_diag(K, sigma_n), 'lower');
if exist('solve_chol') == 3
    solve_chol_ = @solve_chol;
else
    solve_chol_ = @local_solve_chol;
end
invM = solve_chol_(L', eye(N)); % we're inverting a matrix here
invM_y = solve_chol_(L', yy); % or just invM*yy ...
F = 0.5*yy'*invM_y + sum(log(diag(L)));
K_bar = 0.5*(invM - (invM_y*invM_y'));
if (nargin>4) && backprop
    % The back-propagation way (recommended):
    ell_bar = kernel_rev(obj, K_bar);
else
    % A lot of GP code swaps to computing forward derivatives for the kernel,
    % rather than continuing to back-propagate. That leads to a loop over each
    % lengthscale as here, or creating a higher-dimensional array containing
    % dK/d_ell for each lengthscale. Backpropagation seems easier and faster.
    D = numel(ell);
    ell_bar = zeros(size(ell));
    sq_dist = @(x) bsxfun(@minus, x', x).^2;
    for dd = 1:D
        dk_dl = 2*(K.*sq_dist(xx(dd,:)/(ell(dd).^1.5)));
        ell_bar(dd) = K_bar(:)'*dk_dl(:);
    end
end


function [K, obj] = kernel_fwd(xx, ell)
% Compute kernel matrix from input locations and lengthscales,
% using the standard RBF/squared-exponential/Gaussian kernel.
% This is the "forwards" computation.
% Inputs:
%       xx DxN N d-dimensional feature vectors of input locations
%      ell Dx1 lengthscales for the kernel
%
% Outputs:
%        K NxN Kernel matrix
%      obj cached computations for backpropagation
zz = bsxfun(@rdivide, xx, ell(:)); % DxN
z_2 = sum(zz.*zz, 1);
dd = bsxfun(@plus, bsxfun(@minus, z_2', zz'*(2*zz)), z_2);
K = exp(-dd);
if nargout > 1
    obj = struct(); obj.K = K; obj.zz = zz; obj.ell = ell(:);
end


function ell_bar = kernel_rev(obj, K_bar)
% This is the "reverse-mode" computation of derivatives corresponding to
% KERNEL_FWD. Compute derivatives wrt lengthscales from derivatives wrt
% kernel elements.
% Inputs:
%          obj cached computations from kernel_fwd
%        K_bar NxN df/dK for some function f 
%
% Outputs:
%      ell_bar Dx1 df/d_ell for the same function.
K = obj.K; zz = obj.zz; ell = obj.ell;
d_bar = -K_bar .* K;
z_bar = 4*bsxfun(@minus, bsxfun(@times, sum(d_bar, 1), zz), zz*d_bar);
ell_bar = - sum(z_bar.*bsxfun(@rdivide, zz, ell), 2);


function K = plus_diag(K, xx)
% Utility function: add xx onto diagonal of square matrix K
N = size(K, 1);
K(1:(N+1):end) = K(1:(N+1):end) + xx(:)';

function X = local_solve_chol(L, B)
% Ideally get the mex version from the GPML package
X = L\(L'\B);

