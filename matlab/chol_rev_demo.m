% Example setup:
N = 500;
A = cov(randn(3*N, N)); % a symmetric +ve definite matrix
L = chol(A)';           % lower-triangular Cholesky decomposition
Lbar = tril(randn(size(L))); % declare what dF/dL is for some function F

% Push the derivative back through the Cholesky, Abar = dF/d(tril(A))
Abar = chol_rev(L, Lbar);

% Perturb the input, look at the perturbation of the output and check
% for consistency.
dA = cov(randn(3*N, N));
eps = 1e-5;
dL = (chol(A + (eps/2)*dA)' - chol(A - (eps/2)*dA)') / eps;
dF1 = dL(:)' * Lbar(:);
dF2 = dA(:)' * Abar(:);

fprintf('Error: %g\n', dF1 - dF2);
