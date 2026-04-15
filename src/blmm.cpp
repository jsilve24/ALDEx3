// BLMM: batched profiled REML kernel for ALDEx3
//
// The implementation follows the profiled lme4-style Gaussian LMM
// parameterisation in terms of top-level covariance parameters theta and the
// profiled residual scale sigma^2. The anchor objective is genuinely batched:
// one factorisation of the shared q x q system per objective evaluation, one
// multi-right-hand-side solve against the draw block Y, and one average of the
// per-draw profiled REML residual terms across draws.

#include <TMB.hpp>

// TMB's matrix type is not directly compatible with Eigen LLT/TriangularView
// in this setting, so these AD-compatible kernels keep the implementation
// limited to the small q x q and p x p systems used by the profiled LMM.
template<class Type>
matrix<Type> lchol(const matrix<Type>& A) {
  const int n = A.rows();
  matrix<Type> L(n, n);
  L.setZero();
  for (int j = 0; j < n; ++j) {
    Type s = A(j, j);
    for (int k = 0; k < j; ++k) {
      s -= L(j, k) * L(j, k);
    }
    if (!(s > Type(0.0))) {
      matrix<Type> bad(0, 0);
      return bad;
    }
    L(j, j) = sqrt(s);
    for (int i = j + 1; i < n; ++i) {
      Type t = A(i, j);
      for (int k = 0; k < j; ++k) {
        t -= L(i, k) * L(j, k);
      }
      L(i, j) = t / L(j, j);
    }
  }
  return L;
}

template<class Type>
matrix<Type> solve_lower(const matrix<Type>& L, const matrix<Type>& B) {
  const int q = B.rows();
  const int nrhs = B.cols();
  matrix<Type> X(q, nrhs);
  X.setZero();
  for (int rhs = 0; rhs < nrhs; ++rhs) {
    for (int i = 0; i < q; ++i) {
      Type value = B(i, rhs);
      for (int k = 0; k < i; ++k) {
        value -= L(i, k) * X(k, rhs);
      }
      X(i, rhs) = value / L(i, i);
    }
  }
  return X;
}

template<class Type>
Type objective_function<Type>::operator()() {
  DATA_MATRIX(X);
  DATA_MATRIX(Y);
  DATA_ARRAY(basis);
  DATA_IVECTOR(is_log);

  PARAMETER_VECTOR(phi);

  // n = number of observations, p = fixed-effect columns,
  // S = number of MC draws (all sharing the same covariance parameters),
  // q = total dimension of the random-effect vector (sum of group-level blocks)
  const int n = X.rows();
  const int p = X.cols();
  const int S = Y.cols();
  const int n_theta = phi.size();
  const int q = basis.dim(1);

  // Unconstrain variance-component parameters: constrained entries (lower > -Inf)
  // were log-transformed in phi, so we exponentiate them back to get theta.
  vector<Type> theta(n_theta);
  for (int j = 0; j < n_theta; ++j) {
    theta(j) = (is_log(j) == 1) ? exp(phi(j)) : phi(j);
  }

  // LambdaZt (q x n) = Lambda(theta) %*% Zt, the scaled random-effects design
  // matrix. It is linear in theta, so we reconstruct it from the precomputed
  // basis slices: LambdaZt = sum_j theta_j * basis[j,,].
  matrix<Type> LambdaZt(q, n);
  LambdaZt.setZero();
  for (int j = 0; j < n_theta; ++j) {
    for (int r = 0; r < q; ++r) {
      for (int c = 0; c < n; ++c) {
        LambdaZt(r, c) += theta(j) * basis(j, r, c);
      }
    }
  }

  // A = I + LambdaZt %*% t(LambdaZt)  (q x q)
  // This is the inner covariance system that appears in the profiled LMM: under
  // the lme4 Cholesky parameterisation V = sigma^2 (I + Zt' Lambda' Lambda Zt),
  // A is the q x q block that needs to be factored for the profiled likelihood.
  matrix<Type> A = LambdaZt * LambdaZt.transpose();
  for (int i = 0; i < q; ++i) {
    A(i, i) += Type(1.0);
  }

  // L = lower Cholesky of A. ldL2 = log|det(A)| = 2 * sum(log(L_ii)).
  // This term enters the REML log-determinant correction.
  matrix<Type> L = lchol(A);
  if (L.rows() == 0) {
    return Type(1e20);  // A is not positive definite; reject this theta
  }
  Type ldL2 = Type(0.0);
  for (int i = 0; i < q; ++i) {
    ldL2 += log(L(i, i));
  }
  ldL2 *= Type(2.0);

  // RZX = L^{-1} %*% LambdaZt %*% X  (q x p)
  // M   = X'X - RZX'RZX               (p x p, the Schur complement of A in [A, LambdaZtX; ...])
  // LX  = lower Cholesky of M. ldRX2 = log|det(M)|, the fixed-effects REML
  // log-determinant term (equivalent to lme4's "ldRX2" in its notation).
  matrix<Type> LambdaZtX = LambdaZt * X;
  matrix<Type> RZX = solve_lower(L, LambdaZtX);

  matrix<Type> M = X.transpose() * X - RZX.transpose() * RZX;
  matrix<Type> LX = lchol(M);
  if (LX.rows() == 0) {
    return Type(1e20);  // Schur complement is singular; reject
  }
  Type ldRX2 = Type(0.0);
  for (int i = 0; i < p; ++i) {
    ldRX2 += log(LX(i, i));
  }
  ldRX2 *= Type(2.0);

  // RZY (q x S) and Cu_Y (p x S): same block-elimination applied to Y.
  // Cbeta = LX^{-1} %*% Cu_Y solves the profiled fixed-effects system for
  // all S draws simultaneously (multi-right-hand-side solve).
  matrix<Type> LambdaZtY = LambdaZt * Y;
  matrix<Type> RZY = solve_lower(L, LambdaZtY);
  matrix<Type> Cu_Y = X.transpose() * Y - RZX.transpose() * RZY;
  matrix<Type> Cbeta = solve_lower(LX, Cu_Y);

  // PWRSS_s[s] = penalised weighted residual sum of squares for draw s.
  // Using the block-elimination identity:
  //   ||y_s||^2 - ||RZY_s||^2 - ||Cbeta_s||^2
  // This equals y_s' P y_s where P is the hat-matrix complement, and it is
  // the numerator of the profiled sigma^2 estimate for draw s.
  const Type np = Type(n - p);
  vector<Type> PWRSS_s(S);
  bool valid = true;
  for (int s = 0; s < S; ++s) {
    PWRSS_s(s) = Y.col(s).squaredNorm()
      - RZY.col(s).squaredNorm()
      - Cbeta.col(s).squaredNorm();
    if (!(PWRSS_s(s) > Type(0.0))) {
      valid = false;
    }
  }

  if (!valid) {
    return Type(1e20);
  }

  // Profiled REML negative log-likelihood, averaged over all S draws.
  // Each draw contributes: 0.5 * (ldL2 + ldRX2 + (n-p) * log(PWRSS_s / (n-p)))
  // The log-determinant terms ldL2 and ldRX2 are shared across draws (they
  // depend only on theta), so averaging amounts to averaging only the PWRSS term.
  Type mean_log_PWRSS = Type(0.0);
  for (int s = 0; s < S; ++s) {
    mean_log_PWRSS += log(PWRSS_s(s) / np);
  }
  mean_log_PWRSS /= Type(S);

  const Type nll = Type(0.5) * (ldL2 + ldRX2 + np * mean_log_PWRSS);

  REPORT(PWRSS_s);
  REPORT(mean_log_PWRSS);
  REPORT(L);
  REPORT(RZX);
  REPORT(LX);
  REPORT(LambdaZt);
  REPORT(ldL2);
  REPORT(ldRX2);

  return nll;
}
