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

  const int n = X.rows();
  const int p = X.cols();
  const int S = Y.cols();
  const int n_theta = phi.size();
  const int q = basis.dim(1);

  vector<Type> theta(n_theta);
  for (int j = 0; j < n_theta; ++j) {
    theta(j) = (is_log(j) == 1) ? exp(phi(j)) : phi(j);
  }

  matrix<Type> LambdaZt(q, n);
  LambdaZt.setZero();
  for (int j = 0; j < n_theta; ++j) {
    for (int r = 0; r < q; ++r) {
      for (int c = 0; c < n; ++c) {
        LambdaZt(r, c) += theta(j) * basis(j, r, c);
      }
    }
  }

  matrix<Type> A = LambdaZt * LambdaZt.transpose();
  for (int i = 0; i < q; ++i) {
    A(i, i) += Type(1.0);
  }

  matrix<Type> L = lchol(A);
  if (L.rows() == 0) {
    return Type(1e20);
  }
  Type ldL2 = Type(0.0);
  for (int i = 0; i < q; ++i) {
    ldL2 += log(L(i, i));
  }
  ldL2 *= Type(2.0);

  matrix<Type> LambdaZtX = LambdaZt * X;
  matrix<Type> RZX = solve_lower(L, LambdaZtX);

  matrix<Type> M = X.transpose() * X - RZX.transpose() * RZX;
  matrix<Type> LX = lchol(M);
  if (LX.rows() == 0) {
    return Type(1e20);
  }
  Type ldRX2 = Type(0.0);
  for (int i = 0; i < p; ++i) {
    ldRX2 += log(LX(i, i));
  }
  ldRX2 *= Type(2.0);

  matrix<Type> LambdaZtY = LambdaZt * Y;
  matrix<Type> RZY = solve_lower(L, LambdaZtY);
  matrix<Type> Cu_Y = X.transpose() * Y - RZX.transpose() * RZY;
  matrix<Type> Cbeta = solve_lower(LX, Cu_Y);

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

  Type mean_log_PWRSS = Type(0.0);
  for (int s = 0; s < S; ++s) {
    mean_log_PWRSS += log(PWRSS_s(s) / np);
  }
  mean_log_PWRSS /= Type(S);

  if (!valid) {
    return Type(1e20);
  }

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
