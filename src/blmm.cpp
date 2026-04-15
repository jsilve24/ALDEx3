// BLMM: Batched average profiled REML for ALDEx3
//
// Uses lme4's relative-Cholesky-factor parameterisation via precomputed
// basis matrices:  LambdaZt(theta) = sum_j  theta_j * basis[j,,]
// where basis[j,,] = d(LambdaZt)/d(theta_j)  (q x n, exact linear map).
//
// Unconstrained phi:
//   phi_j = log(theta_j)  if is_log(j) == 1  (bounded theta, e.g. diagonal Chol)
//   phi_j = theta_j        if is_log(j) == 0  (off-diagonal, already unconstrained)
//
// Objective: negative profiled average REML over n x S response block Y.
//   lbar = -1/2 [ldL2 + ldRX2 + (n-p) log(PWRSS_bar / (n-p))]
// where PWRSS_bar = mean_s PWRSS_s, computed via one q x q Cholesky + one
// multi-RHS solve (NOT S separate Cholesky factorisations).
//
// REPORT: PWRSS_s, L, RZX, LX, LambdaZt  (needed for score/fixed-effects)

#include <TMB.hpp>

// -----------------------------------------------------------------------
// Manual lower-Cholesky (compatible with CppAD reverse-mode AD)
// Returns L such that L L' = A.  Assumes A is symmetric positive definite.
// -----------------------------------------------------------------------
template<class Type>
matrix<Type> lchol(const matrix<Type>& A) {
  int n = A.rows();
  matrix<Type> L(n, n);
  L.setZero();
  for (int j = 0; j < n; j++) {
    Type s = A(j, j);
    for (int k = 0; k < j; k++) s -= L(j, k) * L(j, k);
    L(j, j) = sqrt(s);
    for (int i = j + 1; i < n; i++) {
      Type t = A(i, j);
      for (int k = 0; k < j; k++) t -= L(i, k) * L(j, k);
      L(i, j) = t / L(j, j);
    }
  }
  return L;
}

// Forward substitution: solve L * X = B  (L lower triangular, B is q x S)
template<class Type>
matrix<Type> fwd_sub(const matrix<Type>& L, const matrix<Type>& B) {
  int q = B.rows(), S = B.cols();
  matrix<Type> X(q, S);
  X.setZero();
  for (int s = 0; s < S; s++) {
    for (int i = 0; i < q; i++) {
      Type val = B(i, s);
      for (int k = 0; k < i; k++) val -= L(i, k) * X(k, s);
      X(i, s) = val / L(i, i);
    }
  }
  return X;
}

// -----------------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() () {

  // ----- Data -----------------------------------------------------------
  DATA_MATRIX(X);         // n x p  fixed-effects design
  DATA_MATRIX(Y);         // n x S  response block for one feature
  DATA_ARRAY(basis);      // n_theta x q x n  (dLambdaZt / dtheta_j)
  DATA_IVECTOR(is_log);   // 1 if phi_j = log(theta_j), else 0

  // ----- Parameters -----------------------------------------------------
  PARAMETER_VECTOR(phi);  // unconstrained, length = n_theta

  // ----- Dimensions -----------------------------------------------------
  int n       = X.rows();
  int p       = X.cols();
  int S       = Y.cols();
  int n_theta = phi.size();
  int q       = basis.dim(1);   // basis is [n_theta, q, n]

  // ----- phi -> theta ---------------------------------------------------
  vector<Type> theta(n_theta);
  for (int j = 0; j < n_theta; j++)
    theta(j) = (is_log(j) == 1) ? exp(phi(j)) : phi(j);

  // ----- LambdaZt = sum_j theta_j * basis[j,,]  (q x n) ----------------
  matrix<Type> LambdaZt(q, n);
  LambdaZt.setZero();
  for (int j = 0; j < n_theta; j++) {
    for (int r = 0; r < q; r++)
      for (int c = 0; c < n; c++)
        LambdaZt(r, c) += theta(j) * basis(j, r, c);
  }

  // ----- A = LambdaZt LambdaZt' + I_q  (q x q) -------------------------
  matrix<Type> A = LambdaZt * LambdaZt.transpose();
  for (int i = 0; i < q; i++) A(i, i) += Type(1.0);

  // ----- Cholesky of A: L L' = A ----------------------------------------
  matrix<Type> L = lchol(A);

  Type ldL2 = Type(0.0);
  for (int i = 0; i < q; i++) ldL2 += log(L(i, i));
  ldL2 *= Type(2.0);

  // ----- RZX = L^{-1} LambdaZt X  (q x p) ------------------------------
  matrix<Type> LambdaZtX = LambdaZt * X;           // q x p
  matrix<Type> RZX = fwd_sub(L, LambdaZtX);        // q x p

  // ----- M = X'X - RZX'RZX  (p x p) ------------------------------------
  matrix<Type> M = X.transpose() * X - RZX.transpose() * RZX;

  // ----- Cholesky of M: LX LX' = M  (p x p) ----------------------------
  matrix<Type> LX = lchol(M);

  Type ldRX2 = Type(0.0);
  for (int i = 0; i < p; i++) ldRX2 += log(LX(i, i));
  ldRX2 *= Type(2.0);

  // ----- Multi-RHS solve against Y_d ------------------------------------
  matrix<Type> LambdaZtY = LambdaZt * Y;                      // q x S
  matrix<Type> RZY        = fwd_sub(L,  LambdaZtY);           // q x S
  matrix<Type> Cu_Y       = X.transpose() * Y
                            - RZX.transpose() * RZY;           // p x S
  matrix<Type> Cbeta      = fwd_sub(LX, Cu_Y);                // p x S

  // ----- PWRSS_s = ||y_s||^2 - ||RZY_s||^2 - ||Cbeta_s||^2 -----------
  vector<Type> PWRSS_s(S);
  for (int s = 0; s < S; s++) {
    PWRSS_s(s) = Y.col(s).squaredNorm()
               - RZY.col(s).squaredNorm()
               - Cbeta.col(s).squaredNorm();
    if (PWRSS_s(s) < Type(0.0)) PWRSS_s(s) = Type(1e-10);
  }

  // ----- PWRSS_bar = mean(PWRSS_s) --------------------------------------
  Type PWRSS_bar = Type(0.0);
  for (int s = 0; s < S; s++) PWRSS_bar += PWRSS_s(s);
  PWRSS_bar /= Type(S);

  if (PWRSS_bar <= Type(0.0)) return Type(1e10);

  // ----- Profiled average REML (negative, for minimisation) -------------
  Type np  = Type(n - p);
  Type nll = Type(0.5) * (ldL2 + ldRX2 + np * log(PWRSS_bar / np));

  // ----- Report intermediate quantities for score & fixed-effects -------
  REPORT(PWRSS_s);
  REPORT(L);
  REPORT(RZX);
  REPORT(LX);
  REPORT(LambdaZt);
  REPORT(ldL2);
  REPORT(ldRX2);

  return nll;
}
