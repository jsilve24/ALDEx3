// BLMM: batched profiled REML kernel for ALDEx3
//
// The R layer caches design cross-products once. This kernel updates the sparse
// relative covariance factor from theta, constructs the shared random-effects
// system, and profiles fixed effects and residual scale over all draws.

#define TMB_LIB_INIT R_init_ALDEx3
#include <TMB.hpp>

// Small dense Cholesky used for the fixed-effect Schur complement and for the
// explicit dense backend. Returning a 0 x 0 matrix preserves the established
// sentinel behavior for a non-positive-definite system.
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
      Type z = A(i, j);
      for (int k = 0; k < j; ++k) {
        z -= L(i, k) * L(j, k);
      }
      L(i, j) = z / L(j, j);
    }
  }
  return L;
}

template<class Type>
matrix<Type> solve_lower(const matrix<Type>& L, const matrix<Type>& B) {
  matrix<Type> X(B.rows(), B.cols());
  X.setZero();
  for (int rhs = 0; rhs < B.cols(); ++rhs) {
    for (int i = 0; i < B.rows(); ++i) {
      Type z = B(i, rhs);
      for (int k = 0; k < i; ++k) {
        z -= L(i, k) * X(k, rhs);
      }
      X(i, rhs) = z / L(i, i);
    }
  }
  return X;
}

template<class Type>
matrix<Type> solve_upper_from_lower(const matrix<Type>& L,
                                    const matrix<Type>& B) {
  matrix<Type> X(B.rows(), B.cols());
  X.setZero();
  for (int rhs = 0; rhs < B.cols(); ++rhs) {
    for (int i = B.rows() - 1; i >= 0; --i) {
      Type z = B(i, rhs);
      for (int k = i + 1; k < B.rows(); ++k) {
        z -= L(k, i) * X(k, rhs);
      }
      X(i, rhs) = z / L(i, i);
    }
  }
  return X;
}

template<class Type>
matrix<Type> solve_spd_from_chol(const matrix<Type>& L,
                                 const matrix<Type>& B) {
  return solve_upper_from_lower(L, solve_lower(L, B));
}

// Eigen's sparse triangular solver skips a column when the current dense-RHS
// value compares equal to zero. That value-dependent optimization is unsafe
// while CppAD records at a parameter zero: a skipped update remains absent when
// the same tape is later evaluated at a nonzero parameter. Replay the LDLT
// solves explicitly, visiting every stored entry of the fixed L pattern.
template<class Type, class Solver>
matrix<Type> solve_sparse_ldlt_fixed(const Solver& ldlt,
                                     const matrix<Type>& B) {
  const int q = B.rows();
  matrix<Type> X = ldlt.permutationP() * B;
  Eigen::SparseMatrix<Type> L = ldlt.matrixL();

  // Unit-lower forward substitution, column-oriented.
  for (int rhs = 0; rhs < X.cols(); ++rhs) {
    for (int i = 0; i < q; ++i) {
      const Type xi = X(i, rhs);
      for (typename Eigen::SparseMatrix<Type>::InnerIterator it(L, i);
           it; ++it) {
        if (it.row() > i) {
          X(it.row(), rhs) -= xi * it.value();
        }
      }
    }
  }

  const vector<Type> D = ldlt.vectorD();
  for (int rhs = 0; rhs < X.cols(); ++rhs) {
    for (int i = 0; i < q; ++i) {
      X(i, rhs) /= D(i);
    }
  }

  // Unit-upper backward substitution using the same stored columns of L.
  for (int rhs = 0; rhs < X.cols(); ++rhs) {
    for (int i = q - 1; i >= 0; --i) {
      Type xi = X(i, rhs);
      for (typename Eigen::SparseMatrix<Type>::InnerIterator it(L, i);
           it; ++it) {
        if (it.row() > i) {
          xi -= it.value() * X(it.row(), rhs);
        }
      }
      X(i, rhs) = xi;
    }
  }
  return ldlt.permutationPinv() * X;
}

template<class Type>
Type objective_function<Type>::operator()() {
  DATA_SPARSE_MATRIX(Lambdat);
  DATA_IVECTOR(Lind);
  DATA_SPARSE_MATRIX(G);
  DATA_MATRIX(ZX);
  DATA_MATRIX(XX);
  DATA_MATRIX(ZY);
  DATA_MATRIX(XY);
  DATA_VECTOR(yy);
  DATA_INTEGER(n);
  DATA_IVECTOR(is_log);
  DATA_INTEGER(sparse_backend);

  PARAMETER_VECTOR(phi);

  const int q = Lambdat.rows();
  const int p = XX.rows();
  const int S = yy.size();
  const Type np = Type(n - p);

  vector<Type> theta(phi.size());
  for (int j = 0; j < phi.size(); ++j) {
    theta(j) = is_log(j) ? exp(phi(j)) : phi(j);
  }

  // Preserve every structural slot supplied by lme4, including entries whose
  // current parameter value is exactly zero. Lind follows compressed @x order.
  Eigen::SparseMatrix<Type> T = Lambdat;
  if (T.nonZeros() != Lind.size()) {
    return Type(1e20);
  }
  for (int k = 0; k < T.nonZeros(); ++k) {
    T.valuePtr()[k] = theta(Lind(k));
  }

  // A = I + T G T'. The sparse products are shared by both backends; no prune
  // is applied, so a parameter zero at tape construction does not remove its
  // symbolic contribution. The diagonal insertion guarantees the SPD pattern.
  Eigen::SparseMatrix<Type> A_sparse = T * G * T.transpose();
  for (int i = 0; i < q; ++i) {
    A_sparse.coeffRef(i, i) += Type(1.0);
  }
  A_sparse.makeCompressed();

  matrix<Type> BX = T * ZX;
  matrix<Type> BY = T * ZY;
  matrix<Type> rhs(q, p + S);
  rhs << BX, BY;

  matrix<Type> solved(q, p + S);
  Type ldL2 = Type(0.0);

  if (sparse_backend == 1) {
    // During CppAD recording Eigen chooses a permutation and symbolic LDLT
    // pattern from A_sparse. Subsequent fn/gr/he evaluations replay the same
    // sparse arithmetic with refreshed parameter values.
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Type> > ldlt;
    ldlt.compute(A_sparse);
    if (ldlt.info() != Eigen::Success) {
      return Type(1e20);
    }
    vector<Type> D = ldlt.vectorD();
    for (int i = 0; i < q; ++i) {
      if (!(D(i) > Type(0.0))) {
        return Type(1e20);
      }
      ldL2 += log(D(i));
    }
    solved = solve_sparse_ldlt_fixed(ldlt, rhs);
    if (ldlt.info() != Eigen::Success) {
      return Type(1e20);
    }
  } else {
    matrix<Type> A = A_sparse;
    matrix<Type> L = lchol(A);
    if (L.rows() == 0) {
      return Type(1e20);
    }
    for (int i = 0; i < q; ++i) {
      ldL2 += Type(2.0) * log(L(i, i));
    }
    solved = solve_spd_from_chol(L, rhs);
  }

  matrix<Type> AX = solved.leftCols(p);
  matrix<Type> AY = solved.rightCols(S);

  matrix<Type> M_raw = XX - BX.transpose() * AX;
  matrix<Type> M = Type(0.5) * (M_raw + M_raw.transpose());
  matrix<Type> Cu = XY - BX.transpose() * AY;

  matrix<Type> LX = lchol(M);
  if (LX.rows() == 0) {
    return Type(1e20);
  }
  Type ldRX2 = Type(0.0);
  for (int i = 0; i < p; ++i) {
    ldRX2 += Type(2.0) * log(LX(i, i));
  }

  matrix<Type> Cbeta = solve_lower(LX, Cu);
  vector<Type> PWRSS_s(S);
  Type mean_log_PWRSS = Type(0.0);
  bool valid = true;
  for (int s = 0; s < S; ++s) {
    PWRSS_s(s) = yy(s) - BY.col(s).dot(AY.col(s))
      - Cbeta.col(s).squaredNorm();
    if (!(PWRSS_s(s) > Type(0.0))) {
      valid = false;
    }
  }
  if (!valid) {
    return Type(1e20);
  }
  for (int s = 0; s < S; ++s) {
    mean_log_PWRSS += log(PWRSS_s(s) / np);
  }
  mean_log_PWRSS /= Type(S);

  const Type nll = Type(0.5) *
    (ldL2 + ldRX2 + np * mean_log_PWRSS);

  if (isDouble<Type>::value) {
    matrix<Type> beta = solve_upper_from_lower(LX, Cbeta);
    REPORT(beta);
  }
  REPORT(PWRSS_s);
  REPORT(LX);
  REPORT(ldL2);
  REPORT(ldRX2);
  return nll;
}
