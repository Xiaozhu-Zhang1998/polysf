// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]
Rcpp::List lpfun(Eigen::VectorXd f_obj, 
                 Eigen::MatrixXd f_con,
                 Rcpp::StringVector f_dir,
                 Eigen::VectorXd f_rhs) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("lpSolve");
  Rcpp::Function f = pkg["lp"];
  return f(Rcpp::Named("objective.in", f_obj),
           Rcpp::Named("const.mat", f_con),
           Rcpp::Named("const.dir", f_dir),
           Rcpp::Named("const.rhs", f_rhs));
}


// [[Rcpp::export]]
Rcpp::IntegerVector im_equ(Eigen::MatrixXd A, Eigen::MatrixXd b) {
  int m = A.rows();
  int n = A.cols();
  // create B
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m, m+n);
  B.leftCols(n) = A;
  B.rightCols(m) = Eigen::MatrixXd::Identity(m, m);
  // create I
  Rcpp::IntegerVector I = {};
  // create f_obj
  Eigen::VectorXd f_obj = Eigen::VectorXd::Ones(m);
  // create f_con
  Eigen::MatrixXd f_con = Eigen::MatrixXd::Zero(m+n+1, m);
  f_con.topRows(m+n) = B.transpose();
  f_con.row(m+n) = b.transpose();
  // create f_dir
  Rcpp::StringVector f_dir(m+n+1);
  for (auto k=0; k< m+n+1 ; ++k) {
    if (k<n) {
      f_dir[k] = "=";
    } else if (k<m+n){
      f_dir[k] = ">=";
    } else{
      f_dir[k] = "=";
    }
  }
  // create f_rhs
  Eigen::VectorXd f_rhs = Eigen::VectorXd::Zero(m+n+1);
  // facial reduction
  for (auto j=n; j<n+m; ++j) {
    f_dir[j] = ">";
    f_rhs[j] = 0.5;
    Rcpp::List fit = lpfun(f_obj, f_con, f_dir, f_rhs);
    Rcpp::IntegerVector status = fit["status"];
    if (status[0] == 0) {
      Eigen::VectorXd soln = fit["solution"];
      Eigen::VectorXd z = B.transpose() * soln;
      for(auto k=n; k<m+n; ++k) {
        if( abs(z[k]) > 1e-7 ) {
          I.push_back(k);
        }
      }
    }
    f_dir[j] = ">=";
    f_rhs[j] = 0;
  }
  Rcpp::IntegerVector I_ = Rcpp::unique(I).sort();
  return(I_);
}


// [[Rcpp::export]]
Rcpp::IntegerVector lasso_fr(Eigen::MatrixXd Omega, Eigen::MatrixXd pi) {
  int m = Omega.rows();
  int n = Omega.cols();
  // create I
  Rcpp::IntegerVector I = {};
  // create f_obj
  Eigen::VectorXd f_obj = Eigen::VectorXd::Ones(m);
  // create f_con
  Eigen::MatrixXd f_con = Eigen::MatrixXd::Zero(n+1, m);
  f_con.topRows(n) = Omega.transpose();
  f_con.row(n) = pi.transpose();
  // create f_dir
  Rcpp::StringVector f_dir(n+1);
  for (auto k=0; k< n+1 ; ++k) {
    if (k<n) {
      f_dir[k] = ">=";
    } else {
      f_dir[k] = "==";
    } 
  }
  // create f_rhs
  Eigen::VectorXd f_rhs = Eigen::VectorXd::Zero(n+1);
  // facial reduction
  for (auto j=(n/2); j<n; ++j) {
    f_dir[j] = ">";
    f_rhs[j] = 0.5;
    Rcpp::List fit = lpfun(f_obj, f_con, f_dir, f_rhs);
    Rcpp::IntegerVector status = fit["status"];
    if (status[0] == 0) {
      Eigen::VectorXd soln = fit["solution"];
      Eigen::VectorXd z = Omega.transpose() * soln;
      for(auto k=(n/2); k<n; ++k) {
        if( abs(z[k]) > 1e-7 ) {
          I.push_back(k);
        }
      }
    }
    f_dir[j] = ">=";
    f_rhs[j] = 0;
    // if (j % 100 == 0) 
    //   Rcpp::Rcout << "finished " << j << " out of " << m << std::endl;
  } 
  Rcpp::IntegerVector I_ = Rcpp::unique(I).sort();
  return I_;
} 



Rcpp::IntegerVector MRank(Eigen::MatrixXd M) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("Matrix");
  Rcpp::Function f = pkg["rankMatrix"];
  Rcpp::IntegerVector val = f(M);
  return val;
}


//' @name ind_rows
//' @title Select independent rows of a matrix
//' @description Given any matrix \code{M}, this function selects the first \code{r=rank(M)} 
//' independent rows in \code{M} which form a row basis of \code{M}. 
//' @param M A matrix.
//' @return The function \code{ind_rows} returns a list containing the following components:
//'   \item{N}{A matrix with \code{r=rank(M)} rows. Its rows come from the selected independent rows in \code{M}.}
//'   \item{idx}{A logical vector. An element is \code{TRUE} if the corresponding row in \code{M} is selected, \code{FALSE} otherwise.}
//' @examples
//' # generate a matrix M of rank 2 with dependent rows
//' M = matrix(rnorm(16), nrow = 4)
//' M[2,] = rnorm(1) * M[1,] + rnorm(1) * M[3,]
//' M[4,] = rnorm(1) * M[2,] + rnorm(1) * M[3,]
//' 
//' # select the first 2 independent rows
//' ind_rows(M)
// [[Rcpp::export]]
Rcpp::List ind_rows(Eigen::MatrixXd M) {
  Rcpp::IntegerVector r = MRank(M);
  int m = M.rows();
  int n = M.cols();
  Rcpp::LogicalVector idx(m);
  idx.fill(0);
  int pt1 = 0; // pt
  int pt2 = 0;
  Eigen::MatrixXd M_cur = Eigen::MatrixXd::Zero(m,n);
  Rcpp::IntegerVector r_cur = MRank(M_cur);
  while(r_cur[0] < r[0]) {
    M_cur.row(pt1) = M.row(pt2);
    Rcpp::IntegerVector r_att = MRank(M_cur.topRows(pt1+1));
    if (r_att[0] > r_cur[0]) {
      r_cur = r_att;
      idx[pt2] = 1;
      pt1++;
    } 
    pt2++;
  }
  
  return Rcpp::List::create(Rcpp::Named("N") = M_cur.topRows(pt1),
                            Rcpp::Named("idx") = idx);
}



//' @name add_to_full
//' @title Extend a full-row-rank matrix to a full-rank square matrix.
//' @description Given any full-row-rank matrix \code{N} with \code{n=ncol(N)} and \code{r=rank(N)}. 
//' This function adds \code{n-r} independent rows to the matrix \code{N}, 
//' such that the resulting new matrix is a full-rank square matrix.  
//' @param N A full-row-rank matrix.
//' @return A full-rank square matrix whose first \code{r} rows come from \code{N}.
//' @examples 
//' # generate a full-row-rank matrix
//' N = matrix(rnorm(10), nrow = 2)
//' 
//' # add 3 independent rows to the matrix to form a full-rank square matrix
//' add_to_full(N)
// [[Rcpp::export]]
Eigen::MatrixXd add_to_full(Eigen::MatrixXd N) {
  Rcpp::IntegerVector r = MRank(N);
  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(N.cols(), N.cols());
  Eigen::MatrixXd N_cur = Eigen::MatrixXd::Zero(N.cols(), N.cols());
  N_cur.topRows(N.rows()) = N;
  Rcpp::IntegerVector r_cur = r;
  int count = 0;
  int pt = r[0];
  int idx = 0;
  while(count < N.cols() - r[0]) {
    N_cur.row(pt) = P.row(idx);
    Rcpp::IntegerVector r_att = MRank(N_cur.topRows(pt+1));
    if (r_att[0] > r_cur[0]) {
      pt++;
      r_cur = r_att;
      count++;
    } 
    idx++;
  }
  return N_cur;
}
