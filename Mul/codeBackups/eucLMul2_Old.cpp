#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat eucL2(arma::mat x, int m) {
  int nY = x.n_rows - m;
  int d = x.n_cols; //数据维数
  // 声明 Y 矩阵（作用域覆盖整个函数）
  arma::mat Y;
  
  // 处理 m=0 的情况
  if (m == 0) {
     Y= x;  
    // 直接返回原矩阵（或复制前nY行，但此时nY=x.n_rows）
  }
  else {
    Y.set_size(nY, 2*d);  // 重新分配内存
    for (int i = 0; i < nY; i++) {
      // 检查索引是否越界
      //提取x的i以及i+m元素，横过来作为Y的第i行
      Y.row(i)=arma::join_rows(x.row(i), x.row(i+m));
    }
  }
  
  // 预分配距离矩阵
  arma::mat Euc(nY, nY, fill::zeros);
  for(int j = 0; j < nY; j++) {
    for(int k = j; k < nY; k++) {  // 利用对称性
      double dist = sum(square(Y.row(j) - Y.row(k)));
      Euc(j, k) = dist;
      if(j != k) {
        Euc(k, j) = dist;  // 对称赋值
      }
    }
  }
  
  
  return Euc;
}    