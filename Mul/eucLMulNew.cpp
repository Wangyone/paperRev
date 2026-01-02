#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/**
 * 计算时间序列重构后的欧式距离平方矩阵
 *
 * 对多元时间序列进行重构后，计算重构向量间的欧式距离平方矩阵。
 * 使用矩阵运算优化避免循环，显著提高计算效率。
 *
 * @param x 输入的时间序列矩阵，维度为 N × d
 * @param lag 滞后参数，决定重构时的滞后步长
 * @return 欧式距离平方矩阵，维度为 (N-lag) × (N-lag)
 */
// [[Rcpp::export]]
arma::mat eucL(arma::mat x, int lag)
{
  // 获取数据维度
  int N = x.n_rows; // 时间点总数
  int d = x.n_cols; // 原始维度

  // 计算重构向量数量
  int nY = N - lag; // 重构向量数量

  // 边界情况处理
  if (nY <= 0)
  {
    return arma::mat(); // 返回空矩阵
  }

  // 构建重构矩阵 Y
  arma::mat Y;

  if (lag == 0)
  {
    // lag=0: 直接使用原始数据
    Y = x;
  }
  else
  {
    // lag>0: 拼接当前点和滞后点
    Y.set_size(nY, 2 * d);

    // 前半部分: 复制x的前nY行
    Y.cols(0, d - 1) = x.rows(0, nY - 1);

    // 后半部分: 复制x的第lag到nY+lag-1行
    Y.cols(d, 2 * d - 1) = x.rows(lag, nY + lag - 1);
  }

  // 计算欧式距离平方矩阵
  // 公式: ||Y_i - Y_j||^2 = ||Y_i||^2 + ||Y_j||^2 - 2⟨Y_i, Y_j⟩

  // 1. 计算每行的平方和
  arma::vec row_sq_sum = sum(square(Y), 1);

  // 2. 计算点积矩阵
  arma::mat dot_prod = Y * Y.t();

  // 3. 使用广播机制计算距离矩阵
  arma::mat broad_row = row_sq_sum * arma::ones<arma::rowvec>(nY);
  arma::mat broad_col = arma::ones<arma::colvec>(nY) * row_sq_sum.t();

  arma::mat Euc = broad_row + broad_col - 2.0 * dot_prod;

  // 确保对角线为0
  Euc.diag().zeros();

  return Euc;
}