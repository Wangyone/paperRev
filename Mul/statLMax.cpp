#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/**
 * @brief 计算变点监测统计量 (直接返回最大值)
 * @author Wang
 */
// [[Rcpp::export]]
double statL_max(arma::mat x, int T, int m, double alphaV, 
                 arma::vec stat_weight,
                 String kern, double kern_para)
{
    // 数据维度
    int N = x.n_rows;
    int d = x.n_cols;
    int nY = N - m;
    
    // 边界检查，如果没有足够的数据点，返回 0 或 -Inf
    if (nY <= 0) {
        return 0.0; 
    }
    
    // ============================================================
    // 1. 构建重构矩阵 Y
    // ============================================================
    arma::mat Y;
    
    if (m == 0) {
        Y = x;
    } else {
        Y.set_size(nY, 2 * d);
        Y.cols(0, d - 1) = x.rows(0, nY - 1);
        Y.cols(d, 2 * d - 1) = x.rows(m, nY + m - 1);
    }
    
    // ============================================================
    // 2. 距离矩阵 (逻辑不变)
    // ============================================================
    arma::mat Wjk(nY, nY);
    
    if (kern == "energy") {
        arma::vec row_sq_sum = sum(square(Y), 1);
        arma::mat dot_prod = Y * Y.t();
        arma::mat Euc_sq = repmat(row_sq_sum, 1, nY) +
            repmat(row_sq_sum.t(), nY, 1) -
            2.0 * dot_prod;
        Euc_sq = arma::max(Euc_sq, arma::zeros<arma::mat>(nY, nY));
        Wjk = pow(Euc_sq, alphaV / 2.0);
        
    } else if (kern == "gauss") {
        arma::vec row_sq_sum = sum(square(Y), 1);
        arma::mat dot_prod = Y * Y.t();
        arma::mat Euc_sq = repmat(row_sq_sum, 1, nY) +
            repmat(row_sq_sum.t(), nY, 1) -
            2.0 * dot_prod;
        Euc_sq = arma::max(Euc_sq, arma::zeros<arma::mat>(nY, nY));
        Wjk = exp(-Euc_sq / (4.0 * kern_para));
        
    } else if (kern == "quadexp") {
        double d1 = 2.0 * kern_para;
        double d2 = 2.0 * d1; 
        
        for (int j = 0; j < nY; j++) {
            for (int k = 0; k < nY; k++) {
                arma::rowvec diff = Y.row(j) - Y.row(k);
                arma::vec v2 = -square(diff.t());
                arma::vec v3 = exp(v2 / d2) % (d1 + v2) / d1;
                double prod_val = 1.0;
                for (unsigned int i = 0; i < v3.n_elem; i++) {
                    prod_val *= v3(i);
                }
                Wjk(j, k) = prod_val;
            }
        }
    } else {
        stop("Unknown kernel type.");
    }
    
    // ============================================================
    // 3. 计算统计量序列的最大值
    // ============================================================
    int Tm = T - m;
    int stat_len = nY - Tm;
    
    if (stat_len <= 0) return 0.0;
    
    // 使用累积和优化
    arma::mat Wjk_cum = arma::cumsum(arma::cumsum(Wjk, 0), 1);
    arma::mat Wjk_cum_pad(nY + 1, nY + 1, arma::fill::zeros);
    Wjk_cum_pad.submat(1, 1, nY, nY) = Wjk_cum;
    
    auto submatrix_sum = [&Wjk_cum_pad](int i, int j) -> double {
        return Wjk_cum_pad(i + 1, j + 1);
    };
    
    double Tm_sq_inv = 1.0 / (Tm * Tm);
    
    // 【修改点】：初始化最大值为极小值
    double max_stat_val = -arma::datum::inf; 
    
    for (int tt = 0; tt < stat_len; tt++) {
        int Tmt = Tm + tt + 1;
        
        double sum_A = submatrix_sum(Tmt - 1, Tm - 1);
        double sum_B = submatrix_sum(Tm - 1, Tm - 1);
        double sum_C = submatrix_sum(Tmt - 1, Tmt - 1);
        
        double stat_val = 0.0;
        double Tmt_inv = 1.0 / Tmt;
        double Tmt_sq_inv = 1.0 / (Tmt * Tmt);
        
        if (kern == "energy") {
            stat_val = 2.0 * sum_A * Tmt_inv / Tm -
                sum_B * Tm_sq_inv -
                sum_C * Tmt_sq_inv;
        } else if (kern == "gauss") {
            stat_val = -2.0 * sum_A * Tmt_inv / Tm +
                sum_B * Tm_sq_inv +
                sum_C * Tmt_sq_inv;
            stat_val *= pow(datum::pi / kern_para, d);
        } else if (kern == "quadexp") {
            stat_val = 2.0 * sum_A * Tmt_inv / Tm -
                sum_B * Tm_sq_inv -
                sum_C * Tmt_sq_inv;
            stat_val *= -1.0; 
        }
        
        // 【修改点】：直接计算加权后的值并更新最大值
        double current_weighted_val = stat_weight(tt) * stat_val;
        
        if (current_weighted_val > max_stat_val) {
            max_stat_val = current_weighted_val;
        }
    }
    
    return max_stat_val;
}