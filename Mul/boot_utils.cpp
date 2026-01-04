// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>

using namespace Rcpp;

// 1. 极速 Stationary Bootstrap 索引生成器
// [[Rcpp::export]]
IntegerVector gen_sb_indices(int n, int block_size) {
    // 预分配内存，略多于 n 以防万一
    std::vector<int> indices;
    indices.reserve(n + block_size);
    
    // 随机数生成器
    // 使用 std::random_device 获得种子，或者你可以传入种子以复现
    std::random_device rd;
    std::mt19937 gen(rd());
    std::geometric_distribution<> geom_dist(1.0 / block_size);
    std::uniform_int_distribution<> unif_dist(0, n - 1);
    
    while(indices.size() < n) {
        int len = geom_dist(gen) + 1; // 几何分布长度
        int start = unif_dist(gen);   // 随机起始点 (0-based)
        
        for(int i = 0; i < len; ++i) {
            // 循环索引 (circular block bootstrap)
            indices.push_back((start + i) % n);
        }
    }
    
    // 截断到 n，并转为 R 的 1-based 索引 (方便 R 使用，或者直接在 C++ 用 0-based)
    // 这里返回 1-based 以兼容你现有的 R 代码逻辑 (y[idx, ])
    IntegerVector res(n);
    for(int i = 0; i < n; ++i) {
        res[i] = indices[i] + 1;
    }
    
    return res;
}

// 2. 极速 Fisher 组合统计量 (矩阵版)
// 替代 R 的 apply(matrix, 1, calc...)
// [[Rcpp::export]]
NumericVector calc_fisher_mat(NumericMatrix p_mat, NumericVector weights) {
    int n_rows = p_mat.nrow();
    int n_cols = p_mat.ncol();
    NumericVector results(n_rows);
    
    for(int i = 0; i < n_rows; ++i) {
        double sum_log = 0.0;
        bool has_na = false;
        
        for(int j = 0; j < n_cols; ++j) {
            double p = p_mat(i, j);
            
            if (NumericVector::is_na(p)) {
                has_na = true; 
                break;
            }
            
            // 边界保护
            if (p < 1e-10) p = 1e-10;
            if (p > 1.0 - 1e-10) p = 1.0 - 1e-10;
            
            sum_log += weights[j] * std::log(p);
        }
        
        if (has_na) {
            results[i] = NA_REAL;
        } else {
            results[i] = -2.0 * sum_log;
        }
    }
    return results;
}