#pragma once

#include"../Eigen/Sparse"
#include"../Eigen/Dense"
#include"LU_decomposition.cpp"
#include<iostream>

template<class T>
class Condition_Estimate {
    using SparseMatrix = Eigen::SparseMatrix<T>;
    using DenseMatrix = Eigen::Matrix<T,-1,-1>;
    using Vector = Eigen::Matrix<T, -1, 1>;

    SparseMatrix (*assembleMatrix)(int);
    int iterations;

    T infinity_norm(const Vector& vec) {
        return vec.cwiseAbs().maxCoeff();
    }

    T find_sub_condition(const SparseMatrix SUB_A) {
        LU_decomposition<T> lu = LU_decomposition<T>(SUB_A);

        T lambda_max = power_iteration(SUB_A);
        T lambda_min = inverse_iteration(lu, SUB_A);
        
        return lambda_max/lambda_min;
    }

    T inverse_iteration(LU_decomposition<T>& lu, const SparseMatrix& A) {
        int N = A.rows();
        Vector x = Vector::Random(N);
        T norm = 1;
        int counter = 0;
        Vector x_prev = x;

        while(norm>pow(10,-4) && counter<iterations) {
            auto tmp = lu.solve(x);
            x_prev = x;

            if (!(counter%10)) {
                x = tmp / tmp.norm();
                norm = infinity_norm(x-x_prev/x_prev.norm());
            } else {
                x = tmp;
            }

            counter++;
        }

        return (x.dot(A*x))/(x.dot(x));
    }

    T power_iteration(const SparseMatrix& A) {
        int N = A.rows();
        Vector x = Vector::Random(N);
        T norm = 1;
        int counter = 0;
        Vector x_prev = x;

        while(norm>pow(10,-4) && counter<iterations) {
            auto tmp = A*x;
            x_prev = x;

            if (!(counter%10)) {
                x = tmp / tmp.norm();
                norm = infinity_norm(x-x_prev/x_prev.norm());
            } else {
                x = tmp;
            }

            counter++;
        }

        // for(int i=0;i<iterations;i++) {
        //     auto tmp = A*x;
        //     x = tmp / tmp.norm();
        // }
        return (x.dot(A*x))/(x.dot(x));
    }

public:
    Condition_Estimate(SparseMatrix (*A)(int), int it) : assembleMatrix(A), iterations(it) {
    }

    // T find_condition() {
    //     lu = LU_decomposition<T>(DenseMatrix(A));
    //     T lambda_max = power_iteration();
    //     T lambda_min = inverse_iteration();
    //     cout << "Max: " << lambda_max << ", Min: " << lambda_min << endl; 
    //     return lambda_max/lambda_min;
    // }

    bool has_large_condition() {
        T cond[5];
        for(int p=4;p<9;p++) {
            cond[p-4] = log(find_sub_condition(assembleMatrix(pow(2,p)+1)));
        }

        T first_entry = cond[0];
        T last_entry = cond[4]-first_entry;
        for(int i=0;i<5;i++) {
            cond[i] = (cond[i] - first_entry)/last_entry;
        }

        T curve[3];
        for(int i=0;i<3;i++) {
            curve[i] = abs(cond[i] - 2*cond[i+1] + cond[i+2]);
            if (curve[i] > 0.1) {
                return false;
            }
        }
        return true;

    }
};