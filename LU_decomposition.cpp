#pragma once

#include"../Eigen/Dense"
#include<omp.h>

template<class T>
class LU_decomposition {
    using DenseMatrix = Eigen::Matrix<T,-1,-1>;
    using Vector = Eigen::Matrix<T, -1, 1>;

    DenseMatrix A;
    DenseMatrix L;
    DenseMatrix U;
    DenseMatrix P;
    int N;

    void compute() {
        U = A;
        L = DenseMatrix::Identity(N,N);

        for(int i=0;i<N;i++) {  // Diagonals
            T diag = U(i,i);
            for(int j=i+1;j<N;j++) {  // Everything down from current element
                T factor = U(j,i)/diag;
                for(int k=i;k<N;k++) {  // Diagonal and everything right of it 
                    U(j,k) -= factor * U(i,k);
                }
                L(j,i) = factor;
            }
        }      
    }

public:

    LU_decomposition(DenseMatrix A) : A(A) {
        N = A.rows();
        compute();
    }

    LU_decomposition() {

    }

    void set_matrix(const DenseMatrix& M) {
        A = M;
        N = A.rows();
        compute();
    }


    Vector solve(const Vector& b) {
        //Solving Ly = Pb
        Vector y = Vector::Zero(N);
        
        #pragma omp simd
        for (int i=0; i<N; i++) {
            T lx_sum = 0;
            for (int j=0; j<i; j++) {
                lx_sum += y(j) * L(i,j);
            }
            y(i) = (b(i)- lx_sum)/L(i,i);
        }

        // Solving Ub = y
        Vector x = Vector::Zero(N);

        #pragma omp simd
        for (int i=N-1; i>=0; i--) {
            T lx_sum = 0;
            for (int j=N-1; j>i; j--) {
                lx_sum += x(j) * U(i,j);
            }
            x(i) = (y(i)- lx_sum)/U(i,i);           
        }

        return x;
    }



};