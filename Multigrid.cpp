#pragma once

#include"../Eigen/Sparse"
#include"../Eigen/Dense"
#include"LU_decomposition.cpp"
#include"Condition_Estimate.cpp"
#include<iostream>
#include<omp.h>
#include<time.h>
#include<map>
#include<tuple>
#include<chrono>
#include<map>
#include <math.h>

using namespace std;
using Eigen::VectorXi;


template<class T>
class Multigrid {
    using SparseMatrix = Eigen::SparseMatrix<T>;
    using DenseMatrix = Eigen::Matrix<T,-1,-1>;
    using Vector = Eigen::Matrix<T, -1, 1>;

    int cutoff;
    int size;
    T epsilon;
    int iterations;
    T post_norm;
    LU_decomposition<T> lu;

    std::chrono::_V2::system_clock::time_point startTime;
    std::chrono::_V2::system_clock::time_point endTime;

    SparseMatrix (*assemble_Matrix)(int);
    Vector b;
    Vector x;

    T infinity_norm(const Vector& vec) {
        return vec.cwiseAbs().maxCoeff();
    }

    Vector Jacobi(const SparseMatrix& A, const Vector& b, const Vector& x, double w=2/3.) {
        int N = A.rows();

        Vector x_new = A*x - b;
        auto diag = A.diagonal();

        for(int i=0;i<N;i++) {
            x_new(i) = x(i) - w*x_new(i)/diag(i);
        }
        
        return x_new;
    }


    Vector Landweber(const SparseMatrix& A, const Vector& b, const Vector& x, double w=2/3.) {
        double omega = 0.01;
        
        return x - omega*A.transpose() * (A*x - b);
    }

    Vector restriction(const Vector& error, int target_dim) {
        int size = error.size();
        Vector x = Vector::Zero(size/2 + 1);
        x[0] = error[0]*2/3. + error[1]*1/3.;
        int counter = 1;
        for(int i = 2; i<size-1;i+=2) {
            x[counter] = error[i-1]*0.25 + error[i]*0.5 + error[i+1]*0.25;
            counter++;
        }
        x[counter] = error[size-1]*2/3. + error[size-2]*1/3.;
        return x;
    }

    Vector prolongation(const Vector& error, int target_dim) {
        int size = error.size();
        Vector x = Vector::Zero(2*size - 1);

        x[0] = error[0]*2/3.;
        x[1] = error[0]*1/3. + error[1]*1/2.;

        int counter = 1;

        for(int i = 2; i<2*size-2;i+=2) {
            x[i] = error[counter];
            x[i+1] = (error[counter] + error[counter+1])/2.;
            counter++;
        }

        x[x.size()-2] = error[counter]*1/3. + error[counter-1]*1/2.;
        x[x.size()-1] = error[counter]*2/3.;
        return x;
    }


    Vector W_Cycle(const Vector& b, Vector x, int N, double damping) {

        SparseMatrix A = this->assemble_Matrix(N);
        x = Jacobi(A,b,x);
        Vector residual =  b - A*x;

        Vector restricted_residual = restriction(residual, N/2);

        Vector epsilon = Vector::Zero(restricted_residual.size());

        if(N/2 <= cutoff+1) {
            epsilon = lu.solve(restricted_residual);
            // epsilon = Jacobi(assemble_Matrix(restricted_residual.size()),restricted_residual,epsilon);
        } else {
            epsilon = W_Cycle(restricted_residual, epsilon, restricted_residual.size(), damping);
        }

        x = x + damping * prolongation(epsilon, restricted_residual.size()-1);

        x = Jacobi(A,b,x);

        //Second Phase

        residual =  b - A*x;

        restricted_residual = restriction(residual, N/2);

        if(N/2 <= cutoff+1) {
            epsilon = lu.solve(restricted_residual);
            // epsilon = Jacobi(assemble_Matrix(restricted_residual.size()),restricted_residual,epsilon);
        } else {
            epsilon = W_Cycle(restricted_residual, epsilon, restricted_residual.size(), damping);
        }

        x = x + damping * prolongation(epsilon, restricted_residual.size()-1);

        x = Jacobi(A,b,x);

        
        return x;
    }
 


public:

    Multigrid(SparseMatrix (*A)(int), Vector b, Vector x, int cutoff=256, T eps=pow(10,-3)) : assemble_Matrix(A), b(b), x(x), epsilon(eps), cutoff(cutoff) {
        size = b.size();
        iterations = -1;
        post_norm = -1;
    }

    double damping_inital_guess() {
        int N=17;
        int iterations = 2;

        map<T, double> damp_map;

        auto A = this->assemble_Matrix(N);

        lu.set_matrix(DenseMatrix(this->assemble_Matrix(9)));

        Vector b = Vector::Random(N)*10;

        for(double d = 0.02;d<=2;d+=0.02){

            Vector x = Vector::Zero(N);

            for(int i=0;i<iterations;i++) {
                x = W_Cycle(b,x,N,d);
            }
            // cout << "d: " << d << ", error: " << infinity_norm(A*x-b) << endl;
            damp_map.insert({infinity_norm(A*x-b), d});
        }

        double avg = 0;

        for(auto it=damp_map.begin(); it != std::next(damp_map.begin(), 5); ++it) {
            avg += it->second;
        }

        return avg/5.;
        
    }

    Vector& solve(unsigned int it_max=pow(10,2)) {
        startTime = std::chrono::high_resolution_clock::now();
        Condition_Estimate<T> ce(assemble_Matrix,10000);
        bool has_large_condition = ce.has_large_condition();
        cout << has_large_condition << endl;
        if (!has_large_condition) {
            epsilon = pow(10,-12);
            x = solve_jacobi(false);
        } else {
            solve_damp_control(false);
        }
        endTime = std::chrono::high_resolution_clock::now();
        
        return x;
    }



    Vector& solve_damp_control(bool time = true, unsigned int it_max=pow(2,31)) {
        if(time) startTime = std::chrono::high_resolution_clock::now();

        double damping = damping_inital_guess();
        cout << "Inital guess: " << damping << endl;
        double old_damp = damping/2.;
        
        auto A = this->assemble_Matrix(size);

        lu.set_matrix(DenseMatrix(this->assemble_Matrix(cutoff+1)));

        T error = infinity_norm(A*x - b);

        T rel_error = 1;
        int iteration_counter = 0;

        while (abs(rel_error) > epsilon && iteration_counter<it_max) {
            auto x_new = W_Cycle(b,x,size,damping);

            T tmp_error = infinity_norm(A*x_new - b);
            T tmp_rel = (error - tmp_error) / (tmp_error + error);

            if (tmp_rel>0) {
                x = x_new;
            } else {
                damping *= 0.5;
                cout << "Here " << damping << endl;
                x = W_Cycle(b,x,size,damping);

                tmp_error = infinity_norm(A*x - b);
                tmp_rel = (error - tmp_error) / (tmp_error + error);
            }
            


            T grad = (rel_error - tmp_rel) / (old_damp - damping);

            old_damp = damping;

            if(grad>=0) {
                double factor = (0.5*tanh(-2*tmp_rel - 0.5) + 1.5);
                damping*=factor;
            } else {
                double factor = (0.25*tanh(2*tmp_rel + 0.5) + 0.75);
                damping*=factor;
            }

            if (damping>1.9) {
                damping = 1.9;
            } else if(damping<0.1) {
                damping = 0.1;
            }

            cout <<"Iteration: " << iteration_counter << "\tDamping: " <<damping << "\t Relative error: " << tmp_rel << "\t absolute error: " << tmp_error << endl;
            
            rel_error = tmp_rel;
            error = tmp_error;

            iteration_counter++;
        }
    
        iterations = iteration_counter;
        post_norm = infinity_norm(A*x - b);

        if(time) endTime = std::chrono::high_resolution_clock::now();

        return x;
    }

    Vector& solve_jacobi(bool time = true, unsigned int it_max=pow(2,31)) {
        if(time) startTime = std::chrono::high_resolution_clock::now();
        auto A = this->assemble_Matrix(size);

        int iteration_counter = 0;
        T rel_error = 1;

        while(infinity_norm(A*x-b) > epsilon && iteration_counter < it_max ) {
            x = Jacobi(A,b,x);
            iteration_counter++;
            if (!(iteration_counter%5000)) {
                cout << infinity_norm(A*x - b) << endl;
            }
        }

        iterations = iteration_counter;
        post_norm = infinity_norm(A*x - b);

        if(time) endTime = std::chrono::high_resolution_clock::now();

        return x;
    }

    T guess_condition() {
        auto A = this->assemble_Matrix(size);
        T max_eigenvalue = 0;
        T min_eigenvalue = pow(10,31);

        for (int k=0; k<A.outerSize(); ++k) {
            T diagonal = 0;
            T non_diagonal_sum = 0;
            for (typename SparseMatrix::InnerIterator it(A,k); it; ++it)
            {
                if(it.row() == it.col()) {
                    diagonal = abs(it.value());
                } else {
                    non_diagonal_sum += abs(it.value());
                }
            }

            if ((diagonal - non_diagonal_sum) < min_eigenvalue) {
                min_eigenvalue = (diagonal - non_diagonal_sum);
            }

            if ((diagonal + non_diagonal_sum) > max_eigenvalue) {
                max_eigenvalue = (diagonal + non_diagonal_sum);
            }
        }

        return max_eigenvalue/min_eigenvalue;
    }

    Vector& get_result() {
        return x;
    }

    int get_iterations() {
        return iterations;
    }

    T get_error() {
        return post_norm;
    }

    void print_elapsed_time() {
        auto differenz = endTime-startTime;

        auto minutes = std::chrono::duration_cast< std::chrono::minutes >( differenz );
        differenz -= minutes;

        auto seconds = std::chrono::duration_cast< std::chrono::seconds >( differenz );
        differenz -= seconds;

        auto milliseconds = std::chrono::duration_cast< std::chrono::milliseconds >( differenz );
        differenz -= milliseconds;

        auto microseconds = std::chrono::duration_cast< std::chrono::microseconds >( differenz );

        cout << "Finished in: ";
        std::cout << minutes.count() << " Minutes, "
                    << seconds.count() << " Seconds, "
                    << milliseconds.count() << " Milliseconds, "
                    << microseconds.count() << " Microseconds." << std::endl;
    }








};