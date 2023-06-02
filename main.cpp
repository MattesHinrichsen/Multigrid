#include"Multigrid.cpp"
#include"LU_decomposition.cpp"
#include"Condition_Estimate.cpp"
#include<time.h>
#include"Timer.cpp"

using T = double;

using SparseMatrix = Eigen::SparseMatrix<T>;
using DenseMatrix = Eigen::Matrix<T,-1,-1>;
using Vector = Eigen::Matrix<T, -1, 1>;

Vector assembleRHS(int N) {
    return Vector::Random(N)*10;
}

SparseMatrix assembleMatrix(int N) {
    SparseMatrix m(N,N);
    m.reserve(VectorXi::Constant(N, 5));

    T diagonal_value = 2;
    T up_value = -1;
    T down_value = -1;
    T second = 0;

    for(int i=0; i<N;i++){
        m.coeffRef(i,i) = diagonal_value*(N-1);
        if(i-1>=0) m.coeffRef(i-1,i) = down_value*(N-1);
        if(i+1<N) m.coeffRef(i+1,i) = up_value*(N-1);

        if(i-2>=0 && second != 0) m.coeffRef(i-2,i) = second*(N-1);
        if(i+2<N && second != 0) m.coeffRef(i+2,i) = second*(N-1);
    }
    m.makeCompressed();
    return m;
}

int main() {
    srand (time(NULL));
    // srand(0);  
    int N = pow(2,20)+1;
    auto A = assembleMatrix(N);
    auto b = assembleRHS(N);
    Vector x = Vector::Zero(N);

    //Multigrid

    Multigrid<double> mg(assembleMatrix, b, x);

    mg.solve();

    cout << N << ": " << mg.get_error() << " after " << mg.get_iterations() << " iterations." << endl;
    mg.print_elapsed_time();

   /*  //------------------------------------------------------------------------------
    cout << "2nd Run:" << endl;
    Multigrid<double> mg2(assembleMatrix, b, x);

    mg2.solve_damp_control();

    cout << N << ": " << mg2.get_error() << " after " << mg2.get_iterations() << " iterations." << endl;
    mg2.print_elapsed_time();


    //------------------------------------------------------------------------------
    cout << "3rd Run:" << endl;
    Multigrid<double> mg3(assembleMatrix, b, x, 256, pow(10,-12));

    mg3.solve_jacobi();

    cout << N << ": " << mg3.get_error() << " after " << mg3.get_iterations() << " iterations." << endl;
    mg3.print_elapsed_time();    */

    /*
    //LU decomp
    DenseMatrix A(3,3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 10;
    
    Vector b(3);
    b << 1,2,3;

    LU_decomposition<double> lu(A);

    cout << lu.solve(b) << endl;
   
    */

}