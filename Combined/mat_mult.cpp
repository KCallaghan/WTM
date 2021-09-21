// C++ program to multiply two
// rectangular matrices
#include <bits/stdc++.h>
using namespace std;
 
// Multiplies two matrices mat1[][]
// and mat2[][] and prints result.
// (m1) x (m2) and (n1) x (n2) are
// dimensions of given matrices.
void multiply_comp(int m1, int m2, richdem::Array2D<double> &mat1, int n1, int n2,
              richdem::Array2D<std::complex<double>> &mat2, richdem::Array2D<std::complex<double>> &out_mat)
{
    int x, i, j;
    for (i = 0; i < m1; i++)
    {
        for (j = 0; j < n2; j++)
        {
            out_mat(i,j) = 0;
            for (x = 0; x < m2; x++)
            { 
                out_mat(i,j) += mat1(i,x) * mat2(x,j);
            }
        }
    }
}
 
// Driver code

int mat_mult_comp(int m1,int m2, int n1, int n2, richdem::Array2D<double> &mat1,richdem::Array2D<std::complex<double>> &mat2,richdem::Array2D<std::complex<double>> &out_mat)
{   
    // Function call
    multiply_comp(m1, m2, mat1, n1, n2, mat2, out_mat);
    return 0;
}


     //         mat_mult_comp(1,t_it-2,t_it-2,1,beta_matrix,sdelL_matrix,V_lm_result);
       //       mat_mult     (3,t_it-1,t_it-1,1,arp.my_sdel,arp.my_beta_konly, arp.rotational_V_lm);


// Multiplies two matrices mat1[][]
// and mat2[][] and prints result.
// (m1) x (m2) and (n1) x (n2) are
// dimensions of given matrices.




void multiply(int m1, int m2, richdem::Array2D<double> &mat1, int n1, int n2,
              richdem::Array2D<double> &mat2, richdem::Array2D<double> &out_mat)
{
    int x, i, j;
    for (i = 0; i < m1; i++)
    {
        for (j = 0; j < n2; j++)
        {
            out_mat(i,j) = 0;
            for (x = 0; x < m2; x++)
            { 
                out_mat(i,j) += mat1(i,x) * mat2(x,j);
            }
        }
    }
}
 
// Driver code
int mat_mult(int m1,int m2, int n1, int n2, richdem::Array2D<double> &mat1,richdem::Array2D<double> &mat2,richdem::Array2D<double> &out_mat)
{   
    // Function call
    multiply(m1, m2, mat1, n1, n2, mat2, out_mat);
    return 0;
}
