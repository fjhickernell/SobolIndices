#ifndef BOOLEANMATRIX_HH
#define BOOLEANMATRIX_HH

#include <vector>

class BooleanMatrix{
    /*This class contains a method gaussElimination()
     *that computes the rank of the matrix in Z2
     */

    std::vector< std::vector<bool> > mat; //boolean matrix
    int n, m;           //size of matrix nxm
    int rank;           //rank of the matrix

    public:

    /*Constructor
     * Required Parameters:
     * M ==> boolean matrix
     * n ==> number of rows
     * m ==> number of columns
     */

    BooleanMatrix(){
    this -> n = 0;
    this -> m = 0;
    }

    BooleanMatrix(std::vector< std::vector<bool> > M, int n, int m){
        this -> n = n;
        this -> m = m;
	this -> mat = M;
        gaussElimination();
    }
/*    template <size_t size_m>
    BooleanMatrix(bool M[][size_m], int n, int m){
        this -> n = n;
        this -> m = m;
        for (int i = 0; i < n; i++){
            std::vector<bool> row(m);
            for (int j = 0; j < m; j++) row[j] = M[i][j];
            mat.push_back(row);         
        }
        gaussElimination();
    }
*/
    void set(std::vector< std::vector<bool> > M, int n, int m){
        this -> n = n;
        this -> m = m;
        this -> mat = M;
        gaussElimination();
    }

    /* Does Gauss Elimination with partial pivoting on the matrix */
     void gaussElimination(){
        rank = n;
        for (int i = 0; i < n; i++){ //std::cout << i << std::endl;
            if (!mat[i][i]){
		int j, p;
                for (j = i+1; j < n && !mat[j][i]; j++);
		for (p = i+1; p < m && !mat[i][p]; p++);
                if (j == n && p == m){
                       rank--;
                       continue;
                }
                else if (j < n){
                    for (int k = i; k < m; k++){
                        bool t = mat[i][k];
                        mat[i][k] = mat[j][k];
                        mat[j][k] = t;
                    }
		}
            }
            for (int j = i+1; j < n; j++){
                if (mat[j][i]){
                    for (int k = i; k < m; k++)
                        mat[j][k] = mat[j][k] - mat[i][k];
                }
            }
        }
    }

    void print(){
            for (unsigned pp = 0U; pp < n; pp++){
                for (unsigned pp2 = 0U; pp2 < m; pp2++){
                        std::cout << (unsigned) mat[pp][pp2] << " ";
                }
	    std::cout << std::endl;
            }
    }

    /* Get the row rank of the boolean matrix
     * If you require the rank of the matrix, make sure that n > m.
     * i.e. if n < m, call the constructor over the transpose.
     */
    int getRank(){
        return rank;
    }
};

#endif // BOLEANMATRIX_HH
