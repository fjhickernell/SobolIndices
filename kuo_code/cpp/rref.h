#include <vector>

class BooleanMatrix{
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

    BooleanMatrix(std::vector< std::vector<bool> > M, int n, int m){
        this -> n = n;
        this -> m = m;
	this -> mat = M;
        gaussElimination();
    }
    template <size_t size_m>
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

    /* Does Gauss Elimination with partial pivoting on the matrix */
     void gaussElimination(){
        rank = n;
        for (int i = 0; i < n; i++){
            if (!mat[i][i]){
                int j;
                for (j = i+1; j < n && !mat[j][i]; j++);
                if (j == n){
                       rank--;
                       continue;
                }
                else
                    for (int k = i; k < m; k++){
                        bool t = mat[i][k];
                        mat[i][k] = mat[j][k];
                        mat[j][k] = t;
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

    /* Get the row rank of the boolean matrix
     * If you require the rank of the matrix, make sure that n > m.
     * i.e. if n < m, call the constructor over the transpose.
     */
    int getRank(){
        return rank;
    }
};
