#define A( i,j ) *( ap + (j)*lda + i ) // march through the matrix in F_CONTIGUOUS (column-major) manner
#define B( i,j ) *( bp + (j)*ldb + i ) // march through the matrix in F_CONTIGUOUS (column-major) manner

#define dabs( x ) ( (x) < 0 ? -(x) : x )

double MaxAbsDiff(int m, int n, double *ap, int lda, double *bp, int ldb)
/**
 * Get the maximum absolute difference between two matrices
 * 
 * @param m number of rows in the array
 * @param n number of columns in the array
 * @param ap pointer to the first array
 * @param lda size of columns in bytes
 * @param bp pointer to the second array
 * @param ldb size of an entry in second array
 * 
 * @returns void
 */
{
    double diff = 0.;
    int i, j;
    for ( i=0 ; i<m ; i++ )
        for ( j=0 ; j<n ; j++ )
            if ( dabs( A( i,j ) - B( i,j ) ) > diff )
                diff = dabs( A( i,j ) - B( i,j ) );
    return diff;
}
