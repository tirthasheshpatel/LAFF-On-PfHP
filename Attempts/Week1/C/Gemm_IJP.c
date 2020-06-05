#define alpha( i,j ) A[ (j)*lda + (i) ] // march through the matrix in F_CONTIGUOUS (column-major) manner
#define  beta( i,j ) B[ (j)*ldb + (i) ] // march through the matrix in F_CONTIGUOUS (column-major) manner
#define gamma( i,j ) C[ (j)*ldc + (i) ] // march through the matrix in F_CONTIGUOUS (column-major) manner

void MyGemm(int m, int n, int k, double *A, int lda,
            double *B, int ldb, double *C, int ldc)
{
    int i, j, p;
    for ( i=0 ; i<m ; i++ )
        for ( j=0 ; j<n ; j++ )
            for ( p=0 ; p<k ; p++ )
                gamma( i,j ) += alpha( i,p ) * beta ( p,j );
}
