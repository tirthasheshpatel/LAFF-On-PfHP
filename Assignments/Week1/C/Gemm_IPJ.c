#define alpha( i,j ) A[ (j)*lda + (i) ]
#define  beta( i,j ) B[ (j)*ldb + (i) ]
#define gamma( i,j ) C[ (j)*ldc + (i) ]

void MyGemm(int m, int n, int k,
            double *A, int lda,
            double *B, int ldb,
            double *C, int ldc)
{
    int i, p, j;

    for ( i=0 ; i<m ; i++ )
        for ( p=0 ; p<k ; p++ )
            for ( j=0 ; j<n ; j++ )
                gamma( i,j ) += alpha( i,p ) * beta( p,j );
}
