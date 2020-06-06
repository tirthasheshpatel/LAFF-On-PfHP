#define alpha( i,j ) A[ (j)*lda + (i) ]
#define  beta( i,j ) B[ (j)*ldb + (i) ]
#define gamma( i,j ) C[ (j)*ldc + (i) ]

void MyGemm(int m, int n, int k,
            const double *A, int lda,
            const double *B, int ldb,
                  double *C, int ldc)
{
    int j, p, i;

    for ( j=0 ; j<n ; j++ )
        for ( p=0 ; p<k ; p++ )
            for ( i=0 ; i<m ; i++ )
                gamma( i,j ) += alpha( i,p ) * beta( p,j );
}
