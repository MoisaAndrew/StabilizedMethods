#ifdef __cplusplus
extern "C" {
#endif

#define INTEL_FORTRAN
    /*#define IRIX*/                   /* MipsPro Compiler */
    /*#define SOLARIS*/                /* Sun Workshop Compiler */
    /*#define GCC*/

    /* without underscore, upper case*/
#ifdef INTEL_FORTRAN
#define UPPERCASE
#define FORTRAN_NAME(x) x
#endif

/* with underscore, lower case*/
#ifdef IRIX  
#define FORTRAN_NAME(x) x ## _
#endif

/* with underscore, lower case*/
#ifdef SOLARIS  
#define FORTRAN_NAME(x) x ## _
#endif

/* with underscore, lower case*/
#ifdef GCC
#define FORTRAN_NAME(x) x ## _
#endif


#ifdef UPPERCASE
#define ROCK2F FORTRAN_NAME(ROCK2)
#define ROCK4F FORTRAN_NAME(ROCK4) 

#define RKCF FORTRAN_NAME(RKC)
#define rkcdid_ FORTRAN_NAME(RKCDID)
#else
#define ROCK2F FORTRAN_NAME(rock2)
#define ROCK4F FORTRAN_NAME(rock4) 

#define RKCF FORTRAN_NAME(rkc) 
#define rkcdid_ FORTRAN_NAME(rkcdid)
#endif



/* This is the direct interface to FORTRAN                     */

    void ROCK2(const unsigned* N, double* X, const double* XEND, double* H, double* Y,
        void FCN(const unsigned*, const double*, const double*, double*),
        void RHO(const unsigned*, const double*, const double*, double*),
        const double* ATOL, const double* RTOL,
        double* WORK, unsigned* IWORK, int* IDID);


    void ROCK4(const unsigned* N, double* X, const double* XEND, double* H, double* Y,
        void FCN(const unsigned*, const double*, const double*, double*),
        void RHO(const unsigned*, const double*, const double*, double*),
        const double* ATOL, const double* RTOL,
        double* WORK, unsigned* IWORK, int* IDID);

    void RKC(const unsigned* N,
        void FCN(const unsigned*, const double*, const double*, double*),
        double* Y, double* X, const double* XEND,
        const double* RTOL, const double* ATOL,
        unsigned* INFO, double* WORK, int* IDID);

    extern struct {
        int nfe, nsteps, naccpt, nrejct, nfesig, maxm;
    } rkcdid_;


#ifdef __cplusplus
}
#endif
