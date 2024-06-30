#include <memory.h>
#include <stdlib.h>
#include <chrono>

#include "fortran.h"
#include "c.h"

#include "problems.h"

#include "utils.h"


enum MethodName
{
    ROCK4_F,
    ROCK2_F,
    ROCK2_C,
    RKC_F,
    RKC_C,
    DUMKA3,
    TSRKC2,
    TSRKC3
};


void run_method_test(
    const unsigned n, const enum MethodName method,
    const double fromtolp, const double totolp, const double tolpstep,
    const double x0, const double xend, double* yp, double* y0, double* y, const double h0,
    const FcnEqDiff fcn, Rho rho, const double* modelsolution,
    double* work, unsigned iwork[12], double* report[2],
    const bool printstats, const bool printreport, const unsigned reportlength)
{
    using namespace std::chrono;

    unsigned i;
    double atol, rtol;
    high_resolution_clock::time_point ts, tf;
    double time_span;
    double x, h;

    if (method == RKC_F || method == RKC_C)
    {
        i = iwork[2];
        iwork[1] = iwork[0];
        iwork[2] = iwork[1];
        iwork[0] = 1 - i;
    }

    i = 0;
    for (double tolp = fromtolp; tolp >= totolp; tolp -= tolpstep)
    {
        int idid = 0;

        atol = rtol = pow(10., tolp);

        x = x0;
        h = h0;
        memcpy(y, y0, n * sizeof(double));
        memcpy(yp, y0, n * sizeof(double));

        ts = high_resolution_clock::now();

        switch (method)
        {
        case ROCK4_F:
            ROCK4F(&n, &x, &xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);
            break;
        case ROCK2_F:
            ROCK2F(&n, &x, &xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);
            break;
        case ROCK2_C:
            idid = rock2c(n, x, xend, &h, y, fcn, rho, solout_h, &atol, &rtol, iwork);
            break;
        case RKC_F:
            RKCF(&n, fcn, y, &x, &xend, &rtol, &atol, iwork, work, &idid);
            break;
        case RKC_C:
            idid = rkcc(n, x, xend, y, fcn, rho, solout_h, &atol, rtol, iwork);
            break;
        case DUMKA3:
            idid = dumka3(n, &x, xend, h, atol, rtol, fcn, rho, 
                y, &(work[0]), &(work[n]), &(work[2 * n]), &(work[3 * n]), &(work[4 * n]), iwork);
            break;
        case TSRKC2:
            idid = tsrkc2(n, x, x, xend, &h, yp, y, fcn, rho, solout_h, &atol, &rtol, iwork);
            break;
        case TSRKC3:
            idid = tsrkc3(n, x, x, xend, &h, yp, y, fcn, rho, solout_h, &atol, &rtol, iwork);
            break;
        }
        

        tf = high_resolution_clock::now();

        time_span = 1000 * duration_cast<duration<double>>(tf - ts).count();

        if (idid == 1)
        {
            if (printstats)
            {
                printf("atol = %1.2e, rtol = %1.2e\n", atol, rtol);
                printf("fcn= %i step= %i accpt= %i rejct= %i fspr= %i stg= %i\n",
                    iwork[4], iwork[5], iwork[6], iwork[7], iwork[8], iwork[9]);

                if (modelsolution != NULL)
                {
                    print_error(n, modelsolution, y);
                }

                printf(", time= %.3f ms\n\n", time_span);
            }
            if (printreport && modelsolution != NULL)
            {
                report[0][i] = get_absolute_error(n, modelsolution, y);
                report[1][i] = time_span;
                i++;
            }
        }
        else
        {
            if (printreport && modelsolution != NULL)
            {
                report[0][i] = -1;
                report[1][i] = time_span;
                i++;
            }
        }
    }
    if (printreport && modelsolution != NULL)
    {
        print_report(reportlength, report);
    }

    if (method == RKC_F || method == RKC_C)
    {
        i = iwork[0];
        iwork[0] = iwork[1];
        iwork[1] = iwork[2];
        iwork[2] = 1 - i;
    }
}


int main()
{
    unsigned n;
    FcnEqDiff fcn;
    Rho rho = NULL;
    double x, x0, xend, h, h0;
    double* y, * y0, * yp;

    double rtol, atol;

    // problem to test
    get_raddiff(&n, &fcn, &rho, &x0, &h0, &xend, &y0);
    unsigned iwork[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // = 0 - methods attempt to compute the spectral radius internally;
    // = 1 - rho returns an upper bound of the spectral radius;
    iwork[0] = 1;
    // = 0 - The Jacobian is not constant; = 1 - The Jacobian is constant;
    iwork[1] = 1;
    // = 0 - function solout is called after every successful step;
    // = 1 - function solout is never called
    iwork[2] = 0;
    // = 0 - Atol and rtol are scalar; = 1 - Atol and rtol are arrays of length n
    iwork[3] = 0;

    double* work = (double*)malloc(8 * n * sizeof(double));
    int idid;
    double* modelsolution = NULL;

    x = x0;
    h = h0;
    y = (double*)malloc(n * sizeof(double));
    yp = (double*)malloc(n * sizeof(double));
    memcpy(y, y0, n * sizeof(double));
    memcpy(yp, y0, n * sizeof(double));
    
    rtol = 2.0e-15, atol = 2.0e-15;
    ROCK4F(&n, &x, &xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);

    if (idid == 1)
    {
        modelsolution = (double*)malloc(n * sizeof(double));
        memcpy(modelsolution, y, n * sizeof(double));

        printf("maxspr= %i minspr= %i\n\n", iwork[10], iwork[11]);
    }
    
    // tol_max = 10 ^ fromtolp
    const double fromtolp = -1.0;
    // tol_min = 10 ^ totolp
    const double totolp = -6.0;
    // tol = 10 ^ (fromtolp - i * tolpstep)
    const double tolpstep = 1.0;
    // print stats like number of function evaluations, number of steps, ...
    const bool printstats = false;
    // print work/precision report
    const bool printreport = true;
    const unsigned reportlength = (unsigned)(1.5 + (fromtolp - totolp) / tolpstep);
    double* report[2];
    if (printreport)
    {
        report[0] = (double*)malloc(reportlength * sizeof(double));
        report[1] = (double*)malloc(reportlength * sizeof(double));
    }

    
    // Testing of different methods in a row at one start can lead to a large spread of results.
    // For more accurate time measurements, it is recommended to comment all the methods except one.
    // Do not forget to test speed of the methods in Release configuration only.
    /*
    printf("\n-----------------------------------------rock4f-----------------------------------------\n");

    run_method_test(n, ROCK4_F, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\n-----------------------------------------rock2f-----------------------------------------\n");

    run_method_test(n, ROCK2_F, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\n-----------------------------------------rock2c-----------------------------------------\n");
    
    run_method_test(n, ROCK2_C, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\n------------------------------------------rkcf------------------------------------------\n");

    run_method_test(n, RKC_F, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\n------------------------------------------rkcc------------------------------------------\n");

    run_method_test(n, RKC_C, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    
    /*
    printf("\n-----------------------------------------dumka3-----------------------------------------\n");

    run_method_test(n, DUMKA3, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\n-----------------------------------------tsrkc2-----------------------------------------\n");
    
    run_method_test(n, TSRKC2, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    
    
    printf("\n-----------------------------------------tsrkc3-----------------------------------------\n");
    
    run_method_test(n, TSRKC3, fromtolp, totolp, tolpstep, x0, xend, yp, y0, y, h0, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    

    if (printreport)
    {
        free(report[0]);
        free(report[1]);
    }
    free(work);
    free(y);
    free(y0);
    free(yp);
    free(modelsolution);
}
