#include <memory>
#include <cstdlib>
#include <chrono>

#include "fortran.h"
#include "c.h"

#include "problems.h"
#include "utils.h"


enum MethodName
{
    ROCK4_F,
    ROCK2_F,
    RKC_F,
    RKC_C,
    DUMKA3,
    TSRKC2,
    TSRKC3,
    MONO,
};


static void run_method_test(
    const ProblemParams params, const enum MethodName method,
    const double fromtolp, const double totolp, const double tolpstep,
    const double* y0, double* y,
    const FcnEqDiff fcn, const Rho rho, const double* modelsolution,
    double* work, unsigned iwork[12], double* report[2],
    const bool printstats, const bool printreport, const unsigned reportlength)
{
    using namespace std::chrono;

    unsigned i;
    double atol, rtol;
    high_resolution_clock::time_point ts, tf;
    double time_span;
    double x, h;

    if (method == RKC_F)
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

        x = params.x0;
        h = params.h0;
        memcpy(y, y0, params.nDefault * sizeof(double));

        ts = high_resolution_clock::now();

        switch (method)
        {
        case ROCK4_F:
            ROCK4F(&params.nDefault, &x, &params.xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);
            break;
        case ROCK2_F:
            ROCK2F(&params.nDefault, &x, &params.xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);
            break;
        case RKC_F:
            RKCF(&params.nDefault, fcn, y, &x, &params.xend, &rtol, &atol, iwork, work, &idid);
            break;
        case RKC_C:
            idid = rkc_solver(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork, 2);
            break;
        case DUMKA3:
            idid = dumka3(params.nDefault, &x, params.xend, h, atol, rtol, fcn, rho, y, 
                &(work[0]), &(work[params.nDefault]), &(work[2 * params.nDefault]), 
                &(work[3 * params.nDefault]), &(work[4 * params.nDefault]), iwork);
            break;
        case TSRKC2:
            idid = rkc_solver(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork, 1);
            break;
        case TSRKC3:
            idid = rkc_solver(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork, 0);
            break;
        case MONO:
            idid = mono(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork);
            break;
        }
        

        tf = high_resolution_clock::now();

        time_span = 1000 * duration_cast<duration<double>>(tf - ts).count();

        if (idid == 1)
        {
            if (printstats)
            {
                printf("atol = %1.2e, rtol = %1.2e\r\n", atol, rtol);
                printf("fcn= %i step= %i accpt= %i rejct= %i fspr= %i stg= %i\r\n",
                    iwork[4], iwork[5], iwork[6], iwork[7], iwork[8], iwork[9]);

                if (modelsolution != NULL)
                {
                    print_error(params.nDefault, modelsolution, y);
                }

                printf(", time= %.3f ms\r\n\r\n", time_span);
            }
            if (printreport && modelsolution != NULL)
            {
                report[0][i] = get_absolute_error(params.nDefault, modelsolution, y);
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

    if (method == RKC_F)
    {
        i = iwork[0];
        iwork[0] = iwork[1];
        iwork[1] = iwork[2];
        iwork[2] = 1 - i;
    }
}


int main()
{
    ProblemParams* params;
    FcnEqDiff fcn;
    Rho rho;

    double rtol, atol;

    // problem to test
    get_raddiff(&params, &fcn, &rho);
    unsigned iwork[12] = { params->isRhoDefined, params->isJacConst, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    // print stats like number of function evaluations, number of steps, ...
    const bool printstats = false;
    // print work/precision report
    const bool printreport = true;

    double* work = (double*)malloc(8 * params->nDefault * sizeof(double));
    int idid;
    double* modelsolution = NULL;

    double x = params->x0;
    double h = params->h0;
    const double* y0 = params->y0(params->nDefault);
    
    double* y = (double*)malloc(params->nDefault * sizeof(double));
    memcpy(y, y0, params->nDefault * sizeof(double));
    
    rtol = 2.0e-15, atol = 2.0e-15;
    ROCK4F(&params->nDefault, &x, &params->xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);
    
    if (idid == 1)
    {
        modelsolution = (double*)malloc(params->nDefault * sizeof(double));
        memcpy(modelsolution, y, params->nDefault * sizeof(double));

        printf("maxspr= %i minspr= %i\r\n\r\n", iwork[10], iwork[11]);
    }
    
    // tol_max = 10 ^ fromtolp
    const double fromtolp = -1.0;
    // tol_min = 10 ^ totolp
    const double totolp = -8.0;
    // tol = 10 ^ (fromtolp - i * tolpstep)
    const double tolpstep = 1.0;

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
    printf("\r\n-----------------------------------------rock4f-----------------------------------------\r\n");

    run_method_test(*params, ROCK4_F, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\r\n-----------------------------------------rock2f-----------------------------------------\r\n");

    run_method_test(*params, ROCK2_F, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\r\n------------------------------------------rkcf------------------------------------------\r\n");

    run_method_test(*params, RKC_F, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\r\n------------------------------------------rkcc------------------------------------------\r\n");

    run_method_test(*params, RKC_C, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    
    /*
    printf("\r\n-----------------------------------------dumka3-----------------------------------------\r\n");

    run_method_test(*params, DUMKA3, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\r\n-----------------------------------------tsrkc2-----------------------------------------\r\n");
    
    run_method_test(*params, TSRKC2, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    
    
    printf("\r\n-----------------------------------------tsrkc3-----------------------------------------\r\n");
    
    run_method_test(*params, TSRKC3, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    
    
    printf("\r\n------------------------------------------mono------------------------------------------\r\n");

    run_method_test(*params, MONO, fromtolp, totolp, tolpstep, y0, y, fcn, rho,
        modelsolution, work, iwork, report, printstats, printreport, reportlength);
    

    if (printreport)
    {
        free(report[0]);
        free(report[1]);
    }
    free(params);
    free(work);
    free(y);
    free((void*)y0);
    free(modelsolution);
}
