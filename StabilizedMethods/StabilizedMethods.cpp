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
    VODPK_F,
};


static void run_method_test(
    const ProblemParams params, const enum MethodName method,
    const double fromtolp, const double totolp, const double tolpstep,
    const double* y0, double* y, const double* modelsolution,
    const Fcn fcn, const Rho rho, const Jac jac, const PSol psol,
    double* work, unsigned iwork[31], double* report[2],
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
        iwork[0] = 1;
        iwork[1] = params.isRhoDefined;
        iwork[2] = params.isSpradConst;
        iwork[3] = 0; iwork[4] = 0; iwork[5] = 0;
    }
    else if (method != VODPK_F)
    {
        iwork[0] = params.isRhoDefined;
        iwork[1] = params.isSpradConst;
        iwork[2] = 0; iwork[3] = 0; iwork[4] = 0; iwork[5] = 0;
    }
    else
    {
        iwork[0] = params.nDefault; iwork[1] = 0; 
        iwork[2] = 1; iwork[3] = 1;
        iwork[4] = 0; iwork[5] = 10000;
    }
    iwork[6] = 0; iwork[7] = 0; iwork[8] = 0;

    work[4] = 0; work[5] = 0; work[6] = 0; work[7] = 0;

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
            idid = rkc_solver(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork, work, 2);
            break;
        case DUMKA3:
            idid = dumka3(params.nDefault, &x, params.xend, h, atol, rtol, fcn, rho, y, 
                &(work[0]), &(work[params.nDefault]), &(work[2 * params.nDefault]), 
                &(work[3 * params.nDefault]), &(work[4 * params.nDefault]), iwork);
            break;
        case TSRKC2:
            idid = rkc_solver(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork, work, 1);
            break;
        case TSRKC3:
            idid = rkc_solver(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork, work, 0);
            break;
        case MONO:
            idid = mono(params.nDefault, x, params.xend, y, fcn, rho, solout_h, &atol, rtol, iwork);
            break;
        case VODPK_F:
            const unsigned itol = 1;
            const unsigned itask = 1;
            const unsigned iopt = 1;
            const unsigned liw = 31;
            const unsigned lrw = 61 + 18 * params.nDefault;
            const unsigned mf = 21;
            idid = 1;
            VODPK(fcn, &params.nDefault, y, &x, &params.xend, &itol, &rtol, &atol,
                &itask, &idid, &iopt, work, &lrw, iwork, &liw, jac, psol, &mf);
            break;
        }
        

        tf = high_resolution_clock::now();

        time_span = (params.nDefault < 50000 ? 1000 : 1) * duration_cast<duration<double>>(tf - ts).count();

        if (idid == 1 || (idid == 2 && method == VODPK_F))
        {
            if (printstats)
            {
                printf("atol = %1.2e, rtol = %1.2e\r\n", atol, rtol);
                if (method != VODPK_F)
                {
                    printf("fcn= %i step= %i accpt= %i rejct= %i fspr= %i stg= %i\r\n",
                        iwork[4], iwork[5], iwork[6], iwork[7], iwork[8], iwork[9]);
                }
                else
                {
                    printf("fcn= %i step= %i jac= %i psol= %i nonlin= %i conv= %i\r\n",
                        iwork[11], iwork[10], iwork[12], iwork[23], iwork[19], iwork[20]);
                }

                if (modelsolution != NULL)
                {
                    print_error(params.nDefault, modelsolution, y);
                }

                if (params.nDefault < 50000)
                {
                    printf(", time= %.3f ms\r\n\r\n", time_span);
                }
                else
                {
                    printf(", time= %.3f s\r\n\r\n", time_span);
                }
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
}


int main()
{
    ProblemParams* params;
    Fcn fcn;
    Rho rho;
    Jac jac;
    PSol psol;

    const double rtol = 1.0e-10, atol = 1.0e-10;

    // problem to test
    get_heat3d(&params, &fcn, &rho, &jac, &psol);
    unsigned iwork[31] = { params->isRhoDefined, params->isSpradConst };

    const bool printspatialerrors = false;
    // print stats like number of function evaluations, number of steps, ...
    const bool printstats = false;
    // print work/precision report
    const bool printreport = true;

    double* work = (double*)malloc((61 + 18 * params->nDefault) * sizeof(double));
    int idid;

    double x = params->x0;
    double h = params->h0;

    double* modelsolution = NULL;
    
    if (printspatialerrors)
    {
        modelsolution = params->y_exact(params->nDefault);
        if (modelsolution == NULL)
        {
            double* u = params->y0(params->nFine);
            if (8 * params->nFine > 61 + 18 * params->nDefault)
            {
                work = (double*)realloc(work, 8 * params->nFine * sizeof(double));
            }

            ROCK4F(&params->nFine, &x, &params->xend, &h, u, fcn, rho, &atol, &rtol, work, iwork, &idid);

            modelsolution = params->remap_solution(params->nFine, params->nDefault, u);

            free(u);

            x = params->x0;
            h = params->h0;
        }
    }

    const double* y0 = params->y0(params->nDefault);

    double* y = (double*)malloc(params->nDefault * sizeof(double));
    memcpy(y, y0, params->nDefault * sizeof(double));
    
    ROCK4F(&params->nDefault, &x, &params->xend, &h, y, fcn, rho, &atol, &rtol, work, iwork, &idid);
    
    if (idid == 1)
    {
        if (printspatialerrors)
        {
            print_error(params->nDefault, modelsolution, y);
            printf(", ");
        }
        printf("maxspr= %i minspr= %i\r\n\r\n", iwork[10], iwork[11]);

        modelsolution = (double*)malloc(params->nDefault * sizeof(double));
        memcpy(modelsolution, y, params->nDefault * sizeof(double));
    }
    
    // tol_max = 10 ^ fromtolp
    const double fromtolp = -3.0;
    // tol_min = 10 ^ totolp
    const double totolp = -8.0;
    // tol = 10 ^ (fromtolp - i * tolpstep)
    const double tolpstep = 0.5;

    const unsigned reportlength = (unsigned)(1.5 + (fromtolp - totolp) / tolpstep);
    double* report[2];
    if (printreport)
    {
        report[0] = (double*)malloc(reportlength * sizeof(double));
        report[1] = (double*)malloc(reportlength * sizeof(double));
    }

    
    /*
    printf("\r\n------------------------------------------rkcf------------------------------------------\r\n");

    run_method_test(*params, RKC_F, fromtolp, totolp, tolpstep, y0, y, modelsolution, 
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\r\n------------------------------------------rkcc------------------------------------------\r\n");

    run_method_test(*params, RKC_C, fromtolp, totolp, tolpstep, y0, y, modelsolution, 
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    
    
    printf("\r\n-----------------------------------------rock4f-----------------------------------------\r\n");

    run_method_test(*params, ROCK4_F, fromtolp, totolp, tolpstep, y0, y, modelsolution,
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    
    /*
    printf("\r\n-----------------------------------------rock2f-----------------------------------------\r\n");

    run_method_test(*params, ROCK2_F, fromtolp, totolp, tolpstep, y0, y, modelsolution, 
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\r\n-----------------------------------------dumka3-----------------------------------------\r\n");

    run_method_test(*params, DUMKA3, fromtolp, totolp, tolpstep, y0, y, modelsolution, 
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    */
    /*
    printf("\r\n-----------------------------------------tsrkc2-----------------------------------------\r\n");
    
    run_method_test(*params, TSRKC2, fromtolp, totolp, tolpstep, y0, y, modelsolution, 
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\r\n-----------------------------------------tsrkc3-----------------------------------------\r\n");
    
    run_method_test(*params, TSRKC3, fromtolp, totolp, tolpstep, y0, y, modelsolution,
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    
    /*
    printf("\r\n------------------------------------------mono------------------------------------------\r\n");

    run_method_test(*params, MONO, fromtolp, totolp, tolpstep, y0, y, modelsolution, 
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    */
    
    printf("\r\n-----------------------------------------vodpk------------------------------------------\r\n");

    run_method_test(*params, VODPK_F, fromtolp, totolp, tolpstep, y0, y, modelsolution,
        fcn, rho, jac, psol, work, iwork, report, printstats, printreport, reportlength);
    

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
