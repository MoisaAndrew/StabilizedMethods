#include "methods_common.h"

extern "C" 
{
    int rock2c(const unsigned n, const double x, const double xend, double* h, double* y,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double* rtol,
        unsigned iwork[12]);


    int rkcc(const unsigned n, const double x, const double xend, double* y,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double rtol,
        unsigned iwork[10]);


    /// <summary>
    /// Two-Step RKC method of order 2
    /// </summary>
    /// <param name="n">Number of differential equations of the system</param>
    /// <param name="x0">Initial point of integration</param>
    /// <param name="x1">Additional point of integration. Set x1 = x0 if solution is unknown elsewhere</param>
    /// <param name="xend">End of the interval of integration, may be less than x0</param>
    /// <param name="h">Initial step size guess (usually between 1e-4 and 1e-6)</param>
    /// <param name="y0">Initial value of the solution (array of length n)</param>
    /// <param name="y1">Additional value of the solution (array of length n)</param>
    /// <param name="f">Name (external) of function computing the value of f(x, y)</param>
    /// <param name="rho">Name (external) of a function giving the spectral radius of the Jacobian matrix.
    /// Supply a dummy function if iwork[0] == 1, and TSRKC2 will compute spectral radius internally</param>
    /// <param name="solout">Name (external) of a function providing the numerical solution during integration. 
    /// Supply a dummy function if iwork[2] = 0.</param>
    /// <param name="atol">Absolute error tolerance. Can be scalar or vector of length n</param>
    /// <param name="rtol">Relative error tolerance. Can be scalar or vector of length n</param>
    /// <param name="iwork">Integer array of length 12 that 
    /// gives information on how the problem is to be solved and 
    /// communicates statistics about the integration process.
    /// <para> iwork[0] </para>
    ///	<para> = 0 TSRKC2 attempts to compute the spectral radius internally (rho can be a dummy function); </para>
    ///	<para> = 1 rho returns an upper bound of the spectral radius of the Jacobian matrix </para>
    /// <para> iwork[1] </para>
    ///	<para> = 0 The Jacobian is not constant; </para>
    ///	<para> = 1 The Jacobian is constant </para>
    /// <para> iwork[2] </para>
    ///	<para> = 0 function solout is called after every successful step </para>
    ///	<para> = 1 function solout is never called </para>
    /// <para> iwork[3] </para>
    ///	<para> = 0 Atol and rtol are scalar </para>
    ///	<para> = 1 Atol and rtol are arrays of length n </para>
    /// </param>
    /// <returns>
    /// <para> 1 Successful computation x = xend; </para>
    /// <para> -1 Invalid input parameters; </para>
    /// <para> -2 Stepsize becomes too small; </para>
    /// <para> -3 The method used in TSRKC2 to estimate the spectral radius did not converge. </para>
    /// <para> iwork[4] - Number of function evaluations </para>
    /// <para> iwork[5] - Number of steps </para>
    /// <para> iwork[6] - Number of accepted steps </para>
    /// <para> iwork[7] - Number of rejected steps </para>
    /// <para> iwork[8] - Number of evaluations of f used to estimate the spectral radius (equal to zero if iwork[0] = 1) </para>
    /// <para> iwork[9] - Maximum number of stages used </para>
    /// <para> iwork[10] - Maximum value of the estimated bound for the spectral radius </para>
    /// <para> iwork[11] - Minimum value of the estimated bound for the spectral radius </para>
    /// </returns>
    int tsrkc2(const unsigned n,
        const double x0, const double x1, const double xend,
        double* h, double* y0, double* y1,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double* rtol,
        unsigned iwork[12]);
}