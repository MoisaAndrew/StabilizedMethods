#include "methods_common.h"

extern "C" 
{
    /// <summary>
    /// Solver that includes the following methods:
    /// <para> = 0 TSRKC3 [https://doi.org/10.1016/j.cam.2024.116291]. </para>
    /// <para> = 1 TSRKC2 [https://doi.org/10.1016/j.cam.2024.115868], </para>
    /// <para> = 2 RKC [https://doi.org/10.1016/S0377-0427(97)00219-7], </para>
    /// </summary>
    /// <param name="n">Number of differential equations of the system</param>
    /// <param name="x">Initial point of integration</param>
    /// <param name="xend">End of the interval of integration, may be less than x</param>
    /// <param name="y">Initial value of the solution (array of length n)</param>
    /// <param name="f">Name (external) of function computing the value of f(x, y)</param>
    /// <param name="rho">Name (external) of a function giving the spectral radius of the Jacobian matrix.
    /// Supply a NULL pointer if iwork[0] = 1, and the solver will compute spectral radius internally</param>
    /// <param name="solout">Name (external) of a function providing the numerical solution during integration. 
    /// Supply a NULL pointer if iwork[2] = 0</param>
    /// <param name="atol">Absolute error tolerance. Can be scalar or vector of length n</param>
    /// <param name="rtol">Relative error tolerance</param>
    /// <param name="iwork">Integer array of length 10 that 
    /// gives information on how the problem is to be solved and 
    /// communicates statistics about the integration process.
    /// <para> iwork[0] </para>
    ///	<para> = 0 the solver attempts to compute the spectral radius internally; </para>
    ///	<para> = 1 rho returns an upper bound of the spectral radius of the Jacobian matrix </para>
    /// <para> iwork[1] </para>
    ///	<para> = 0 The Jacobian is not constant; </para>
    ///	<para> = 1 The Jacobian is constant </para>
    /// <para> iwork[2] </para>
    ///	<para> = 0 function solout is never called </para>
    ///	<para> = 1 function solout is called after every successful step </para>
    /// <para> iwork[3] </para>
    ///	<para> = 0 Atol is scalar </para>
    ///	<para> = 1 Atol is array of length n </para>
    /// </param>
    /// <returns>
    /// <para> 1 Successful computation x = xend; </para>
    /// <para> 3 Improper error control; </para>
    /// <para> 4 Stepsize becomes too small; </para>
    /// <para> 5 Invalid input parameters; </para>
    /// <para> 6 The method used in the solver to estimate the spectral radius did not converge. </para>
    /// <para> iwork[4] - Number of function evaluations </para>
    /// <para> iwork[5] - Number of steps </para>
    /// <para> iwork[6] - Number of accepted steps </para>
    /// <para> iwork[7] - Number of rejected steps </para>
    /// <para> iwork[8] - Number of evaluations of f used to estimate the spectral radius (equal to zero if iwork[0] = 1) </para>
    /// <para> iwork[9] - Maximum number of stages used </para>
    /// </returns>
    int rkc_solver(const unsigned n, const double x, const double xend, double* y,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double rtol,
        unsigned iwork[10], const int method);

    /// <summary>
    /// Monotonic stabilized method by B. Faleichik
    /// <para> Parameters are the same as for the rkc_solver </para>
    /// </summary>
    int mono(const unsigned n, const double x, const double xend, double* y,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double rtol,
        unsigned iwork[10]);


    int dumka3(const unsigned neqn, double* time, const double tend, const double h0,
        const double atol, const double rtol,
        const FcnEqDiff f, const Rho cour,
        double* y, double* z0, double* z1, double* z2,
        double* z3, double* oldEigenVector,
        unsigned iwork[12]);
}