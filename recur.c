#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "timer.h"

void integral_recur (int nmin, int nmax, double vals[]);
void integral_gen (int nmin, int nmax, double vals[]);
double ff (double x, void *params);
int adjust_rep_count (int count, double time, double tmin, double tmax);

double ff (double x, void *params)
{
    double n = *(double *) params;
    double f = exp (-x) * pow (x, n);

    return f;
}

void integral_recur (int nmin, int nmax, double vals[])
{
    int i;

    vals[nmax] = 0.;
    for (i = nmax - 1; i >= nmin; i--)
    {
        vals[i] = vals[i + 1] / (i + 1) + 1. / ((i + 1) * M_E);
    }
}

void integral_gen (int nmin, int nmax, double vals[])
{
    int i;
    double result, error;
    double a = 0., b = 1.;
    double abserr = 1.e-9, relerr = 1.e-9;
    double n;
    size_t np = 1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (np);

    gsl_function F;

    F.function = &ff;

    for (i = nmax; i >= nmin; i--)
    {
        n = (double) i;
        F.params = &n;

        gsl_integration_qag (&F, a, b, abserr, relerr, np, GSL_INTEG_GAUSS15,
            w, &result, &error);

        vals[i] = result;
    }


    gsl_integration_workspace_free (w);

}

#define N 100

int main (void)
{
    int i;
    double integral1[N + 1], integral2[N + 1];
    int nmax = N;
    int nmin = 10;

    integral_recur (nmin, nmax, integral1);
    integral_gen (nmin, nmax, integral2);

    for (i = nmin; i <= nmax; i++)
    {
        printf ("%5d    % 14.6f   % 14.6f   %g\n", i, integral1[i],
            integral2[i], fabs (integral1[i] - integral2[i]));
    }

    double time, time1, time2, tmin = 1., tmax = 2.;
    int count = 1000;

    do
    {

        timer_start ();
        for (i = 0; i < count; i++)
        {
            integral_recur (nmin, nmax, integral1);
        }
        time = timer_stop ();

        time1 = time / count;
        printf ("Recur: %10.2f usec %10.6f sec %10d\n", time1 * 1.e6, time,
            count);
        /*
         * adjust count such that cpu time is between
         * tmin and tmax
         */
        count = adjust_rep_count (count, time, tmin, tmax);
    }
    while ((time > tmax) || (time < tmin));

    count = 1000;

    do
    {

        timer_start ();
        for (i = 0; i < count; i++)
        {
            integral_gen (nmin, nmax, integral2);
        }
        time = timer_stop ();

        time2 = time / count;
        printf ("Integ: %10.2f usec %10.6f sec %10d\n", time1 * 1.e6, time,
            count);
        /*
         * adjust count such that cpu time is between
         * tmin and tmax
         */
        count = adjust_rep_count (count, time, tmin, tmax);
    }
    while ((time > tmax) || (time < tmin));

    printf ("Speedup: %.2f\n", time2 / time1);

    return 0;
}
