//#include <tr1/cmath> // fonctionne avec intel oneapi
#include <cmath>


#include "caiFilters.h"

double wk(double* a,
          double* k)
{
    const double ak = *a * *k;
    return 0.5 * ak * ak * std::cyl_bessel_k(2.0,ak );
    //return 0.5 * ak * ak * std::tr1::cyl_bessel_k(2.0,ak ); // fonctionne avec intel oneapi
}

double wtildek(double* a,
               double* k)
{
    return sqrt( wk(a,k)  );
}
