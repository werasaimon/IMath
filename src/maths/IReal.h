#ifndef IREAL_H
#define IREAL_H

#include <cmath>

namespace IMath
{

/* --- Constants --- */

#ifdef GS_REAL_DOUBLE

using Real = double;

#else

using Real = float;

#endif

}

#endif // IREAL_H
