// R builds package DLLs from top-level src files. Keep implementation files
// organized in subdirectories and compile them through this translation unit.

#include <R_ext/Arith.h>

#undef ISNAN
#undef R_FINITE

#include "distributions/DBTLKJCPC.cc"
#include "distributions/DBTInvGamma.cc"
#include "distributions/DBTMoment.cc"
#include "distributions/DBTInvMoment.cc"
#include "functions/BTLKJCholesky.cc"
#include "invgamma/BTInvGammaCore.cc"
#include "lkj/BTLKJCore.cc"
#include "nonlocal/BTNonlocalCore.cc"
