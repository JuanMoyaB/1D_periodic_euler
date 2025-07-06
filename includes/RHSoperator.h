#ifndef _RHS_OPERATOR
#define _RHS_OPERATOR

#include "DataStructs.h"


template<typename T>
void evaluateEulerRHS(DataStruct<T>& rho, DataStruct<T>& rhou, DataStruct<T>& rhoE,
                      DataStruct<T>& drho_dt, DataStruct<T>& drhou_dt, DataStruct<T>& drhoE_dt,
                      DataStruct<T>& mesh, T gamma);

#endif
