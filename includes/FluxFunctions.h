#ifndef _FLUX_FUNCTIONS
#define _FLUX_FUNCTIONS

#include "DataStructs.h"

template<typename T>
class EulerFlux
{
  private:
    T gamma;

  public:
    EulerFlux(T gamma_ = 1.4);  // valor por defecto

    // in: rho, rhou, rhoE; out: flujos
    void computeFlux(DataStruct<T>& rho, DataStruct<T>& rhou, DataStruct<T>& rhoE,
                     DataStruct<T>& flux_rho, DataStruct<T>& flux_rhou, DataStruct<T>& flux_rhoE);
};

#endif // _FLUX_FUNCTIONS
