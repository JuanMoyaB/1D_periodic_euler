#include "FluxFunctions.h"

template<typename T>
EulerFlux<T>::EulerFlux(T gamma_)
{
  gamma = gamma_;
}

template<typename T>
void EulerFlux<T>::computeFlux(DataStruct<T>& rho, DataStruct<T>& rhou, DataStruct<T>& rhoE,
                               DataStruct<T>& flux_rho, DataStruct<T>& flux_rhou, DataStruct<T>& flux_rhoE)
{
  int N = rho.getSize();

  for (int i = 0; i < N; ++i)
  {
    T r = rho[i];
    T ru = rhou[i];
    T rE = rhoE[i];
    T u = ru / r;
    T p = (gamma - 1.0) * (rE - 0.5 * r * u * u);

    flux_rho[i]  = ru;
    flux_rhou[i] = ru * u + p;
    flux_rhoE[i] = u * (rE + p);
  }
}

// Instanciaciones explícitas
template class EulerFlux<float>;
template class EulerFlux<double>;
