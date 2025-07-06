#include "RHSoperator.h"
#include "FluxFunctions.h"

template<typename T>
void evaluateEulerRHS(DataStruct<T>& rho, DataStruct<T>& rhou, DataStruct<T>& rhoE,
                      DataStruct<T>& drho_dt, DataStruct<T>& drhou_dt, DataStruct<T>& drhoE_dt,
                      DataStruct<T>& mesh, T gamma)
{
  int N = rho.getSize();

  DataStruct<T> flux_rho, flux_rhou, flux_rhoE;
  flux_rho.setSize(N);
  flux_rhou.setSize(N);
  flux_rhoE.setSize(N);

  EulerFlux<T> flux_func(gamma);
  flux_func.computeFlux(rho, rhou, rhoE, flux_rho, flux_rhou, flux_rhoE);

  const T* x = mesh.getData();

  for (int i = 0; i < N; ++i)
  {
    int ip = (i + 1) % N;
    int im = (i - 1 + N) % N;

    T dx = x[ip] - x[im];

    drho_dt[i]  = -(flux_rho[ip]  - flux_rho[im])  / dx;
    drhou_dt[i] = -(flux_rhou[ip] - flux_rhou[im]) / dx;
    drhoE_dt[i] = -(flux_rhoE[ip] - flux_rhoE[im]) / dx;
  }
}

// Instanciación explícita
template void evaluateEulerRHS(DataStruct<float>&, DataStruct<float>&, DataStruct<float>&,
                               DataStruct<float>&, DataStruct<float>&, DataStruct<float>&,
                               DataStruct<float>&, float);

template void evaluateEulerRHS(DataStruct<double>&, DataStruct<double>&, DataStruct<double>&,
                               DataStruct<double>&, DataStruct<double>&, DataStruct<double>&,
                               DataStruct<double>&, double);
