#include "rk4.h"
#include <cassert>

template<typename T>
RungeKuttaEuler<T>::RungeKuttaEuler(DataStruct<T> &_rho, DataStruct<T> &_rhou, DataStruct<T> &_rhoE)
: rho(_rho), rhou(_rhou), rhoE(_rhoE)
{
  nSteps = 4;
  currentStep = 0;

  coeffsA = new T[4]{0., 0.5, 0.5, 1.};
  coeffsB = new T[4]{1., 2., 2., 1.};

  frho  = new DataStruct<T>[nSteps];
  frhou = new DataStruct<T>[nSteps];
  frhoE = new DataStruct<T>[nSteps];

  for (int i = 0; i < nSteps; ++i) {
    frho[i].setSize(rho.getSize());
    frhou[i].setSize(rhou.getSize());
    frhoE[i].setSize(rhoE.getSize());
  }

  rho_i.setSize(rho.getSize());
  rhou_i.setSize(rhou.getSize());
  rhoE_i.setSize(rhoE.getSize());
}

template<typename T>
RungeKuttaEuler<T>::~RungeKuttaEuler()
{
  delete[] coeffsA;
  delete[] coeffsB;
  delete[] frho;
  delete[] frhou;
  delete[] frhoE;
}

template<typename T>
int RungeKuttaEuler<T>::getNumSteps() { return nSteps; }

template<typename T>
void RungeKuttaEuler<T>::initRK() { currentStep = 0; }

template<typename T>
void RungeKuttaEuler<T>::stepUi(T dt)
{
  assert(currentStep < nSteps);

  for (int i = 0; i < rho.getSize(); ++i) {
    if (currentStep == 0) {
      rho_i[i]  = rho[i];
      rhou_i[i] = rhou[i];
      rhoE_i[i] = rhoE[i];
    } else {
      rho_i[i]  = rho[i]  + coeffsA[currentStep] * dt * frho[currentStep - 1][i];
      rhou_i[i] = rhou[i] + coeffsA[currentStep] * dt * frhou[currentStep - 1][i];
      rhoE_i[i] = rhoE[i] + coeffsA[currentStep] * dt * frhoE[currentStep - 1][i];
    }
  }
}

template<typename T>
void RungeKuttaEuler<T>::setFi(DataStruct<T> &_frho, DataStruct<T> &_frhou, DataStruct<T> &_frhoE)
{
  for (int i = 0; i < rho.getSize(); ++i) {
    frho[currentStep][i]  = _frho[i];
    frhou[currentStep][i] = _frhou[i];
    frhoE[currentStep][i] = _frhoE[i];
  }

  currentStep++;
}

template<typename T>
void RungeKuttaEuler<T>::finalizeRK(T dt)
{
  for (int i = 0; i < rho.getSize(); ++i) {
    T sumRho  = 0.;
    T sumRhou = 0.;
    T sumRhoE = 0.;

    for (int s = 0; s < nSteps; ++s) {
      sumRho  += coeffsB[s] * frho[s][i];
      sumRhou += coeffsB[s] * frhou[s][i];
      sumRhoE += coeffsB[s] * frhoE[s][i];
    }

    rho[i]  += dt * (1. / 6.) * sumRho;
    rhou[i] += dt * (1. / 6.) * sumRhou;
    rhoE[i] += dt * (1. / 6.) * sumRhoE;
  }
}

template<typename T>
DataStruct<T>* RungeKuttaEuler<T>::currentRho()  { return &rho_i; }

template<typename T>
DataStruct<T>* RungeKuttaEuler<T>::currentRhou() { return &rhou_i; }

template<typename T>
DataStruct<T>* RungeKuttaEuler<T>::currentRhoE() { return &rhoE_i; }

// Instanciaciones
template class RungeKuttaEuler<float>;
template class RungeKuttaEuler<double>;
