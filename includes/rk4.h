#ifndef _RUNGE_KUTTA
#define _RUNGE_KUTTA

#include "DataStructs.h"

template<typename T>
class RungeKuttaEuler
{
private:
  int nSteps;
  int currentStep;

  T *coeffsA, *coeffsB;

  // referencia a las soluciones actuales
  DataStruct<T> &rho, &rhou, &rhoE;

  // soluciones intermedias
  DataStruct<T> rho_i, rhou_i, rhoE_i;

  // RHS por paso
  DataStruct<T> *frho, *frhou, *frhoE;

public:
  RungeKuttaEuler(DataStruct<T> &_rho, DataStruct<T> &_rhou, DataStruct<T> &_rhoE);
  ~RungeKuttaEuler();

  int getNumSteps();
  void initRK();
  void stepUi(T dt);
  void setFi(DataStruct<T> &_frho, DataStruct<T> &_frhou, DataStruct<T> &_frhoE);
  void finalizeRK(T dt);

  DataStruct<T>* currentRho();
  DataStruct<T>* currentRhou();
  DataStruct<T>* currentRhoE();
};

#endif
