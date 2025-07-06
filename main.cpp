#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include "DataStructs.h"

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

const FLOATTYPE gamma_ = 1.4;

Conserved<FLOATTYPE> computeFlux(const Conserved<FLOATTYPE>& U) {
    Conserved<FLOATTYPE> F;
    FLOATTYPE rho = U.rho;
    FLOATTYPE u = U.rhou / rho;
    FLOATTYPE E = U.rhoE / rho;
    FLOATTYPE p = (gamma_ - 1.0) * (U.rhoE - 0.5 * rho * u * u);

    F.rho = U.rhou;
    F.rhou = rho * u * u + p;
    F.rhoE = u * (U.rhoE + p);
    return F;
}

void writeToCSV(const std::string& filename, const DataStruct<FLOATTYPE>& x,
                const DataStruct<Conserved<FLOATTYPE>>& U) {
    std::ofstream file(filename);
    for (int i = 0; i < x.getSize(); ++i) {
        auto u = U[i];
        file << x[i] << "," << u.rho << "," << u.rhou << "," << u.rhoE << "\n";
    }
    file.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        if (rank == 0) std::cout << "Usage: ./euler1d N\n";
        MPI_Finalize();
        return 1;
    }

    int N = std::stoi(argv[1]);
    FLOATTYPE dx = 1.0 / N;
    FLOATTYPE CFL = 0.4;
    FLOATTYPE t_final = 0.2;
    FLOATTYPE time = 0.0;

    DataStruct<FLOATTYPE> x(N);
    DataStruct<Conserved<FLOATTYPE>> U(N), dU(N);

    for (int i = 0; i < N; ++i) {
        x[i] = (i + 0.5) * dx;
        U[i].rho = 1.0 + 0.2 * sin(2.0 * M_PI * x[i]);
        U[i].rhou = 0.0;
        U[i].rhoE = U[i].rho / (gamma_ - 1.0);
    }

    double t0 = MPI_Wtime();
    while (time < t_final) {
        FLOATTYPE dt = CFL * dx;
        if (time + dt > t_final) dt = t_final - time;

        for (int i = 0; i < N; ++i) {
            int ip = (i + 1) % N;
            int im = (i - 1 + N) % N;

            auto Fp = computeFlux(U[ip]);
            auto Fm = computeFlux(U[im]);

            dU[i].rho = -(Fp.rho - Fm.rho) / (2.0 * dx);
            dU[i].rhou = -(Fp.rhou - Fm.rhou) / (2.0 * dx);
            dU[i].rhoE = -(Fp.rhoE - Fm.rhoE) / (2.0 * dx);
        }

        for (int i = 0; i < N; ++i) {
            U[i].rho += dt * dU[i].rho;
            U[i].rhou += dt * dU[i].rhou;
            U[i].rhoE += dt * dU[i].rhoE;
        }

        time += dt;
    }
    double tf = MPI_Wtime() - t0;

    if (rank == 0) {
        std::cout << "Comp. time: " << tf << " sec.\n";
        writeToCSV("solution.csv", x, U);
    }

    MPI_Finalize();
    return 0;
}


