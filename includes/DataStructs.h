#ifndef _DATASTRUCTS
#define _DATASTRUCTS

#include <vector>

template<typename T>
class DataStruct {
    std::vector<T> data;
public:
    // Constructor por defecto
    DataStruct() {}

    // ✅ Constructor que acepta tamaño
    DataStruct(int n) { data.resize(n); }

    void setSize(int n) { data.resize(n); }
    int getSize() const { return data.size(); }

    T* getData() { return data.data(); }
    const T* getData() const { return data.data(); }

    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }
};

// Estructura auxiliar para conservar variables físicas
template<typename T>
struct Conserved {
    T rho;
    T rhou;
    T rhoE;
};

#endif
