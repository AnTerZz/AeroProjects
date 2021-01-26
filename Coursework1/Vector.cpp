#include <cmath>
#include <cstring>
#include <stdexcept>
#include "Vector.h"

using namespace std;

Vector::Vector(const unsigned int pSize) {
    if (pSize == 0) throw std::logic_error("Vector size must be > 0.");

    size = pSize;
    data = new double[pSize];
    memset(data, 0, pSize * sizeof(double));
}

// Copies an existing vector. First we allocate a new block of memory of the
// same size, and then copy the data from pSrc into it.
Vector::Vector(const Vector& pSrc) {
    size = pSrc.size;
    data = new double[pSrc.size];
    memcpy(data, pSrc.data, size * sizeof(double));
}

// "Moves" an existing vector. In this case, we transfer ownership of the
// existing block of data from pSrc to this newly created vector. This avoids
// copying the data. In pSrc, we set the data pointer to nullptr, which is why
// pSrc is not passed as a const parameter.
//
// Note that this is mainly used for efficiency reasons, particular when
// returning temporary objects created in the operator functions below, by
// value.
Vector::Vector(Vector&& pSrc) {
    size = pSrc.size;
    data = pSrc.data;
    pSrc.data = nullptr;
    pSrc.size = 0;
}

// Return a reference to a double so that the calling routine can set values in
// the array, as well as read from them.
double& Vector::operator[](unsigned int i) {
    if (i >= size) throw std::logic_error("Index out of range.");

    return data[i];
}

// Assignment operation. If the vectors are of different sizes, we need to
// deallocate our current storage, allocate new storage of the correct size,
// before finally copying the data across.
Vector& Vector::operator=(const Vector& pSrc) {
    if (pSrc.size != size) {
        delete[] data;
        size = pSrc.size;
        data = new double[size];
    }
    memcpy(data, pSrc.data, size * sizeof(double));
    return *this;
}

// Add two vectors together. First create a new temporary vector to hold the
// result. Perform the addition, element-by-element. Then return the temporary
// object explicitly using move semantics, since the temporary object 'result'
// is about to be destroyed anyway as it goes out of scope.
Vector Vector::operator+(const Vector& pSrc) {
    if (pSrc.size != size) throw std::logic_error("Vectors are different size.");
    Vector result(size);
    for (int i = 0; i < size; ++i) {
        result.data[i] = data[i] + pSrc.data[i];
    }
    return std::move(result);
}


Vector Vector::operator*(const double factor) {
    Vector result(size);
    for (int i = 0; i < size; ++i) {
        result.data[i] = data[i] * factor;
    }
    return std::move(result);
}


Vector Vector::operator/(const double factor) {
    Vector result(size);
    for (int i = 0; i < size; ++i) {
        result.data[i] = data[i] / factor;
    }
    return std::move(result);
}