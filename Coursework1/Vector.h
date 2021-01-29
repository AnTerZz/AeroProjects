//Header file which defines the user-defined Vector class
class Vector {
public:
    Vector() = delete;                  // Do not allow vector without size to be created
    Vector(const unsigned int pSize);   // Create zero vector of given size
    Vector(const Vector& pSrc);         // User-Defined Copy constructor
    Vector(Vector&& pSrc);              // Default move constructor

    double& operator[](unsigned int i); // Access entries in custom vector from this class

    //Operator overloads to allow for vector algebra
    Vector& operator=(const Vector& pSrc); 
    Vector operator+(const Vector& pSrc);
    Vector operator*(const double factor);
    Vector operator/(const double factor);

private:
    double* data = nullptr;           //privately held data and size for each instanciated vector
    unsigned int size = 0;
};
#pragma once //Preprocessor "include guard" replacement
