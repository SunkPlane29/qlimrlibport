#include <iostream>

extern "C" void greet(const double* a, const double* b, int length, double *out) {
    double sum = 0.0;
    double dot = 0.0;

    std::cout << "Received two arrays of length " << length << ":\n";
    for (int i = 0; i < length; ++i) {
        sum += a[i] + b[i];
        dot += a[i] * b[i];
        std::cout << "  a[" << i << "] = " << a[i] << ", b[" << i << "] = " << b[i] << '\n';
    }

    out[0] = sum;
    out[1] = dot;
}
