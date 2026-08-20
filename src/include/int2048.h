#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

namespace sjtu {
class int2048 {
private:
    static const double PI;

    std::vector<int> digits;
    bool neg;

    void trim();
    bool absLess(const int2048& other) const;
    bool absEqual(const int2048& other) const;
    void absAdd(const int2048& other);
    // Assumes |this| >= |other|
    void absSub(const int2048& other);

    static void fft(std::vector<std::complex<double> >& a, bool invert);
    static std::vector<int> multiplyLarge(const std::vector<int>& a, const std::vector<int>& b);
    static std::vector<int> multiplySmall(const std::vector<int>& a, const std::vector<int>& b);
    static std::vector<int> multiply(const std::vector<int>& a, const std::vector<int>& b);

    // Long division with normalization. Assumes a >= b >= 1.
    static void longDivide(const std::vector<int>& a, const std::vector<int>& b,
                           std::vector<int>& q, std::vector<int>& r);

public:
    static const int BASE = 1000;
    static const int BASE_DIGITS = 3;

    int2048();
    int2048(long long);
    int2048(const std::string &);
    int2048(const int2048 &);

    void read(const std::string &);
    void print();

    int2048 &add(const int2048 &);
    friend int2048 add(int2048, const int2048 &);
    int2048 &minus(const int2048 &);
    friend int2048 minus(int2048, const int2048 &);

    int2048 operator+() const;
    int2048 operator-() const;
    int2048 &operator=(const int2048 &);
    int2048 &operator+=(const int2048 &);
    friend int2048 operator+(int2048, const int2048 &);
    int2048 &operator-=(const int2048 &);
    friend int2048 operator-(int2048, const int2048 &);
    int2048 &operator*=(const int2048 &);
    friend int2048 operator*(int2048, const int2048 &);
    int2048 &operator/=(const int2048 &);
    friend int2048 operator/(int2048, const int2048 &);
    int2048 &operator%=(const int2048 &);
    friend int2048 operator%(int2048, const int2048 &);

    friend std::istream &operator>>(std::istream &, int2048 &);
    friend std::ostream &operator<<(std::ostream &, const int2048 &);

    friend bool operator==(const int2048 &, const int2048 &);
    friend bool operator!=(const int2048 &, const int2048 &);
    friend bool operator<(const int2048 &, const int2048 &);
    friend bool operator>(const int2048 &, const int2048 &);
    friend bool operator<=(const int2048 &, const int2048 &);
    friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif
