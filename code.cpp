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
    static const int BASE = 10000;
    static const int BASE_DIGITS = 4;

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

#include <algorithm>
#include <cmath>
#include <iomanip>

namespace sjtu {

const double int2048::PI = 3.14159265358979323846;

void int2048::trim() {
    while (digits.size() > 1 && digits.back() == 0) {
        digits.pop_back();
    }
    if (digits.empty()) {
        digits.push_back(0);
        neg = false;
    }
    if (digits.size() == 1 && digits[0] == 0) {
        neg = false;
    }
}

bool int2048::absLess(const int2048& other) const {
    if (digits.size() != other.digits.size()) {
        return digits.size() < other.digits.size();
    }
    for (int i = (int)digits.size() - 1; i >= 0; --i) {
        if (digits[i] != other.digits[i]) {
            return digits[i] < other.digits[i];
        }
    }
    return false;
}

bool int2048::absEqual(const int2048& other) const {
    return digits == other.digits;
}

// Assumes both are non-negative
void int2048::absAdd(const int2048& other) {
    int carry = 0;
    size_t n = std::max(digits.size(), other.digits.size());
    digits.resize(n, 0);
    for (size_t i = 0; i < n; ++i) {
        int sum = digits[i] + carry + (i < other.digits.size() ? other.digits[i] : 0);
        digits[i] = sum % BASE;
        carry = sum / BASE;
    }
    if (carry) digits.push_back(carry);
}

// Assumes both are non-negative and |this| >= |other|
void int2048::absSub(const int2048& other) {
    int borrow = 0;
    for (size_t i = 0; i < digits.size(); ++i) {
        int sub = digits[i] - borrow - (i < other.digits.size() ? other.digits[i] : 0);
        if (sub < 0) {
            sub += BASE;
            borrow = 1;
        } else {
            borrow = 0;
        }
        digits[i] = sub;
    }
    trim();
}

void int2048::fft(std::vector<std::complex<double> >& a, bool invert) {
    int n = (int)a.size();
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (int j = 0; j < len / 2; ++j) {
                std::complex<double> u = a[i + j], v = a[i + j] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert) {
        for (auto& x : a) x /= n;
    }
}

std::vector<int> int2048::multiplySmall(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.empty() || b.empty()) return {};
    std::vector<long long> res(a.size() + b.size(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            res[i + j] += (long long)a[i] * b[j];
        }
    }
    std::vector<int> result;
    result.reserve(res.size());
    long long carry = 0;
    for (size_t i = 0; i < res.size(); ++i) {
        long long cur = res[i] + carry;
        result.push_back((int)(cur % BASE));
        carry = cur / BASE;
    }
    while (carry) {
        result.push_back((int)(carry % BASE));
        carry /= BASE;
    }
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    if (result.size() == 1 && result[0] == 0) result.clear();
    return result;
}

std::vector<int> int2048::multiplyLarge(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.empty() || b.empty()) return {};
    std::vector<std::complex<double> > fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < (int)(a.size() + b.size())) n <<= 1;
    fa.resize(n); fb.resize(n);
    fft(fa, false); fft(fb, false);
    for (int i = 0; i < n; ++i) fa[i] *= fb[i];
    fft(fa, true);
    std::vector<long long> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = (long long)std::round(fa[i].real());
    }
    long long carry = 0;
    for (int i = 0; i < n; ++i) {
        long long cur = res[i] + carry;
        res[i] = cur % BASE;
        carry = cur / BASE;
    }
    std::vector<int> result;
    result.reserve(n);
    for (int i = 0; i < n; ++i) {
        result.push_back((int)res[i]);
    }
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    if (result.size() == 1 && result[0] == 0) result.clear();
    return result;
}

std::vector<int> int2048::multiply(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.empty() || b.empty()) return {};
    if (std::min(a.size(), b.size()) < 64) {
        return multiplySmall(a, b);
    }
    return multiplyLarge(a, b);
}

// Check if b * d <= current (both positive, little-endian)
static bool prodLeCurrent(const std::vector<int>& b, int d, const std::vector<int>& current) {
    static const int BASE = int2048::BASE;
    if (d == 0) return true;
    // Compute b * d
    std::vector<long long> prod(b.size() + 1, 0);
    long long c = 0;
    for (size_t j = 0; j < b.size(); ++j) {
        long long p = (long long)b[j] * d + c;
        prod[j] = p % BASE;
        c = p / BASE;
    }
    prod[b.size()] = c;
    // Compare with current
    int psize = (int)prod.size();
    while (psize > 1 && prod[psize-1] == 0) psize--;
    int csize = (int)current.size();
    while (csize > 1 && current[csize-1] == 0) csize--;
    if (psize > csize) return false;
    if (psize < csize) return true;
    for (int j = csize - 1; j >= 0; --j) {
        if (prod[j] > current[j]) return false;
        if (prod[j] < current[j]) return true;
    }
    return true;
}

// Subtract b * d from current. Assumes b * d <= current.
static void subProd(std::vector<int>& current, const std::vector<int>& b, int d) {
    static const int BASE = int2048::BASE;
    if (d == 0) return;
    long long borrow = 0;
    for (size_t j = 0; j < current.size(); ++j) {
        long long sub = current[j] - borrow - (j < b.size() ? (long long)b[j] * d : 0);
        if (sub < 0) {
            long long b2 = (-sub + BASE - 1) / BASE;
            sub += b2 * BASE;
            borrow = b2;
        } else {
            borrow = 0;
        }
        current[j] = (int)sub;
    }
    while (current.size() > 1 && current.back() == 0) current.pop_back();
    if (current.size() == 1 && current[0] == 0) { /* keep single zero */ }
}

// Long division for positive a, b. Computes q = floor(a/b), r = a mod b.
void int2048::longDivide(const std::vector<int>& a, const std::vector<int>& b,
                          std::vector<int>& q, std::vector<int>& r) {
    if (b.empty() || (b.size() == 1 && b[0] == 0)) return;
    if (a.empty()) { q.clear(); r.clear(); r.push_back(0); return; }

    // Compare a and b from most significant digit
    bool aLessB = false;
    if (a.size() < b.size()) aLessB = true;
    else if (a.size() == b.size()) {
        for (int i = (int)a.size() - 1; i >= 0; --i) {
            if (a[i] != b[i]) { aLessB = a[i] < b[i]; break; }
        }
    }
    if (aLessB) {
        q.clear();
        r = a;
        return;
    }

    q.clear();
    std::vector<int> current;

    // Process from most significant digit to least
    for (int i = (int)a.size() - 1; i >= 0; --i) {
        current.insert(current.begin(), a[i]);
        while (current.size() > 1 && current.back() == 0) current.pop_back();

        // Binary search for the largest d in [0, BASE-1] such that b * d <= current
        int lo = 0, hi = BASE - 1;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (prodLeCurrent(b, mid, current)) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        int d = lo;
        q.push_back(d);
        subProd(current, b, d);
    }

    // q was built from MSB to LSB, reverse to get little-endian
    std::reverse(q.begin(), q.end());
    while (q.size() > 1 && q.back() == 0) q.pop_back();
    if (q.size() == 1 && q[0] == 0) q.clear();

    r = current;
    if (r.size() == 1 && r[0] == 0) r.clear();
}

int2048::int2048() : neg(false) {
    digits.push_back(0);
}

int2048::int2048(long long v) : neg(false) {
    if (v < 0) {
        neg = true;
        v = -v;
    }
    if (v == 0) {
        digits.push_back(0);
    } else {
        while (v > 0) {
            digits.push_back((int)(v % BASE));
            v /= BASE;
        }
    }
}

int2048::int2048(const std::string& s) : neg(false) {
    read(s);
}

int2048::int2048(const int2048& other) : digits(other.digits), neg(other.neg) {}

void int2048::read(const std::string& s) {
    digits.clear();
    neg = false;
    size_t start = 0;
    if (!s.empty() && s[0] == '-') {
        neg = true;
        start = 1;
    }
    for (int i = (int)s.size() - 1; i >= (int)start; i -= BASE_DIGITS) {
        int val = 0;
        int base = 1;
        for (int j = 0; j < BASE_DIGITS && i - j >= (int)start; ++j) {
            val += (s[i - j] - '0') * base;
            base *= 10;
        }
        digits.push_back(val);
    }
    trim();
    if (digits.size() == 1 && digits[0] == 0) neg = false;
}

void int2048::print() {
    if (neg) std::cout << '-';
    if (digits.empty()) {
        std::cout << '0';
        return;
    }
    std::cout << digits.back();
    for (int i = (int)digits.size() - 2; i >= 0; --i) {
        std::cout << std::setw(BASE_DIGITS) << std::setfill('0') << digits[i];
    }
}

int2048& int2048::add(const int2048& other) {
    if (!neg && !other.neg) {
        absAdd(other);
    } else if (neg && other.neg) {
        absAdd(other);
    } else if (!neg && other.neg) {
        int2048 tmp = other;
        tmp.neg = false;
        if (absLess(tmp)) {
            tmp.absSub(*this);
            *this = tmp;
            neg = true;
        } else {
            absSub(tmp);
        }
    } else {
        // -|this| + |other|
        int2048 tmp = other;
        tmp.neg = false;
        if (absLess(tmp)) {
            tmp.absSub(*this);
            *this = tmp;
            neg = false;
        } else {
            absSub(tmp);
            neg = true;
        }
    }
    trim();
    return *this;
}

int2048 add(int2048 a, const int2048& b) {
    a.add(b);
    return a;
}

int2048& int2048::minus(const int2048& other) {
    return add(-other);
}

int2048 minus(int2048 a, const int2048& b) {
    a.minus(b);
    return a;
}

int2048 int2048::operator+() const {
    return *this;
}

int2048 int2048::operator-() const {
    int2048 result = *this;
    if (digits.empty() || (digits.size() == 1 && digits[0] == 0)) {
        result.neg = false;
    } else {
        result.neg = !result.neg;
    }
    return result;
}

int2048& int2048::operator=(const int2048& other) {
    digits = other.digits;
    neg = other.neg;
    return *this;
}

int2048& int2048::operator+=(const int2048& other) {
    return add(other);
}

int2048 operator+(int2048 a, const int2048& b) {
    a += b;
    return a;
}

int2048& int2048::operator-=(const int2048& other) {
    return minus(other);
}

int2048 operator-(int2048 a, const int2048& b) {
    a -= b;
    return a;
}

int2048& int2048::operator*=(const int2048& other) {
    if (digits.empty() || (digits.size() == 1 && digits[0] == 0) ||
        other.digits.empty() || (other.digits.size() == 1 && other.digits[0] == 0)) {
        digits.clear();
        digits.push_back(0);
        neg = false;
        return *this;
    }
    bool new_neg = (neg != other.neg);
    digits = multiply(digits, other.digits);
    neg = new_neg;
    trim();
    return *this;
}

int2048 operator*(int2048 a, const int2048& b) {
    a *= b;
    return a;
}

int2048& int2048::operator/=(const int2048& other) {
    if (digits.empty() || (digits.size() == 1 && digits[0] == 0)) {
        return *this;
    }
    bool new_neg = (neg != other.neg);
    std::vector<int> ad = digits;
    std::vector<int> bd = other.digits;
    std::vector<int> q, r;
    longDivide(ad, bd, q, r);
    digits = q;
    neg = new_neg;
    // Adjust for floor division
    bool has_remainder = !r.empty() && !(r.size() == 1 && r[0] == 0);
    if (new_neg && has_remainder) {
        int carry = 1;
        for (size_t i = 0; i < digits.size() && carry; ++i) {
            int sum = digits[i] + carry;
            digits[i] = sum % BASE;
            carry = sum / BASE;
        }
        if (carry) digits.push_back(carry);
    }
    trim();
    return *this;
}

int2048 operator/(int2048 a, const int2048& b) {
    a /= b;
    return a;
}

int2048& int2048::operator%=(const int2048& other) {
    int2048 q = *this;
    q /= other;
    int2048 prod = q * other;
    *this -= prod;
    return *this;
}

int2048 operator%(int2048 a, const int2048& b) {
    a %= b;
    return a;
}

std::istream& operator>>(std::istream& is, int2048& n) {
    std::string s;
    is >> s;
    n.read(s);
    return is;
}

std::ostream& operator<<(std::ostream& os, const int2048& n) {
    if (n.neg) os << '-';
    if (n.digits.empty()) {
        os << '0';
        return os;
    }
    os << n.digits.back();
    for (int i = (int)n.digits.size() - 2; i >= 0; --i) {
        os << std::setw(int2048::BASE_DIGITS) << std::setfill('0') << n.digits[i];
    }
    return os;
}

bool operator==(const int2048& a, const int2048& b) {
    return a.neg == b.neg && a.digits == b.digits;
}

bool operator!=(const int2048& a, const int2048& b) {
    return !(a == b);
}

bool operator<(const int2048& a, const int2048& b) {
    if (a.neg != b.neg) return a.neg;
    if (a.neg) return b.absLess(a);
    return a.absLess(b);
}

bool operator>(const int2048& a, const int2048& b) {
    return b < a;
}

bool operator<=(const int2048& a, const int2048& b) {
    return !(b < a);
}

bool operator>=(const int2048& a, const int2048& b) {
    return !(a < b);
}

} // namespace sjtu
