#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
private:
  int sgn;             // 1 for >= 0, -1 for < 0
  std::vector<int> a;  // Digits stored in base 10^8 (little-endian)

  static const int BASE = 100000000;
  static const int BASE_DIGITS = 8;

  void trim();

  static int compare_abs(const int2048 &x, const int2048 &y);
  static int2048 add_abs(const int2048 &x, const int2048 &y);
  static int2048 sub_abs(const int2048 &x, const int2048 &y);
  static int2048 mul_schoolbook(const int2048 &x, const int2048 &y);
  static int2048 mul_ntt(const int2048 &x, const int2048 &y);

  static int2048 shift_limbs_left(const int2048 &x, size_t k);
  static int2048 shift_limbs_right(const int2048 &x, size_t k);
  static int2048 compute_reciprocal(const int2048 &B);

  static void div_mod_small(const int2048 &x, int v, int2048 &q, int2048 &r);
  static void div_mod_knuth(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r);
  static void div_mod_newton(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r);
  static void div_mod_abs(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r);
  static void div_mod(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r);

public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);
  int2048(int2048 &&) noexcept;

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);
  int2048 &operator=(int2048 &&) noexcept;

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

// Non-member declarations in namespace sjtu
int2048 add(int2048, const int2048 &);
int2048 minus(int2048, const int2048 &);
int2048 operator+(int2048, const int2048 &);
int2048 operator-(int2048, const int2048 &);
int2048 operator*(int2048, const int2048 &);
int2048 operator/(int2048, const int2048 &);
int2048 operator%(int2048, const int2048 &);

std::istream &operator>>(std::istream &, int2048 &);
std::ostream &operator<<(std::ostream &, const int2048 &);

bool operator==(const int2048 &, const int2048 &);
bool operator!=(const int2048 &, const int2048 &);
bool operator<(const int2048 &, const int2048 &);
bool operator>(const int2048 &, const int2048 &);
bool operator<=(const int2048 &, const int2048 &);
bool operator>=(const int2048 &, const int2048 &);

} // namespace sjtu

#endif
