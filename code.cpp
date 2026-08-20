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
  std::vector<int> a;  // Digits stored in base 10^9 (little-endian)

  static const int BASE = 1000000000;
  static const int BASE_DIGITS = 9;

  void trim();

  static int compare_abs(const int2048 &x, const int2048 &y);
  static int2048 add_abs(const int2048 &x, const int2048 &y);
  static int2048 sub_abs(const int2048 &x, const int2048 &y);
  static int2048 mul_schoolbook(const int2048 &x, const int2048 &y);
  static int2048 mul_ntt(const int2048 &x, const int2048 &y);
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

namespace {

template <unsigned int MOD, unsigned int G>
struct NTTCalculator {
  static unsigned int qpow(unsigned int base, unsigned int exp) {
    unsigned int res = 1;
    while (exp > 0) {
      if (exp & 1) {
        res = (unsigned long long)res * base % MOD;
      }
      base = (unsigned long long)base * base % MOD;
      exp >>= 1;
    }
    return res;
  }

  static unsigned int mod_inv(unsigned int val) {
    return qpow(val, MOD - 2);
  }

  static void ntt(std::vector<unsigned int> &poly, bool invert) {
    int n = (int)poly.size();
    for (int i = 1, j = 0; i < n; ++i) {
      int bit = n >> 1;
      for (; j & bit; bit >>= 1) {
        j ^= bit;
      }
      j ^= bit;
      if (i < j) {
        unsigned int tmp = poly[i];
        poly[i] = poly[j];
        poly[j] = tmp;
      }
    }

    std::vector<unsigned int> w_table(n >> 1);
    for (int len = 2; len <= n; len <<= 1) {
      unsigned int wlen = qpow(G, invert ? (MOD - 1 - (MOD - 1) / len) : ((MOD - 1) / len));
      int half = len >> 1;
      w_table[0] = 1;
      for (int j = 1; j < half; ++j) {
        w_table[j] = (unsigned long long)w_table[j - 1] * wlen % MOD;
      }

      for (int i = 0; i < n; i += len) {
        for (int j = 0; j < half; ++j) {
          unsigned int u = poly[i + j];
          unsigned int v = (unsigned long long)poly[i + j + half] * w_table[j] % MOD;
          poly[i + j] = (u + v >= MOD ? u + v - MOD : u + v);
          poly[i + j + half] = (u < v ? u + MOD - v : u - v);
        }
      }
    }

    if (invert) {
      unsigned int n_inv = mod_inv(n);
      for (int i = 0; i < n; ++i) {
        poly[i] = (unsigned long long)poly[i] * n_inv % MOD;
      }
    }
  }
};

} // namespace

void int2048::trim() {
  while (!a.empty() && a.back() == 0) {
    a.pop_back();
  }
  if (a.empty()) {
    sgn = 1;
  }
}

int2048::int2048() : sgn(1) {}

int2048::int2048(long long v) {
  if (v < 0) {
    sgn = -1;
    unsigned long long uv = 0ULL - (unsigned long long)v;
    while (uv > 0) {
      a.push_back((int)(uv % BASE));
      uv /= BASE;
    }
  } else if (v > 0) {
    sgn = 1;
    unsigned long long uv = (unsigned long long)v;
    while (uv > 0) {
      a.push_back((int)(uv % BASE));
      uv /= BASE;
    }
  } else {
    sgn = 1;
  }
}

int2048::int2048(const std::string &s) {
  read(s);
}

int2048::int2048(const int2048 &other) : sgn(other.sgn), a(other.a) {}

int2048::int2048(int2048 &&other) noexcept : sgn(other.sgn), a(std::move(other.a)) {
  other.sgn = 1;
}

int2048 &int2048::operator=(const int2048 &other) {
  if (this != &other) {
    sgn = other.sgn;
    a = other.a;
  }
  return *this;
}

int2048 &int2048::operator=(int2048 &&other) noexcept {
  if (this != &other) {
    sgn = other.sgn;
    a = std::move(other.a);
    other.sgn = 1;
  }
  return *this;
}

void int2048::read(const std::string &s) {
  a.clear();
  sgn = 1;
  if (s.empty()) {
    return;
  }

  size_t start = 0;
  while (start < s.size() && (s[start] == ' ' || s[start] == '\t' || s[start] == '\n' || s[start] == '\r')) {
    ++start;
  }
  if (start >= s.size()) {
    return;
  }

  int sign = 1;
  if (s[start] == '-') {
    sign = -1;
    ++start;
  } else if (s[start] == '+') {
    sign = 1;
    ++start;
  }

  while (start < s.size() && s[start] == '0') {
    ++start;
  }

  size_t end = s.size();
  while (end > start && (s[end - 1] == ' ' || s[end - 1] == '\t' || s[end - 1] == '\n' || s[end - 1] == '\r')) {
    --end;
  }

  if (start >= end) {
    sgn = 1;
    return;
  }

  sgn = sign;
  for (int i = (int)end; i > (int)start; i -= BASE_DIGITS) {
    int chunk_start = (i - BASE_DIGITS < (int)start) ? (int)start : i - BASE_DIGITS;
    int val = 0;
    for (int j = chunk_start; j < i; ++j) {
      val = val * 10 + (s[j] - '0');
    }
    a.push_back(val);
  }
  trim();
}

void int2048::print() {
  if (a.empty()) {
    std::putchar('0');
    return;
  }
  if (sgn == -1) {
    std::putchar('-');
  }
  std::printf("%d", a.back());
  for (int i = (int)a.size() - 2; i >= 0; --i) {
    std::printf("%09d", a[i]);
  }
}

int int2048::compare_abs(const int2048 &x, const int2048 &y) {
  if (x.a.size() != y.a.size()) {
    return x.a.size() < y.a.size() ? -1 : 1;
  }
  for (int i = (int)x.a.size() - 1; i >= 0; --i) {
    if (x.a[i] != y.a[i]) {
      return x.a[i] < y.a[i] ? -1 : 1;
    }
  }
  return 0;
}

int2048 int2048::add_abs(const int2048 &x, const int2048 &y) {
  int2048 res;
  size_t n = (x.a.size() > y.a.size()) ? x.a.size() : y.a.size();
  res.a.reserve(n + 1);
  long long carry = 0;
  for (size_t i = 0; i < n || carry; ++i) {
    if (i < x.a.size()) {
      carry += x.a[i];
    }
    if (i < y.a.size()) {
      carry += y.a[i];
    }
    res.a.push_back((int)(carry % BASE));
    carry /= BASE;
  }
  res.trim();
  return res;
}

int2048 int2048::sub_abs(const int2048 &x, const int2048 &y) {
  int2048 res;
  res.a.resize(x.a.size(), 0);
  long long borrow = 0;
  for (size_t i = 0; i < x.a.size(); ++i) {
    borrow += x.a[i];
    if (i < y.a.size()) {
      borrow -= y.a[i];
    }
    if (borrow < 0) {
      res.a[i] = (int)(borrow + BASE);
      borrow = -1;
    } else {
      res.a[i] = (int)borrow;
      borrow = 0;
    }
  }
  res.trim();
  return res;
}

int2048 int2048::mul_schoolbook(const int2048 &x, const int2048 &y) {
  int2048 res;
  if (x.a.empty() || y.a.empty()) {
    return res;
  }
  size_t n = x.a.size(), m = y.a.size();
  std::vector<unsigned long long> c(n + m, 0);
  for (size_t i = 0; i < n; ++i) {
    unsigned long long xi = x.a[i];
    if (xi == 0) {
      continue;
    }
    unsigned long long carry = 0;
    for (size_t j = 0; j < m; ++j) {
      unsigned __int128 cur = c[i + j] + (unsigned __int128)xi * y.a[j] + carry;
      c[i + j] = (unsigned long long)(cur % BASE);
      carry = (unsigned long long)(cur / BASE);
    }
    c[i + m] += carry;
  }
  res.a.resize(n + m + 1, 0);
  unsigned long long carry = 0;
  for (size_t i = 0; i < n + m || carry; ++i) {
    if (i < n + m) {
      carry += c[i];
    }
    res.a[i] = (int)(carry % BASE);
    carry /= BASE;
  }
  res.sgn = 1;
  res.trim();
  return res;
}

int2048 int2048::mul_ntt(const int2048 &x, const int2048 &y) {
  std::vector<unsigned int> va;
  va.reserve(x.a.size() * 3);
  for (size_t i = 0; i < x.a.size(); ++i) {
    int v = x.a[i];
    va.push_back(v % 1000);
    v /= 1000;
    va.push_back(v % 1000);
    va.push_back(v / 1000);
  }
  while (!va.empty() && va.back() == 0) {
    va.pop_back();
  }

  std::vector<unsigned int> vb;
  vb.reserve(y.a.size() * 3);
  for (size_t i = 0; i < y.a.size(); ++i) {
    int v = y.a[i];
    vb.push_back(v % 1000);
    v /= 1000;
    vb.push_back(v % 1000);
    vb.push_back(v / 1000);
  }
  while (!vb.empty() && vb.back() == 0) {
    vb.pop_back();
  }

  if (va.empty() || vb.empty()) {
    return int2048();
  }

  size_t need = va.size() + vb.size() - 1;
  size_t K = 1;
  while (K < need) {
    K <<= 1;
  }

  constexpr unsigned int MOD1 = 998244353;
  constexpr unsigned int MOD2 = 1004535809;
  constexpr unsigned int G = 3;

  std::vector<unsigned int> a1(K, 0), b1(K, 0);
  std::vector<unsigned int> a2(K, 0), b2(K, 0);

  for (size_t i = 0; i < va.size(); ++i) {
    a1[i] = a2[i] = va[i];
  }
  for (size_t i = 0; i < vb.size(); ++i) {
    b1[i] = b2[i] = vb[i];
  }

  bool is_square = (&x == &y) || (x.a == y.a);

  NTTCalculator<MOD1, G>::ntt(a1, false);
  if (is_square) {
    for (size_t i = 0; i < K; ++i) {
      a1[i] = (unsigned long long)a1[i] * a1[i] % MOD1;
    }
  } else {
    NTTCalculator<MOD1, G>::ntt(b1, false);
    for (size_t i = 0; i < K; ++i) {
      a1[i] = (unsigned long long)a1[i] * b1[i] % MOD1;
    }
  }
  NTTCalculator<MOD1, G>::ntt(a1, true);

  NTTCalculator<MOD2, G>::ntt(a2, false);
  if (is_square) {
    for (size_t i = 0; i < K; ++i) {
      a2[i] = (unsigned long long)a2[i] * a2[i] % MOD2;
    }
  } else {
    NTTCalculator<MOD2, G>::ntt(b2, false);
    for (size_t i = 0; i < K; ++i) {
      a2[i] = (unsigned long long)a2[i] * b2[i] % MOD2;
    }
  }
  NTTCalculator<MOD2, G>::ntt(a2, true);

  unsigned int inv_mod1_mod2 = NTTCalculator<MOD2, G>::mod_inv(MOD1 % MOD2);

  std::vector<unsigned long long> base10_3(K + 10, 0);
  unsigned long long carry = 0;
  for (size_t i = 0; i < need || carry; ++i) {
    if (i < need) {
      unsigned long long r1 = a1[i];
      unsigned long long r2 = a2[i];
      unsigned long long diff = (r2 >= r1) ? (r2 - r1) : (r2 + MOD2 - r1);
      unsigned long long k = (diff % MOD2) * inv_mod1_mod2 % MOD2;
      unsigned long long val = r1 + k * (unsigned long long)MOD1;
      carry += val;
    }
    base10_3[i] = carry % 1000;
    carry /= 1000;
  }

  size_t base10_3_len = K + 10;
  while (base10_3_len > 0 && base10_3[base10_3_len - 1] == 0) {
    --base10_3_len;
  }

  int2048 res;
  size_t n_chunks = (base10_3_len + 2) / 3;
  res.a.reserve(n_chunks);
  for (size_t i = 0; i < n_chunks; ++i) {
    int chunk = 0;
    if (3 * i < base10_3_len) {
      chunk += (int)base10_3[3 * i];
    }
    if (3 * i + 1 < base10_3_len) {
      chunk += (int)base10_3[3 * i + 1] * 1000;
    }
    if (3 * i + 2 < base10_3_len) {
      chunk += (int)base10_3[3 * i + 2] * 1000000;
    }
    res.a.push_back(chunk);
  }
  res.sgn = 1;
  res.trim();
  return res;
}

void int2048::div_mod_abs(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r) {
  if (compare_abs(x, y) < 0) {
    q = int2048();
    r = x;
    r.sgn = 1;
    return;
  }

  if (y.a.size() == 1) {
    int v = y.a[0];
    q.a.resize(x.a.size(), 0);
    unsigned long long rem = 0;
    for (int i = (int)x.a.size() - 1; i >= 0; --i) {
      rem = rem * BASE + x.a[i];
      q.a[i] = (int)(rem / v);
      rem %= v;
    }
    q.sgn = 1;
    q.trim();
    r.a.clear();
    if (rem > 0) {
      r.a.push_back((int)rem);
    }
    r.sgn = 1;
    return;
  }

  int m = (int)y.a.size();
  int n = (int)x.a.size() - m;
  int d = BASE / (y.a[m - 1] + 1);

  std::vector<unsigned long long> u(x.a.size() + 1, 0);
  std::vector<unsigned long long> v(m, 0);

  unsigned long long carry = 0;
  for (int i = 0; i < m; ++i) {
    unsigned long long cur = (unsigned long long)y.a[i] * d + carry;
    v[i] = cur % BASE;
    carry = cur / BASE;
  }

  carry = 0;
  for (size_t i = 0; i < x.a.size(); ++i) {
    unsigned long long cur = (unsigned long long)x.a[i] * d + carry;
    u[i] = cur % BASE;
    carry = cur / BASE;
  }
  u[x.a.size()] = carry;

  q.a.resize(n + 1, 0);

  for (int j = n; j >= 0; --j) {
    unsigned __int128 num = ((unsigned __int128)u[j + m] * BASE) + u[j + m - 1];
    unsigned long long q_hat = (unsigned long long)(num / v[m - 1]);
    unsigned long long r_hat = (unsigned long long)(num % v[m - 1]);

    while (q_hat >= BASE || (unsigned __int128)q_hat * v[m - 2] > ((unsigned __int128)r_hat * BASE + u[j + m - 2])) {
      --q_hat;
      r_hat += v[m - 1];
      if (r_hat >= BASE) {
        break;
      }
    }

    unsigned long long sub_carry = 0;
    for (int i = 0; i < m; ++i) {
      unsigned __int128 prod = (unsigned __int128)q_hat * v[i] + sub_carry;
      unsigned long long p_low = (unsigned long long)(prod % BASE);
      sub_carry = (unsigned long long)(prod / BASE);
      if (u[j + i] < p_low) {
        u[j + i] += BASE - p_low;
        sub_carry++;
      } else {
        u[j + i] -= p_low;
      }
    }

    if (u[j + m] < sub_carry) {
      u[j + m] += BASE - sub_carry;
      --q_hat;
      unsigned long long add_carry = 0;
      for (int i = 0; i < m; ++i) {
        unsigned long long sum = u[j + i] + v[i] + add_carry;
        u[j + i] = sum % BASE;
        add_carry = sum / BASE;
      }
      u[j + m] = (u[j + m] + add_carry) % BASE;
    } else {
      u[j + m] -= sub_carry;
    }
    q.a[j] = (int)q_hat;
  }

  q.sgn = 1;
  q.trim();

  r.a.resize(m, 0);
  unsigned long long rem_carry = 0;
  for (int i = m - 1; i >= 0; --i) {
    unsigned long long cur = rem_carry * BASE + u[i];
    r.a[i] = (int)(cur / d);
    rem_carry = cur % d;
  }
  r.sgn = 1;
  r.trim();
}

void int2048::div_mod(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r) {
  if (x.a.empty()) {
    q = int2048();
    r = int2048();
    return;
  }

  int2048 q_abs, r_abs;
  div_mod_abs(x, y, q_abs, r_abs);

  if (r_abs.a.empty()) {
    q = q_abs;
    q.sgn = (x.sgn == y.sgn) ? 1 : -1;
    q.trim();
    r = int2048();
    return;
  }

  if (x.sgn == y.sgn) {
    q = q_abs;
    q.sgn = 1;
    q.trim();
    r = r_abs;
    r.sgn = y.sgn;
    r.trim();
  } else {
    int2048 one(1);
    q = add_abs(q_abs, one);
    q.sgn = -1;
    q.trim();

    r = sub_abs(y, r_abs);
    r.sgn = y.sgn;
    r.trim();
  }
}

int2048 int2048::operator+() const {
  return *this;
}

int2048 int2048::operator-() const {
  if (a.empty()) {
    return *this;
  }
  int2048 res(*this);
  res.sgn = -res.sgn;
  return res;
}

int2048 &int2048::operator+=(const int2048 &other) {
  if (other.a.empty()) {
    return *this;
  }
  if (a.empty()) {
    *this = other;
    return *this;
  }
  if (sgn == other.sgn) {
    int s = sgn;
    *this = add_abs(*this, other);
    sgn = s;
  } else {
    int cmp = compare_abs(*this, other);
    if (cmp == 0) {
      a.clear();
      sgn = 1;
    } else if (cmp > 0) {
      int s = sgn;
      *this = sub_abs(*this, other);
      sgn = s;
    } else {
      int s = other.sgn;
      *this = sub_abs(other, *this);
      sgn = s;
    }
  }
  trim();
  return *this;
}

int2048 operator+(int2048 a, const int2048 &b) {
  return a += b;
}

int2048 &int2048::operator-=(const int2048 &other) {
  if (other.a.empty()) {
    return *this;
  }
  if (a.empty()) {
    *this = -other;
    return *this;
  }
  if (sgn != other.sgn) {
    int s = sgn;
    *this = add_abs(*this, other);
    sgn = s;
  } else {
    int cmp = compare_abs(*this, other);
    if (cmp == 0) {
      a.clear();
      sgn = 1;
    } else if (cmp > 0) {
      int s = sgn;
      *this = sub_abs(*this, other);
      sgn = s;
    } else {
      int s = -sgn;
      *this = sub_abs(other, *this);
      sgn = s;
    }
  }
  trim();
  return *this;
}

int2048 operator-(int2048 a, const int2048 &b) {
  return a -= b;
}

int2048 &int2048::operator*=(const int2048 &other) {
  if (a.empty() || other.a.empty()) {
    a.clear();
    sgn = 1;
    return *this;
  }
  int target_sgn = (sgn == other.sgn) ? 1 : -1;
  if (a.size() * other.a.size() <= 1024) {
    *this = mul_schoolbook(*this, other);
  } else {
    *this = mul_ntt(*this, other);
  }
  sgn = target_sgn;
  trim();
  return *this;
}

int2048 operator*(int2048 a, const int2048 &b) {
  return a *= b;
}

int2048 &int2048::operator/=(const int2048 &other) {
  int2048 q, r;
  div_mod(*this, other, q, r);
  *this = std::move(q);
  return *this;
}

int2048 operator/(int2048 a, const int2048 &b) {
  return a /= b;
}

int2048 &int2048::operator%=(const int2048 &other) {
  int2048 q, r;
  div_mod(*this, other, q, r);
  *this = std::move(r);
  return *this;
}

int2048 operator%(int2048 a, const int2048 &b) {
  return a %= b;
}

int2048 &int2048::add(const int2048 &other) {
  return *this += other;
}

int2048 add(int2048 a, const int2048 &b) {
  return a += b;
}

int2048 &int2048::minus(const int2048 &other) {
  return *this -= other;
}

int2048 minus(int2048 a, const int2048 &b) {
  return a -= b;
}

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s;
  if (is >> s) {
    x.read(s);
  }
  return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.a.empty()) {
    os << 0;
    return os;
  }
  if (x.sgn == -1) {
    os << '-';
  }
  os << x.a.back();
  char buf[16];
  for (int i = (int)x.a.size() - 2; i >= 0; --i) {
    std::snprintf(buf, sizeof(buf), "%09d", x.a[i]);
    os << buf;
  }
  return os;
}

bool operator==(const int2048 &x, const int2048 &y) {
  if (x.a.empty() && y.a.empty()) {
    return true;
  }
  if (x.sgn != y.sgn || x.a.size() != y.a.size()) {
    return false;
  }
  for (size_t i = 0; i < x.a.size(); ++i) {
    if (x.a[i] != y.a[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const int2048 &x, const int2048 &y) {
  return !(x == y);
}

bool operator<(const int2048 &x, const int2048 &y) {
  if (x.a.empty() && y.a.empty()) {
    return false;
  }
  if (x.sgn != y.sgn) {
    return x.sgn < y.sgn;
  }
  int cmp = int2048::compare_abs(x, y);
  return (x.sgn == 1) ? (cmp < 0) : (cmp > 0);
}

bool operator>(const int2048 &x, const int2048 &y) {
  return y < x;
}

bool operator<=(const int2048 &x, const int2048 &y) {
  return !(y < x);
}

bool operator>=(const int2048 &x, const int2048 &y) {
  return !(x < y);
}

} // namespace sjtu

#endif
