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

namespace {

constexpr unsigned int MOD1 = 998244353;
constexpr unsigned int MOD2 = 1004535809;
constexpr unsigned int G = 3;
constexpr int MAX_NTT_SIZE = 1 << 19;

template <unsigned int MOD, unsigned int G_VAL>
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

  static const unsigned int *get_roots(bool invert) {
    static bool init = false;
    static unsigned int ws_fwd[MAX_NTT_SIZE];
    static unsigned int ws_inv[MAX_NTT_SIZE];
    if (!init) {
      init = true;
      for (int len = 2; len <= MAX_NTT_SIZE; len <<= 1) {
        int half = len >> 1;
        unsigned int w_fwd = qpow(G_VAL, (MOD - 1) / len);
        unsigned int w_inv = qpow(G_VAL, MOD - 1 - (MOD - 1) / len);
        unsigned int cur_fwd = 1, cur_inv = 1;
        for (int j = 0; j < half; ++j) {
          ws_fwd[half + j] = cur_fwd;
          ws_inv[half + j] = cur_inv;
          cur_fwd = (unsigned long long)cur_fwd * w_fwd % MOD;
          cur_inv = (unsigned long long)cur_inv * w_inv % MOD;
        }
      }
    }
    return invert ? ws_inv : ws_fwd;
  }

  static void ntt(unsigned int *poly, int n, bool invert) {
    static int rev[MAX_NTT_SIZE];
    static int last_n = 0;
    if (n != last_n) {
      last_n = n;
      int bit = 0;
      while ((1 << bit) < n) {
        ++bit;
      }
      for (int i = 0; i < n; ++i) {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
      }
    }

    for (int i = 0; i < n; ++i) {
      if (i < rev[i]) {
        unsigned int tmp = poly[i];
        poly[i] = poly[rev[i]];
        poly[rev[i]] = tmp;
      }
    }

    const unsigned int *ws = get_roots(invert);

    for (int len = 2; len <= n; len <<= 1) {
      int half = len >> 1;
      const unsigned int *w = ws + half;
      for (int i = 0; i < n; i += len) {
        for (int j = 0; j < half; ++j) {
          unsigned int u = poly[i + j];
          unsigned int v = (unsigned long long)poly[i + j + half] * w[j] % MOD;
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

static unsigned int ntt_a1[MAX_NTT_SIZE];
static unsigned int ntt_b1[MAX_NTT_SIZE];
static unsigned int ntt_a2[MAX_NTT_SIZE];
static unsigned int ntt_b2[MAX_NTT_SIZE];
static unsigned long long ntt_res[MAX_NTT_SIZE + 10];
static unsigned int va_buf[MAX_NTT_SIZE];
static unsigned int vb_buf[MAX_NTT_SIZE];

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
    std::printf("%08d", a[i]);
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
  size_t va_len = 0;
  for (size_t i = 0; i < x.a.size(); ++i) {
    int v = x.a[i];
    va_buf[va_len++] = v % 10000;
    va_buf[va_len++] = v / 10000;
  }
  while (va_len > 0 && va_buf[va_len - 1] == 0) {
    --va_len;
  }

  size_t vb_len = 0;
  for (size_t i = 0; i < y.a.size(); ++i) {
    int v = y.a[i];
    vb_buf[vb_len++] = v % 10000;
    vb_buf[vb_len++] = v / 10000;
  }
  while (vb_len > 0 && vb_buf[vb_len - 1] == 0) {
    --vb_len;
  }

  if (va_len == 0 || vb_len == 0) {
    return int2048();
  }

  size_t need = va_len + vb_len - 1;
  size_t K = 1;
  while (K < need) {
    K <<= 1;
  }

  std::memset(ntt_a1, 0, K * sizeof(unsigned int));
  std::memset(ntt_a2, 0, K * sizeof(unsigned int));
  std::memcpy(ntt_a1, va_buf, va_len * sizeof(unsigned int));
  std::memcpy(ntt_a2, va_buf, va_len * sizeof(unsigned int));

  bool is_square = (&x == &y) || (x.a == y.a);

  if (!is_square) {
    std::memset(ntt_b1, 0, K * sizeof(unsigned int));
    std::memset(ntt_b2, 0, K * sizeof(unsigned int));
    std::memcpy(ntt_b1, vb_buf, vb_len * sizeof(unsigned int));
    std::memcpy(ntt_b2, vb_buf, vb_len * sizeof(unsigned int));
  }

  NTTCalculator<MOD1, G>::ntt(ntt_a1, (int)K, false);
  if (is_square) {
    for (size_t i = 0; i < K; ++i) {
      ntt_a1[i] = (unsigned long long)ntt_a1[i] * ntt_a1[i] % MOD1;
    }
  } else {
    NTTCalculator<MOD1, G>::ntt(ntt_b1, (int)K, false);
    for (size_t i = 0; i < K; ++i) {
      ntt_a1[i] = (unsigned long long)ntt_a1[i] * ntt_b1[i] % MOD1;
    }
  }
  NTTCalculator<MOD1, G>::ntt(ntt_a1, (int)K, true);

  NTTCalculator<MOD2, G>::ntt(ntt_a2, (int)K, false);
  if (is_square) {
    for (size_t i = 0; i < K; ++i) {
      ntt_a2[i] = (unsigned long long)ntt_a2[i] * ntt_a2[i] % MOD2;
    }
  } else {
    NTTCalculator<MOD2, G>::ntt(ntt_b2, (int)K, false);
    for (size_t i = 0; i < K; ++i) {
      ntt_a2[i] = (unsigned long long)ntt_a2[i] * ntt_b2[i] % MOD2;
    }
  }
  NTTCalculator<MOD2, G>::ntt(ntt_a2, (int)K, true);

  static const unsigned int inv_mod1_mod2 = NTTCalculator<MOD2, G>::mod_inv(MOD1 % MOD2);

  std::memset(ntt_res, 0, (K + 10) * sizeof(unsigned long long));
  unsigned long long carry = 0;
  for (size_t i = 0; i < need || carry; ++i) {
    if (i < need) {
      unsigned long long r1 = ntt_a1[i];
      unsigned long long r2 = ntt_a2[i];
      unsigned long long diff = (r2 >= r1) ? (r2 - r1) : (r2 + MOD2 - r1);
      unsigned long long k = (diff % MOD2) * inv_mod1_mod2 % MOD2;
      unsigned long long val = r1 + k * (unsigned long long)MOD1;
      carry += val;
    }
    ntt_res[i] = carry % 10000;
    carry /= 10000;
  }

  size_t res_len = K + 10;
  while (res_len > 0 && ntt_res[res_len - 1] == 0) {
    --res_len;
  }

  int2048 res;
  size_t n_chunks = (res_len + 1) / 2;
  res.a.reserve(n_chunks);
  for (size_t i = 0; i < n_chunks; ++i) {
    int chunk = (int)ntt_res[2 * i];
    if (2 * i + 1 < res_len) {
      chunk += (int)ntt_res[2 * i + 1] * 10000;
    }
    res.a.push_back(chunk);
  }
  res.sgn = 1;
  res.trim();
  return res;
}

int2048 int2048::shift_limbs_left(const int2048 &x, size_t k) {
  if (x.a.empty() || k == 0) {
    return x;
  }
  int2048 res;
  res.a.resize(x.a.size() + k, 0);
  for (size_t i = 0; i < x.a.size(); ++i) {
    res.a[i + k] = x.a[i];
  }
  res.sgn = x.sgn;
  return res;
}

int2048 int2048::shift_limbs_right(const int2048 &x, size_t k) {
  if (k >= x.a.size()) {
    return int2048();
  }
  int2048 res;
  res.a.assign(x.a.begin() + k, x.a.end());
  res.sgn = x.sgn;
  res.trim();
  return res;
}

int2048 int2048::compute_reciprocal(const int2048 &B) {
  size_t m = B.a.size();
  if (m <= 64) {
    int2048 num;
    num.a.assign(2 * m + 1, 0);
    num.a[2 * m] = 1;
    num.sgn = 1;
    int2048 q, r;
    div_mod_knuth(num, B, q, r);
    return q;
  }

  size_t k = (m + 2) / 2;
  int2048 B_high = shift_limbs_right(B, m - k);
  int2048 X0 = compute_reciprocal(B_high);

  int2048 T = add_abs(X0, X0);
  T = shift_limbs_left(T, m - k);

  int2048 X0_sq = X0 * X0;
  int2048 Prod = B * X0_sq;
  int2048 S = shift_limbs_right(Prod, 2 * k);

  int2048 X1;
  if (compare_abs(T, S) >= 0) {
    X1 = sub_abs(T, S);
  } else {
    X1 = int2048();
  }

  int2048 num;
  num.a.assign(2 * m + 1, 0);
  num.a[2 * m] = 1;
  num.sgn = 1;

  int2048 P = X1 * B;
  if (compare_abs(num, P) >= 0) {
    int2048 rem = sub_abs(num, P);
    while (compare_abs(rem, B) >= 0) {
      rem = sub_abs(rem, B);
      X1 = add_abs(X1, int2048(1));
    }
  } else {
    int2048 rem = sub_abs(P, num);
    while (!rem.a.empty()) {
      if (compare_abs(rem, B) <= 0) {
        X1 = sub_abs(X1, int2048(1));
        break;
      } else {
        rem = sub_abs(rem, B);
        X1 = sub_abs(X1, int2048(1));
      }
    }
  }
  X1.trim();
  return X1;
}

void int2048::div_mod_small(const int2048 &x, int v, int2048 &q, int2048 &r) {
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
}

void int2048::div_mod_knuth(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r) {
  int m = (int)y.a.size();
  int n = (int)x.a.size() - m;
  if (n < 0) {
    q = int2048();
    r = x;
    r.sgn = 1;
    return;
  }

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

void int2048::div_mod_newton(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r) {
  size_t m = y.a.size();
  size_t n = x.a.size();
  if (n < m) {
    q = int2048();
    r = x;
    r.sgn = 1;
    return;
  }

  int2048 Inv = compute_reciprocal(y);
  size_t M = m;
  size_t k = (n + M - 1) / M;

  int2048 cur_rem;
  std::vector<int2048> q_blocks(k);

  for (int i = (int)k - 1; i >= 0; --i) {
    int2048 block_i;
    size_t start_idx = i * M;
    size_t end_idx = (i + 1) * M < n ? (i + 1) * M : n;
    if (start_idx < end_idx) {
      block_i.a.assign(x.a.begin() + start_idx, x.a.begin() + end_idx);
      block_i.sgn = 1;
      block_i.trim();
    }

    int2048 Cur = add_abs(shift_limbs_left(cur_rem, M), block_i);

    if (compare_abs(Cur, y) < 0) {
      q_blocks[i] = int2048();
      cur_rem = Cur;
    } else {
      int2048 Prod = Cur * Inv;
      int2048 q_est = shift_limbs_right(Prod, 2 * M);

      int2048 q_y = q_est * y;
      if (compare_abs(Cur, q_y) >= 0) {
        int2048 rem = sub_abs(Cur, q_y);
        while (compare_abs(rem, y) >= 0) {
          rem = sub_abs(rem, y);
          q_est = add_abs(q_est, int2048(1));
        }
        cur_rem = rem;
      } else {
        int2048 rem = sub_abs(q_y, Cur);
        while (!rem.a.empty()) {
          if (compare_abs(rem, y) <= 0) {
            q_est = sub_abs(q_est, int2048(1));
            cur_rem = sub_abs(y, rem);
            break;
          } else {
            rem = sub_abs(rem, y);
            q_est = sub_abs(q_est, int2048(1));
          }
        }
      }
      q_blocks[i] = q_est;
    }
  }

  q = int2048();
  for (size_t i = 0; i < k; ++i) {
    if (!q_blocks[i].a.empty()) {
      q = add_abs(q, shift_limbs_left(q_blocks[i], i * M));
    }
  }
  q.sgn = 1;
  q.trim();
  r = cur_rem;
  r.sgn = 1;
  r.trim();
}

void int2048::div_mod_abs(const int2048 &x, const int2048 &y, int2048 &q, int2048 &r) {
  int cmp = compare_abs(x, y);
  if (cmp < 0) {
    q = int2048();
    r = x;
    r.sgn = 1;
    return;
  }
  if (cmp == 0) {
    q = int2048(1);
    r = int2048();
    return;
  }

  if (y.a.size() == 1) {
    div_mod_small(x, y.a[0], q, r);
    return;
  }

  if (y.a.size() <= 64) {
    div_mod_knuth(x, y, q, r);
    return;
  }

  div_mod_newton(x, y, q, r);
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
  if (a.size() * other.a.size() <= 256) {
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
    std::snprintf(buf, sizeof(buf), "%08d", x.a[i]);
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
