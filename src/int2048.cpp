#include "int2048.h"

namespace sjtu {

namespace {

constexpr unsigned int MOD1 = 998244353;
constexpr unsigned int MOD2 = 1004535809;
constexpr unsigned int G = 3;
constexpr int MAX_NTT_SIZE = 1 << 19;
constexpr int MAX_STATIC_LIMBS = 131072;

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

static unsigned int knuth_u[MAX_STATIC_LIMBS];
static unsigned int knuth_v[MAX_STATIC_LIMBS];
static unsigned long long schoolbook_c[MAX_STATIC_LIMBS];

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
  size_t num_digits = end - start;
  a.reserve((num_digits + BASE_DIGITS - 1) / BASE_DIGITS);
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
  char limb_buf[9];
  limb_buf[8] = '\0';
  for (int i = (int)a.size() - 2; i >= 0; --i) {
    int cur = a[i];
    for (int d = 7; d >= 0; --d) {
      limb_buf[d] = '0' + (cur % 10);
      cur /= 10;
    }
    std::fputs(limb_buf, stdout);
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
    if (carry >= BASE) {
      res.a.push_back((int)(carry - BASE));
      carry = 1;
    } else {
      res.a.push_back((int)carry);
      carry = 0;
    }
  }
  res.trim();
  return res;
}

int2048 int2048::sub_abs(const int2048 &x, const int2048 &y) {
  int2048 res;
  res.a.resize(x.a.size(), 0);
  int borrow = 0;
  for (size_t i = 0; i < x.a.size(); ++i) {
    int cur = x.a[i] - borrow;
    if (i < y.a.size()) {
      cur -= y.a[i];
    }
    if (cur < 0) {
      res.a[i] = cur + BASE;
      borrow = 1;
    } else {
      res.a[i] = cur;
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
  size_t total_len = n + m;
  std::memset(schoolbook_c, 0, (total_len + 1) * sizeof(unsigned long long));

  for (size_t i = 0; i < n; ++i) {
    unsigned long long xi = (unsigned long long)x.a[i];
    if (xi == 0) {
      continue;
    }
    for (size_t j = 0; j < m; ++j) {
      schoolbook_c[i + j] += xi * (unsigned long long)y.a[j];
    }
  }

  res.a.reserve(total_len + 1);
  unsigned long long carry = 0;
  for (size_t i = 0; i < total_len || carry; ++i) {
    if (i < total_len) {
      carry += schoolbook_c[i];
    }
    res.a.push_back((int)(carry % 100000000ULL));
    carry /= 100000000ULL;
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

  std::memcpy(ntt_a1, va_buf, va_len * sizeof(unsigned int));
  std::memset(ntt_a1 + va_len, 0, (K - va_len) * sizeof(unsigned int));
  std::memcpy(ntt_a2, va_buf, va_len * sizeof(unsigned int));
  std::memset(ntt_a2 + va_len, 0, (K - va_len) * sizeof(unsigned int));

  bool is_square = (&x == &y) || (x.a == y.a);

  if (!is_square) {
    std::memcpy(ntt_b1, vb_buf, vb_len * sizeof(unsigned int));
    std::memset(ntt_b1 + vb_len, 0, (K - vb_len) * sizeof(unsigned int));
    std::memcpy(ntt_b2, vb_buf, vb_len * sizeof(unsigned int));
    std::memset(ntt_b2 + vb_len, 0, (K - vb_len) * sizeof(unsigned int));
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

  std::memset(ntt_res, 0, (need + 10) * sizeof(unsigned long long));
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

  size_t res_len = need + 10;
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
  if (m <= 400) {
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
  unsigned long long uv = (unsigned long long)v;
  for (int i = (int)x.a.size() - 1; i >= 0; --i) {
    rem = rem * 100000000ULL + (unsigned long long)x.a[i];
    q.a[i] = (int)(rem / uv);
    rem %= uv;
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

  unsigned long long d = 100000000ULL / ((unsigned long long)y.a[m - 1] + 1);

  unsigned long long carry = 0;
  for (int i = 0; i < m; ++i) {
    unsigned long long cur = (unsigned long long)y.a[i] * d + carry;
    carry = cur / 100000000ULL;
    knuth_v[i] = (unsigned int)(cur - carry * 100000000ULL);
  }

  carry = 0;
  for (size_t i = 0; i < x.a.size(); ++i) {
    unsigned long long cur = (unsigned long long)x.a[i] * d + carry;
    carry = cur / 100000000ULL;
    knuth_u[i] = (unsigned int)(cur - carry * 100000000ULL);
  }
  knuth_u[x.a.size()] = (unsigned int)carry;

  q.a.assign(n + 1, 0);

  for (int j = n; j >= 0; --j) {
    unsigned long long num = (unsigned long long)knuth_u[j + m] * 100000000ULL + knuth_u[j + m - 1];
    unsigned long long q_hat = num / knuth_v[m - 1];
    unsigned long long r_hat = num - q_hat * knuth_v[m - 1];

    while (q_hat >= 100000000ULL || q_hat * knuth_v[m - 2] > r_hat * 100000000ULL + knuth_u[j + m - 2]) {
      --q_hat;
      r_hat += knuth_v[m - 1];
      if (r_hat >= 100000000ULL) {
        break;
      }
    }

    unsigned long long sub_carry = 0;
    for (int i = 0; i < m; ++i) {
      unsigned long long prod = q_hat * knuth_v[i] + sub_carry;
      sub_carry = prod / 100000000ULL;
      unsigned long long p_low = prod - sub_carry * 100000000ULL;
      if (knuth_u[j + i] < p_low) {
        knuth_u[j + i] += (unsigned int)(100000000ULL - p_low);
        sub_carry++;
      } else {
        knuth_u[j + i] -= (unsigned int)p_low;
      }
    }

    if (knuth_u[j + m] < sub_carry) {
      knuth_u[j + m] += (unsigned int)(100000000ULL - sub_carry);
      --q_hat;
      unsigned long long add_carry = 0;
      for (int i = 0; i < m; ++i) {
        unsigned long long sum = (unsigned long long)knuth_u[j + i] + knuth_v[i] + add_carry;
        add_carry = sum / 100000000ULL;
        knuth_u[j + i] = (unsigned int)(sum - add_carry * 100000000ULL);
      }
      knuth_u[j + m] = (unsigned int)((knuth_u[j + m] + add_carry) % 100000000ULL);
    } else {
      knuth_u[j + m] -= (unsigned int)sub_carry;
    }
    q.a[j] = (int)q_hat;
  }

  q.sgn = 1;
  q.trim();

  r.a.assign(m, 0);
  unsigned long long rem_carry = 0;
  for (int i = m - 1; i >= 0; --i) {
    unsigned long long cur = rem_carry * 100000000ULL + knuth_u[i];
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

  q.a.clear();
  q.a.resize(n - m + 2, 0);
  for (size_t i = 0; i < k; ++i) {
    size_t offset = i * M;
    for (size_t j = 0; j < q_blocks[i].a.size(); ++j) {
      if (offset + j < q.a.size()) {
        q.a[offset + j] += q_blocks[i].a[j];
      } else {
        q.a.push_back(q_blocks[i].a[j]);
      }
    }
  }
  long long carry = 0;
  for (size_t i = 0; i < q.a.size() || carry; ++i) {
    if (i < q.a.size()) {
      carry += q.a[i];
      q.a[i] = (int)(carry % BASE);
      carry /= BASE;
    } else {
      q.a.push_back((int)(carry % BASE));
      carry /= BASE;
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

  size_t n = x.a.size();
  size_t m = y.a.size();

  if (m <= 800 || (n - m) <= 100 || (n - m + 1) * m <= 2000000) {
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
  size_t sa = a.size();
  size_t sb = other.a.size();
  if (sa * sb <= 2048 || sa <= 16 || sb <= 16) {
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
  std::string s;
  s.reserve(x.a.size() * int2048::BASE_DIGITS + 2);
  if (x.sgn == -1) {
    s.push_back('-');
  }
  int high = x.a.back();
  char high_buf[16];
  int high_len = std::snprintf(high_buf, sizeof(high_buf), "%d", high);
  s.append(high_buf, high_len);

  for (int i = (int)x.a.size() - 2; i >= 0; --i) {
    int cur = x.a[i];
    char limb_buf[8];
    for (int d = 7; d >= 0; --d) {
      limb_buf[d] = (char)('0' + (cur % 10));
      cur /= 10;
    }
    s.append(limb_buf, 8);
  }
  os << s;
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
