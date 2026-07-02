// circulant_mt_fast.cpp
//
// Multithreaded exact search for optimal circulant graphs C(N; 1, s2, ..., sk).
// Optimality criterion: minimize average shortest path length (sum of distances)
// first, then minimize diameter as a tie-breaker. Arbitrary degree k >= 2
// (up to KMAX = 16). All optimal graphs are reported (every isomorphic copy).
//
// Build:
//   g++ -O3 -march=native -funroll-loops -std=c++17 -pthread circulant_mt.cpp -o circulant_mt
//
// Optional profile-guided build (typically 10-15% faster):
//   g++ -O3 -march=native -fprofile-generate -std=c++17 -pthread circulant_mt.cpp -o cm_pgo
//   ./cm_pgo -n 300 -k 5 -t 16 -o /tmp/prof
//   g++ -O3 -march=native -fprofile-use -std=c++17 -pthread circulant_mt.cpp -o circulant_mt
//
// Usage:
//   circulant_mt -n N|a-b -k K|a-b [-m 0|1] [-t threads] [-o dir | -f file]
//     -n   number of vertices, single value or range
//     -k   degree parameter (half-degree), single value or range
//     -m   lower-bound mode: 1 = tight sphere-capacity bound (default), 0 = simple bound
//     -t   thread count (default: hardware concurrency)
//     -o   output directory, one file {N}_{k}.csv per pair
//     -f   single output file for the whole range
//   With neither -o nor -f the result is printed to stdout.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <limits>
#include <chrono>
#include <algorithm>
#include <set>
#include <numeric>
#include <filesystem>

using u32 = uint32_t;
using u64 = uint64_t;

static const int KMAX = 16;          // maximum supported degree parameter
static const int REREAD = 256;       // how often a worker refreshes the shared record

// Modular inverse via extended Euclid; returns -1 if a is not invertible mod n.
static inline int mod_inv(long long a, long long n) {
    long long t = 0, nt = 1, r = n, nr = ((a % n) + n) % n;
    while (nr) { long long q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp; }
    if (r != 1) return -1;
    if (t < 0) t += n;
    return (int)t;
}

// Barrett reduction used inside the canonicality test: replaces the modulo in
// sig[j] * m mod N by a multiply and a shift, avoiding integer division on the
// hot path. mu = floor(2^64 / N). Inputs satisfy x = sig*m < N^2 < 2^64 for N < 2^32.
struct Barrett {
    u64 N, mu;
    void init(u64 n) { N = n; mu = (n ? (~0ULL) / n : 0); }
    inline u32 fold(u64 x) const {
        unsigned __int128 t = (unsigned __int128)x * mu;
        u64 q = (u64)(t >> 64);
        u64 r = x - q * N;
        if (r >= N) r -= N;                 // single correction is sufficient
        u64 nr = N - r;
        return (u32)(nr < r ? nr : r);      // fold to min(r, N - r)
    }
};

// Number of lattice points at L1-distance t in k dimensions.
static u64 sphere_cap_formula(int k, int t) {
    if (t == 0) return 1;
    u64 s = 0, ck = 1, ct = 1; int top = std::min(k, t);
    for (int i = 1; i <= top; ++i) {
        ck = ck * (u64)(k - i + 1) / (u64)i;
        if (i == 1) ct = 1; else ct = ct * (u64)(t - i + 1) / (u64)(i - 1);
        s += ((u64)1 << i) * ck * ct;
    }
    return s;
}

// Pack (diameter, distance sum) into one 64-bit key so the record is a single atomic.
static inline u64 pack(u32 diam, u64 ds) { return (ds << 16) | (diam & 0xFFFFu); }
static inline u32 unpack_diam(u64 key) { return (u32)(key & 0xFFFFu); }
static inline u64 unpack_ds(u64 key) { return key >> 16; }

// Canonicality test in k-1 normalizations.
// Two signatures are equivalent if one is obtained from the other by multiplying
// every generator by an invertible element of Z_N. We normalize by m = inv[sig[gi]]
// for gi = 1..k-1; gi = 0 gives sig[0] = 1 -> m = inv[1] = 1, the identity, skipped.
// The N-g sign case is redundant because fold(-x) = fold(x) already accounts for
// the reflection. A signature is canonical if no normalization yields a smaller tuple.
static bool is_canonical(int N, const int* __restrict__ sig, int k, const int* __restrict__ inv, const Barrett& br) {
    int img[KMAX];
    for (int gi = 1; gi < k; ++gi) {
        int g = sig[gi];
        int m = inv[g];
        if (m == 0) continue;              // not invertible
        // Fast path. The image always contains 1 (fold(sig[gi]*m) = 1) and the
        // signature starts with 1, so the lexicographic comparison is decided by the
        // second smallest element: m2 = min over j != gi of fold(sig[j]*m), versus sig[1].
        // The full sort is needed only when they are equal.
        int m2 = N;
        for (int j = 0; j < k; ++j) {
            if (j == gi) { img[j] = 1; continue; }
            int v = (int)br.fold((u64)sig[j] * (u64)m);
            img[j] = v;
            if (v < m2) m2 = v;
        }
        if (m2 < sig[1]) return false;     // image strictly smaller
        if (m2 > sig[1]) continue;         // image strictly larger, normalization passed
        for (int a = 1; a < k; ++a) { int v = img[a], b = a - 1;
            while (b >= 0 && img[b] > v) { img[b + 1] = img[b]; --b; } img[b + 1] = v; }
        for (int j = 0; j < k; ++j) { if (img[j] < sig[j]) return false; if (img[j] > sig[j]) break; }
    }
    return true;
}

// Pruned BFS over Z_N. State is packed into a single array vals[] instead of a pair
// (dist[], seen[]): vals[v] = (epoch << DIST_BITS) | dist. A vertex is visited iff the
// high bits equal the current epoch tag; the distance is the low DIST_BITS bits. This
// halves the random memory traffic of marking (one store instead of two) without
// changing the traversal order or the pruning schedule. Central symmetry d(v)=d(N-v)
// lets one probe mark a pair of vertices. At each layer boundary a lower bound on the
// distance sum is compared against the record; the search aborts as soon as it cannot win.
// Limitation: diameter < 2^DIST_BITS (= 1024), ample for the practical range of N.
static const int DIST_BITS = 10;
static const u32 DIST_MASK = (1u << DIST_BITS) - 1;

static bool bfs(int N, int k, const int* __restrict__ s, int mode, u64 best_ds,
                u32* __restrict__ vals, u32 epoch, int* __restrict__ queue,
                const u64* caps, int capsN, u32& out_diam, u64& out_ds) {
    const u32 tag = epoch << DIST_BITS;
    vals[0] = tag; queue[0] = 0;
    u32 qr = 0, qw = 1, vc = 1, curd = 0;
    u64 dsum = 0;
    while (vc < (u32)N && qr != qw) {
        const int u = queue[qr++];
        const u32 d = (vals[u] & DIST_MASK) + 1;
        if (curd != d) {
            curd = d;
            u64 lb;
            if (mode == 0) {
                lb = dsum + (u64)(N - vc) * d;
            } else {
                u64 rem = (u64)(N - vc), extra = 0; int t = (int)d;
                while (rem) {
                    u64 cap = (t < capsN) ? caps[t] : sphere_cap_formula(k, t);
                    u64 take = rem < cap ? rem : cap; extra += (u64)t * take; rem -= take; ++t;
                }
                lb = dsum + extra;
            }
            if (lb > best_ds) return false;
        }
        const u32 wd = tag | d;
        for (int i = 0; i < k; ++i) {
            const int sj = s[i];
            int a = u + sj; if (a >= N) a -= N;
            int b = u - sj; if (b < 0) b += N;
            if ((vals[a] & ~DIST_MASK) != tag) { int ma = N - a; vals[a] = wd;
                if (ma != a) vals[ma] = wd;
                queue[qw++] = a; int add = (ma != a) ? 2 : 1; vc += add; dsum += (u64)d * add; }
            if ((vals[b] & ~DIST_MASK) != tag) { int mb = N - b; vals[b] = wd;
                if (mb != b) vals[mb] = wd;
                queue[qw++] = b; int add = (mb != b) ? 2 : 1; vc += add; dsum += (u64)d * add; }
        }
    }
    if (vc < (u32)N) return false;
    out_diam = curd; out_ds = dsum;
    return true;
}

// Moore lower bound on the diameter: smallest D whose cumulative sphere capacity >= N.
static int moore_diam(int N, const u64* caps, int capsN, int k) {
    u64 cum = 1; int t = 0;
    while (cum < (u64)N) { ++t; u64 c = (t < capsN) ? caps[t] : sphere_cap_formula(k, t); cum += c; }
    return t;
}

// Shared best record, cache-line aligned and padded to avoid false sharing.
struct alignas(64) Best { std::atomic<u64> key; char pad[64 - sizeof(std::atomic<u64>)]; Best() { key.store(~0ULL); } };
struct Found { std::vector<int> sig; u32 diam; u64 ds; };

static void run_one(int N, int k, int mode, int threads, std::FILE* out) {
    if (k < 2 || k > KMAX || k > N / 2) return;
    int h = N / 2;

    // precomputation
    std::vector<int> inv(N, 0);
    for (int i = 1; i < N; ++i) { int m = mod_inv(i, N); if (m > 0) inv[i] = m; }
    Barrett br; br.init((u64)N);
    std::vector<u64> caps; { u64 cum = 1; int t = 0; caps.push_back(1);
        while (cum < (u64)N + 1) { ++t; u64 c = sphere_cap_formula(k, t); caps.push_back(c); cum += c; } }
    int capsN = (int)caps.size();
    int Dmoore = moore_diam(N, caps.data(), capsN, k);
    (void)Dmoore;

    Best best;
    std::vector<std::vector<Found>> tloc(threads);

    // Work distribution. For k >= 3 threads pull linearized (s2, s3) pairs from an
    // atomic counter, which gives finer granularity and better load balance than
    // distributing by s2 alone. For k == 2 there is only s2.
    std::atomic<long long> next_pair{ 0 };
    std::atomic<int> next_s2{ 2 };

    std::vector<std::pair<int,int>> pairs;
    if (k >= 3) {
        for (int s2 = 2; s2 <= h - (k - 2); ++s2)
            for (int s3 = s2 + 1; s3 <= h - (k - 3); ++s3)
                pairs.emplace_back(s2, s3);
    }
    long long npairs = (long long)pairs.size();

    auto worker = [&](int tid) {
        std::vector<u32> vals(N, 0);      // packed (epoch << DIST_BITS) | dist
        u32 epoch = 0;
        std::vector<int> queue(N);
        int sig[KMAX]; sig[0] = 1;
        std::vector<Found>& mine = tloc[tid];

        u64 local_best = best.key.load(std::memory_order_relaxed);
        int since = 0;

        auto eval = [&]() {
            if (!is_canonical(N, sig, k, inv.data(), br)) return;
            if (++since >= REREAD) { local_best = best.key.load(std::memory_order_relaxed); since = 0; }
            u32 D; u64 ds;
            // epoch occupies 32 - DIST_BITS bits; clear the array on wraparound
            if (++epoch == (1u << (32 - DIST_BITS))) { std::fill(vals.begin(), vals.end(), 0u); epoch = 1; }
            if (bfs(N, k, sig, mode, unpack_ds(local_best), vals.data(), epoch,
                    queue.data(), caps.data(), capsN, D, ds)) {
                u64 nk = pack(D, ds);
                u64 cur = best.key.load(std::memory_order_relaxed);
                if (nk < cur) {
                    while (nk < cur && !best.key.compare_exchange_weak(cur, nk,
                                std::memory_order_relaxed)) {}
                    local_best = nk;
                    mine.clear(); mine.push_back({ std::vector<int>(sig, sig + k), D, ds });
                } else if (nk == cur) {
                    mine.push_back({ std::vector<int>(sig, sig + k), D, ds });
                }
            }
        };

        if (k == 2) {
            while (true) {
                int s2 = next_s2.fetch_add(1, std::memory_order_relaxed);
                if (s2 > h - (k - 2)) break;
                sig[1] = s2; eval();
            }
            return;
        }

        // k >= 3: take an (s2, s3) pair, then iterate the odometer over positions 4..k-1
        while (true) {
            long long pi = next_pair.fetch_add(1, std::memory_order_relaxed);
            if (pi >= npairs) break;
            int s2 = pairs[pi].first, s3 = pairs[pi].second;
            sig[1] = s2; sig[2] = s3;

            if (k == 3) { eval(); continue; }

            int idx[KMAX];
            for (int j = 3; j < k; ++j) idx[j] = s3 + (j - 2);
            bool done = (idx[k - 1] > h);
            while (!done) {
                for (int j = 3; j < k; ++j) sig[j] = idx[j];
                eval();
                int p = k - 1;
                while (p >= 3) {
                    if (idx[p] < h - (k - 1 - p)) { ++idx[p];
                        for (int q = p + 1; q < k; ++q) idx[q] = idx[q - 1] + 1;
                        break; }
                    --p;
                }
                if (p < 3) done = true;
            }
        }
    };

    std::vector<std::thread> pool;
    for (int i = 0; i < threads; ++i) pool.emplace_back(worker, i);
    for (auto& t : pool) t.join();

    u64 bk = best.key.load();
    std::vector<Found> all;
    for (auto& v : tloc) for (auto& f : v)
        if (pack(f.diam, f.ds) == bk) all.push_back(f);
    std::sort(all.begin(), all.end(), [](const Found& a, const Found& b) { return a.sig < b.sig; });
    all.erase(std::unique(all.begin(), all.end(),
              [](const Found& a, const Found& b) { return a.sig == b.sig; }), all.end());

    // expand every optimal class to all of its isomorphic copies, sorted lexicographically
    std::set<std::vector<int>> allsigs;
    for (auto& f : all) {
        for (int g : f.sig) {
            int m = inv[g]; if (m == 0) continue;
            std::vector<int> img; for (int x : f.sig) img.push_back((int)br.fold((u64)x * (u64)m));
            std::sort(img.begin(), img.end()); allsigs.insert(img);
        }
    }
    u32 D = all.empty() ? 0 : all.front().diam;
    u64 ds = all.empty() ? 0 : all.front().ds;
    double aspl = (double)ds / (N - 1);
    for (const auto& r : allsigs) {
        int deg = 2 * k; for (int x : r) if (2 * x == N) { deg -= 1; break; }
        long long edges = (long long)N * deg / 2;
        std::string s = "C(" + std::to_string(N) + ";1";
        for (int j = 1; j < (int)r.size(); ++j) s += ";" + std::to_string(r[j]);
        s += ")";
        std::fprintf(out, "%d,%d,%s,%u,%.17g,%lld\n", N, k, s.c_str(), D, aspl, edges);
    }
    std::fflush(out);
}

static bool parse_range(const std::string& s, int& a, int& b) {
    auto d = s.find('-');
    if (d == std::string::npos) { a = b = std::atoi(s.c_str()); return a > 0; }
    a = std::atoi(s.substr(0, d).c_str()); b = std::atoi(s.substr(d + 1).c_str());
    return a > 0 && b >= a;
}

int main(int argc, char** argv) {
    int ns = 0, ne = 0, ks = 0, ke = 0, mode = 1;
    int threads = (int)std::thread::hardware_concurrency();
    if (threads < 1) threads = 1;
    std::string outdir, outfile;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&]() { return (i + 1 < argc) ? argv[++i] : ""; };
        if (a == "-n") { if (!parse_range(need(), ns, ne)) { fprintf(stderr, "bad -n\n"); return 1; } }
        else if (a == "-k") { if (!parse_range(need(), ks, ke)) { fprintf(stderr, "bad -k\n"); return 1; } }
        else if (a == "-m") mode = std::atoi(need());
        else if (a == "-t") threads = std::atoi(need());
        else if (a == "-o") outdir = need();
        else if (a == "-f") outfile = need();
        else if (a == "-h" || a == "--help") {
            printf("usage: %s -n N|a-b -k K|a-b [-m 0|1] [-t threads] [-o dir | -f file]\n", argv[0]); return 0; }
    }
    if (ns <= 0 || ks <= 0) { fprintf(stderr, "need -n and -k\n"); return 1; }
    if (threads < 1) threads = 1;

    auto t0 = std::chrono::steady_clock::now();
    std::FILE* single = nullptr;
    if (!outfile.empty()) {
        single = std::fopen(outfile.c_str(), "w");
        std::fprintf(single, "N,K,S,diameter,averageShortestPathLength,edges\n");
    } else if (outdir.empty()) {
        printf("N,K,S,diameter,averageShortestPathLength,edges\n");
    }
    for (int N = ns; N <= ne; ++N)
        for (int k = ks; k <= ke; ++k) {
            if (k > N / 2) continue;
            std::FILE* out = stdout;
            if (single) out = single;
            else if (!outdir.empty()) {
                std::filesystem::create_directories(outdir);
                std::string path = outdir + "/" + std::to_string(N) + "_" + std::to_string(k) + ".csv";
                out = std::fopen(path.c_str(), "w");
                std::fprintf(out, "N,K,S,diameter,averageShortestPathLength,edges\n");
            }
            run_one(N, k, mode, threads, out);
            if (out != stdout && out != single) std::fclose(out);
        }
    if (single) std::fclose(single);
    auto t1 = std::chrono::steady_clock::now();
    fprintf(stderr, "time: %.3f s (threads=%d, mode=%d)\n",
            std::chrono::duration<double>(t1 - t0).count(), threads, mode);
    return 0;
}
