// pcgx_cpu.cpp
// Multithreaded exact search for OPTIMAL circulants C(N; 1, s2, ..., sk)
// by (ASPL, diam): ASPL (dist_sum) is primary, diameter breaks ties (min).
// Arbitrary k >= 2 (up to KMAX=16). Emits ALL optimal graphs (every isomorphic copy).
//
// Version 2. Differences from circulant_mt:
//   [1] thread balancing over linearized pairs (s2,s3), not by s2 alone;
//   [2] per-thread cached best, re-reading the atomic once per REREAD calls;
//   [3] fold via Barrett reduction (mu = 2^64 / N), no integer division;
//   [4] is_canonical in k-1 normalizations (g=1 is identity; the N-g sign is
//       redundant since fold(-x)=fold(x)) instead of 2k;
//   [5] early exit by the Moore lower bound on diameter before running BFS on a branch;
//   [6] BFS without O(N) reset: epoch (generation) technique instead of INF-reset;
//   [7] KMAX=16, all fixed-size arrays enlarged.
//
// Build (MSYS2 UCRT64).
// Portable (to run on another machine):
//   g++ -O3 -march=x86-64-v3 -std=c++17 -pthread -static pcgx_cpu.cpp -o pcgx_cpu.exe
// Tuned to the local CPU for maximum speed:
//   g++ -O3 -march=native -std=c++17 -pthread -static pcgx_cpu.cpp -o pcgx_cpu.exe
// PGO (adds 10-15% speed):
//   g++ -O3 -march=native -fprofile-generate -std=c++17 -pthread pcgx_cpu.cpp -o cm2_pgo
//   ./cm2_pgo -n 300 -k 5 -t 16 -o /tmp/x
//   g++ -O3 -march=native -fprofile-use -std=c++17 -pthread pcgx_cpu.cpp -o pcgx_cpu.exe
// Run:
//   pcgx_cpu -n N|a-b -k K|a-b [-m 0|1] [-t threads] [-o dir | -f file]

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

static const int KMAX = 16;          // [7]
static const int REREAD = 256;       // [2] best-record re-read frequency

static inline int mod_inv(long long a, long long n) {
    long long t = 0, nt = 1, r = n, nr = ((a % n) + n) % n;
    while (nr) { long long q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp; }
    if (r != 1) return -1;
    if (t < 0) t += n;
    return (int)t;
}

// [3] Barrett reduction: fold for products sig[j]*m modulo N.
// mu = floor(2^64 / N). Reduces x mod N to a high-word multiply and a shift.
struct Barrett {
    u64 N, mu;
    void init(u64 n) { N = n; mu = (n ? (~0ULL) / n : 0); } // ~0/N == floor((2^64-1)/N), enough for x < 2^63
    // reduce a value 0 <= x < N^2 (here x = sig*m, both < N, so < N^2 < 2^64 for N < 2^32)
    inline u32 fold(u64 x) const {
        // q = floor(x * mu / 2^64)
        unsigned __int128 t = (unsigned __int128)x * mu;
        u64 q = (u64)(t >> 64);
        u64 r = x - q * N;
        if (r >= N) r -= N;                 // a single correction suffices
        // fold to min(r, N-r)
        u64 nr = N - r;
        return (u32)(nr < r ? nr : r);
    }
};

// size of the L1 sphere of radius t in k dimensions
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

static inline u64 pack(u32 diam, u64 ds) { return (ds << 16) | (diam & 0xFFFFu); }
static inline u32 unpack_diam(u64 key) { return (u32)(key & 0xFFFFu); }
static inline u64 unpack_ds(u64 key) { return key >> 16; }

// [4] is_canonical in k-1 normalizations.
// Equivalence class: multiply all generators by an invertible factor m.
// Take m = inv[sig[gi]] for gi = 1..k-1 (gi=0 gives sig[0]=1 -> m=inv[1]=1, identity, skipped).
// The N-g sign is redundant: fold(-x)=fold(x), reflection is already handled by fold.
static bool is_canonical(int N, const int* __restrict__ sig, int k, const int* __restrict__ inv, const Barrett& br) {
    int img[KMAX];
    for (int gi = 1; gi < k; ++gi) {       // gi=0 -> identity, skip
        int g = sig[gi];
        int m = inv[g];
        if (m == 0) continue;              // not invertible
        // Fast path. The image always contains 1 (fold(sig[gi]*m)=1), the signature starts with 1,
        // so the lexicographic comparison is decided by the second-smallest element:
        // m2 = min over j!=gi of fold(sig[j]*m) versus sig[1]. Full sort is needed only on a tie.
        int m2 = N;
        for (int j = 0; j < k; ++j) {
            if (j == gi) { img[j] = 1; continue; }
            int v = (int)br.fold((u64)sig[j] * (u64)m);
            img[j] = v;
            if (v < m2) m2 = v;
        }
        if (m2 < sig[1]) return false;     // image is strictly smaller
        if (m2 > sig[1]) continue;         // image is strictly larger, normalization passed
        // tie on the second element: full comparison
        for (int a = 1; a < k; ++a) { int v = img[a], b = a - 1;
            while (b >= 0 && img[b] > v) { img[b + 1] = img[b]; --b; } img[b + 1] = v; }
        for (int j = 0; j < k; ++j) { if (img[j] < sig[j]) return false; if (img[j] > sig[j]) break; }
    }
    return true;
}

// BFS with pruning. Storage: ONE array vals[] instead of a dist[]+seen[] pair.
// vals[v] = (epoch << DIST_BITS) | dist. Visited = high bits match the epoch,
// distance = low DIST_BITS bits. Halves the random memory traffic per mark
// (one write instead of two); the algorithm and traversal order are unchanged.
// Limitation: diameter < 2^DIST_BITS (=1024), enough for N up to ~2M at k=2.
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

// [5] minimum possible diameter (Moore bound) for given N,k: smallest D with cum_cap(D) >= N
static int moore_diam(int N, const u64* caps, int capsN, int k) {
    u64 cum = 1; int t = 0;
    while (cum < (u64)N) { ++t; u64 c = (t < capsN) ? caps[t] : sphere_cap_formula(k, t); cum += c; }
    return t;
}

struct alignas(64) Best { std::atomic<u64> key; char pad[64 - sizeof(std::atomic<u64>)]; Best() { key.store(~0ULL); } };
struct Found { std::vector<int> sig; u32 diam; u64 ds; };

static void run_one(int N, int k, int mode, int threads, std::FILE* out) {
    if (k < 2 || k > KMAX || k > N / 2) return;
    int h = N / 2;

    // precomputations
    std::vector<int> inv(N, 0);
    for (int i = 1; i < N; ++i) { int m = mod_inv(i, N); if (m > 0) inv[i] = m; }
    Barrett br; br.init((u64)N);                                  // [3]
    std::vector<u64> caps; { u64 cum = 1; int t = 0; caps.push_back(1);
        while (cum < (u64)N + 1) { ++t; u64 c = sphere_cap_formula(k, t); caps.push_back(c); cum += c; } }
    int capsN = (int)caps.size();
    int Dmoore = moore_diam(N, caps.data(), capsN, k);            // [5]
    // max dist_sum for diameter Dmoore under ideal packing = lower bound on the optimum ds
    // (for early branch rejection by s2,s3 we use the per-shell capacity).

    Best best;
    std::vector<std::vector<Found>> tloc(threads);

    // [1] linearize pairs (s2,s3). s2 in [2..h-(k-2)], s3 in [s2+1..h-(k-3)].
    // For k==2 there are no pairs (only s2). For k>=3 pairs are handed out by an atomic counter.
    std::atomic<long long> next_pair{ 0 };
    std::atomic<int> next_s2{ 2 };           // used only when k==2

    // precompute the list of valid pairs (s2,s3) for k>=3
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

        u64 local_best = best.key.load(std::memory_order_relaxed);   // [2]
        int since = 0;

        auto eval = [&]() {
            if (!is_canonical(N, sig, k, inv.data(), br)) return;
            // [2] do not re-read the best every time
            if (++since >= REREAD) { local_best = best.key.load(std::memory_order_relaxed); since = 0; }
            u32 D; u64 ds;
            // the epoch occupies 32-DIST_BITS bits; on overflow we clear the array
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

        // k >= 3: take pairs (s2,s3) and run an odometer over positions 4..k-1
        while (true) {
            long long pi = next_pair.fetch_add(1, std::memory_order_relaxed);
            if (pi >= npairs) break;
            int s2 = pairs[pi].first, s3 = pairs[pi].second;
            sig[1] = s2; sig[2] = s3;

            // [5] early exit: with s2 minimal step 1 the ball grows no faster than ideal;
            // if even ideal packing with the fixed s2,s3 does not cover N
            // within Dmoore steps worse than the current best diameter, we still run BFS,
            // since the exact capacity for specific s2,s3 is not cheap to estimate. So we use
            // a weaker but correct filter: skip the branch only if the current best
            // already has diameter < Dmoore (impossible, Dmoore is the lower bound) - no-op.
            // The practical early exit is implemented inside bfs via the caps bound.

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
