# optimal-circulants-s1

Multithreaded exact search for optimal circulant graphs
**C(N; 1, s₂, …, s_k)** with the first generator fixed to 1.

For a given number of vertices `N` and degree parameter `k`, the program enumerates
every signature `(s₂, …, s_k)` and reports the graphs that are optimal under the
criterion: **minimize the average shortest path length first, then the diameter**.
All optimal graphs are output, including every isomorphic copy, in lexicographic order.

## Build

```
g++ -O3 -march=native -funroll-loops -std=c++17 -pthread circulant_mt_fast.cpp -o circulant_mt
```

Optional profile-guided build (typically 10–15% faster on a fixed workload):

```
g++ -O3 -march=native -fprofile-generate -std=c++17 -pthread circulant_mt_fast.cpp -o cm_pgo
./cm_pgo -n 300 -k 5 -t 16 -o prof
g++ -O3 -march=native -fprofile-use -std=c++17 -pthread circulant_mt_fast.cpp -o circulant_mt
```

The code is standard C++17 and builds with g++ or clang on Linux, and with MSYS2
UCRT64 g++ on Windows.

## Usage

```
circulant_mt -n N|a-b -k K|a-b [-m 0|1] [-t threads] [-o dir | -f file]
```

| Flag | Meaning |
|------|---------|
| `-n` | number of vertices, single value or range `a-b` |
| `-k` | degree parameter (half-degree), single value or range |
| `-m` | lower-bound mode: `1` tight sphere-capacity bound (default), `0` simple bound |
| `-t` | thread count (default: hardware concurrency) |
| `-o` | output directory, one file `{N}_{k}.csv` per pair |
| `-f` | single output file for the whole range |

With neither `-o` nor `-f`, results are printed to stdout.

Output columns: `N, K, S, diameter, averageShortestPathLength, edges`,
where `S` is the signature written as `C(N;1;s2;...;sk)`.

Example:

```
circulant_mt -n 100-600 -k 5 -t 16 -f result.csv
```

## Method

The search is an exact exhaustive enumeration over the signature space
`C(⌊N/2⌋−1, k−1)`, with the following techniques:

- **Isomorphism reduction.** Equivalent signatures (multiplication of all
  generators by an invertible element of Z_N) are filtered by a canonicality test
  that performs `k−1` normalizations rather than `2k`, using the identity
  `fold(−x) = fold(x)` and the trivial unit generator.
- **Barrett reduction.** The modular multiplication in the canonicality test uses
  a precomputed constant `μ = ⌊2⁶⁴/N⌋` instead of integer division.
- **Pruned BFS.** A lower bound on the distance sum, based on L1-sphere capacities,
  aborts the breadth-first search at the first layer where the candidate cannot beat
  the current record.
- **Central symmetry.** The relation `d(v) = d(N−v)` lets one probe mark a pair of
  vertices, halving the traversal.
- **Generation-tagged distances.** The distance array is not cleared between
  candidates; a per-vertex epoch tag marks validity, removing the O(N) reset.
- **Thread balancing.** Worker threads pull linearized `(s₂, s₃)` pairs from an
  atomic counter for fine-grained load balance; the shared record is cache-line
  aligned to avoid false sharing.

The asymptotic order of the exact search is unchanged by these techniques; they
reduce the constant factor and the cost of deduplication.

## Performance

Benchmarked against `pcgpp` (Moskalenko) on an AMD Ryzen 9 5980HX (8 cores / 16
threads), 16 threads, exact mode (eps = 0), single fixed generator, output to file,
median of 7 runs. Outputs of both programs match up to isomorphism on the whole grid.

Relative speedup grows with the degree parameter `k`, because the share of work spent
on the canonicality test grows with the number of generators, which is where the
optimizations concentrate:

| k | speedup vs pcgpp |
|---|------------------|
| 2 | 0.74 – 0.78 (slower; BFS-dominated, overhead not amortized) |
| 3 | 0.75 – 1.19 (near parity) |
| 4 | up to ≈ 2.8 (N = 400–750) |
| 5 | 1.8 – 2.7 (N = 200–750) |

Compared with an older single-threaded baseline, the speedup reaches one to two
orders of magnitude at large `N`, reflecting the absence of pruning and symmetry in
that baseline rather than a property of the present work.

Scaling from 8 to 16 threads is about 1.5×, consistent with 8 physical cores plus
simultaneous multithreading. At full thread occupancy the bottleneck is memory
bandwidth rather than per-thread compute, so the exact search is near its
architectural ceiling on this class of hardware.

## License

See [LICENSE](LICENSE).
