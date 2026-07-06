# TSP

A header-only C++ library that reads TSPLIB instances, builds a greedy initial
tour, improves it with a Lin-Kernighan-style local search (LKH), and escapes
local optima with Chained Lin-Kernighan.

## Structure

All code is header-only and lives under the `TSP` namespace. Including
`TSP/TSP.h` pulls in the rest in order.

| File | Role |
|---|---|
| [TSP.h](TSP/TSP.h) | `problem` struct definition; includes the remaining headers |
| [TSPLIBloader.h](TSP/TSPLIBloader.h) | Parse TSPLIB `.tsp` / `.opt.tour` files, extract the archive |
| [weight.h](TSP/weight.h) | Distance functions (`EUC`, `MAX`, `MAN`, `CEIL_2D`, `GEO`, `ATT`) and weight-matrix construction |
| [disjoint_set.h](TSP/disjoint_set.h) | Union-Find (used by the greedy construction) |
| [greedy.h](TSP/greedy.h) | Greedy edge-matching initial tour |
| [LKH.h](TSP/LKH.h) | Lin-Kernighan local search, Iterative/Chained LKH |

`TSP/main.cpp` is only a quick Xcode compile check, not the entry point of the
library itself.

## Algorithm

1. **Greedy construction** — greedily matches the cheapest edges first to build
   an initial tour (Christofides-style greedy edge heuristic).
2. **LKH** — repeatedly applies sequential edge exchanges (2.5-opt family, depth
   up to 5) to descend to a local optimum. `LKHbacktracking` implements this
   recursively.
3. **IterativeLKH** — the same search as `LKHbacktracking` with the recursion
   unrolled onto an explicit stack. Its result is always meant to match `LKH`,
   and the two have been confirmed to agree on several TSPLIB instances.
4. **ChainedLKH** — Chained Lin-Kernighan with the double-bridge kick of
   Martin/Otto/Felten (the approach used by Applegate/Cook/Rohe's `linkern` and
   Concorde's default LK). After reaching a local optimum, it perturbs the tour
   with a double-bridge (4-opt) kick and descends again, continuing the chain
   from the kicked tour only when the result is no worse. Reference paper:
   `tr331-04.pdf` (Fischer & Merz, *Embedding a Chained Lin-Kernighan Algorithm
   into a Distributed Algorithm*).

## Preparing the data (TSPLIB)

`ALL_tsp.tar.gz` contains [TSPLIB](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/)
instances. Unpack it as below so that `ReadFile` can read them directly.

```cpp
TSP::PATH = "/path/to/repo/root/"; // directory containing ALL_tsp.tar.gz, trailing '/' required
TSP::makeFile();                  // one-time: unpack into the ALL_tsp/ folder
```

If it is already unpacked, you can call `ReadFile` directly without `makeFile()`.

## Usage example

```cpp
#include "TSP/TSP.h"
using namespace TSP;

int main(){
    TSP::PATH = "/path/to/repo/root/";

    problem X;
    ReadFile("berlin52", X); // ALL_tsp/berlin52.tsp (+ .opt.tour fills X.optTour if present)

    std::vector<int> rout(X.dimension);
    greedy(rout.data(), X);          // initial tour
    ChainedLKH(rout.data(), X, 30);  // improve with 30 double-bridge kicks

    printf("tour length: %f\n", loss(rout, X));
}
```

Build (plain g++, without Xcode):

```sh
g++ -std=c++17 -O2 -I. your_driver.cpp -o your_driver
```

## Notes and limitations

- Distances follow the TSPLIB specification: `EUC/MAX/MAN` round, `CEIL_2D`
  rounds up, and `GEO`/`ATT` use their respective TSPLIB formulas. `GEO`
  distances have been verified to reproduce the published optimal tour lengths
  exactly (e.g. `ulysses16` = 6859, `gr202` = 40160).
- `LKH` keeps only the 5 nearest neighbors per node for its candidate list, so
  the neighbor structure is `O(n)` memory rather than `O(n^2)`. Even so, each
  `LKHstep` performs an `O(n)` scan per improvement, which makes very large
  instances (tens of thousands of cities and up) slow.
- `ChainedLKH`/`IterativeLKH` are single-threaded, single-node implementations.
  The distributed (multi-node) version discussed in the reference paper is not
  implemented.
