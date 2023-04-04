# Area Optimization

Compilation:
```
cgal_create_CMakeLists -s optimal_polygon #(έχει εκτελεστεί από εμάς)
cmake -DCGAL_DIR=/usr/include/CGAL -DCMAKE_BUILD_TYPE=Release
make
```
How to Run:
```
./optimal_polygon -i <point set input file> –o <output file> –algorithm <local_search or simulated_annealing> -initialization_algorithm <incremental or convex_hull> -L <int> <-max or -min> -threshold <double | only in local_search> –annealing <"local" or "global" or "subdivision" | only in simulated annealing> 
```
e.g.
```
./optimal_polygon -i instances/data/images/euro-night-0000100.instance -o output.txt -algorithm local_search -initialization_algorithm incremental -L 5 -max -threshold 100000

./optimal_polygon -i instances/data/images/euro-night-0000100.instance -o output.txt -algorithm local_search -initialization_algorithm convex_hull -L 5 -min -threshold 100000

./optimal_polygon -i instances/data/images/euro-night-0000100.instance -o output.txt -algorithm simulated_annealing -initialization_algorithm convex_hull -L 5000 -min -annealing local

./optimal_polygon -i instances/data/images/euro-night-0000100.instance -o output.txt -algorithm simulated_annealing -initialization_algorithm convex_hull -L 5000 -max -annealing subdivision
```
Clean Executables:
```
make clean
```

- Where:   
    - `-i [InputFile]`
    - `-ο [OutputFile]`
    - `algorithm [local_search ΟR simulated_annealing]`
    - `initialization_algorithm [incremental OR convex_hull]`
    - `L [int]` Local Search: Maximum number of consecutive points that can be used to replace an edge / Simulated Annealing: Maximum number of iterations
    - `max or min` Maximize or Minimize area
    - `threshold [double]` area changes below this threshold are considered negligible and the algorithm stops
    - `annealing [local OR global OR subdivision]`

## Results for point sets of different size (in ms):

### Local Search:

`Instance: euro-night / max, threshold = 100000:`

| Number of points  | L = 5  | L = 10  | Area Improvement (best, if applicable) |
| :---              | :----: |  :---:  |         :---:                          |
| 10                | 1      | 1       | 3455020                                |
| 100               | 2734   | 5480    | 8128348                                |
| 200               | 46923  | 95608   | 8476742                                |
| 300               | 112328 | 229700  | 9870532                                |
| 500               | 725600 | 1486033 | 6244068 (L=10)                         |

`Instance: euro-night / min, threshold = 100000:`

| Number of points  | L = 5  | L = 10 | Area Improvement (best) |
| :---              | :----: |  :---: |         :---:           |
| 10                | 2      | 3      | -111900                 |
| 100               | 2352   | 5062   | -3902290                |
| 200               | 42414  | 90914  | -11022264 (L=10)        |
| 300               | 118977 | 253350 | -7786334                |
| 500               | 350494 | 736229 | -2177374                |

`Instance: euro-night / max, threshold = 10000:`

| Number of points  | L = 5  | L = 10  | Area Improvement (best) |
| :---              | :----: |  :---:  |         :---:           |
| 10                | 1      | 1       | 3455020                 |
| 100               | 5365   | 10549   | 9369808                 |
| 200               | 106463 | 215080  | 10306934                |
| 300               | 362285 | 736075  | 11511030                |

`Instance: euro-night / min, threshold = 10000:`

| Number of points  | L = 5  | L = 10 | Area Improvement (best) |
| :---              | :----: |  :---: |         :---:           |
| 10                | 2      | 3      | -111900                 |
| 100               | 4239   | 9087   | -4512850                |
| 200               | 128093 | 256280 | -13683222 (L=10)        |
| 300               | 118977 | 253350 | -7786334                |

`Instance: uniform-1 / max, threshold = 1000:`

| Number of points  | L = 5  | L = 10  | Area Improvement (best) |
| :---              | :----: |  :---:  |         :---:           |
| 10                | 1      | 2       | 16918                   |
| 100               | 6646   | 12824   | 4472248                 |
| 200               | 106758 | 213945  | 16159108                |

`Instance: uniform-1 / min, threshold = 1000:`

| Number of points  | L = 5  | L = 10 | Area Improvement (best) |
| :---              | :----: |  :---: |         :---:           |
| 10                | 2      | 3      | -12908                  |
| 100               | 9021   | 20457  | -4361692                |
| 200               | 97073  | 212873 | -15432382               |

### Simulated Annealing:

(Total Time: Time only for `simulated_annealing_algorithm` + Time only for `convex_hull_algorithm`)

`Instance: euro-night / max, local Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time  |    Area Improvement     |Total Time: |
| :---              | :----: |       :---:             |:---:       |
| 10                | 8 ms   | 0                       | 8 ms       |
| 100               | 40 ms  | 67762                   | 88 ms      |
| 200               | 74 ms  | 207854                  | 612 ms     |
| 300               | 116 ms | 90254                   | 2680 ms    |
| 500               | 231 ms | 137746                  | 17455 ms   |

`Instance: euro-night / min, local Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time  |    Area Improvement     | Total Time: |
| :---              | :----: |       :---:             |:---:        |
| 10                | 6 ms   | 0                       | 6 ms        |
| 100               | 38 ms  | -1997462                | 55 ms       |
| 200               | 77 ms  | -2196300                | 448 ms      |
| 300               | 116 ms | -3929434                | 1625 ms     |
| 500               | 244 ms | -2038794                | 13230 ms    |

`Instance: euro-night / max, global Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time |    Area Improvement     |Total Time:     |
| :---              | :----:|       :---:             |:---:           |
| 10                | 8ms   | 936122                  | 8 ms           |
| 100               | 26ms  | 1402862                 | 62 ms          |
| 200               | 32ms  | 520892                  | 573 ms         |
| 300               | 51ms  | 91406                   | 2562 ms        |
| 500               | 73ms  | 73414                   | 18258 ms       |

`Instance: euro-night / min, global Step (init_algorithm: convexHull):, -L=5000`

| Number of points  |  Time |    Area Improvement     | Total Time: |
| :---              | :----:|       :---:             |     :---:   |
| 10                | 8 ms  | 0                       | 8 ms        |
| 100               | 20 ms | -6349356                | 36 ms       |
| 200               | 32 ms | -4058084                | 403 ms      |
| 300               | 48 ms | -3461766                | 1525 ms     |
| 500               | 62 ms | -790746                 | 12463 ms    |

Οι παρακάτω χρόνοι αποτελούν τους συνολικούς χρόνους μόνο για τον `simulated_annealing_algorithm_subdivision`.

`Instance: euro-night / max, subdivision(init_algorithm: convexHull), -L=5000:`

| Number of points   | Time   |    Area Improvement     |
| :---               | :----: |         :---:           |
| 100                | 238 ms | 8161480                 |
| 500                | 556 ms | 2284198                 |
| 1000               | 1620 ms| 2588488                 |
| 2000               | 7854 ms| 3682464                 |

`Instance: euro-night / min, subdivision(init_algorithm: convexHull), -L=5000:`

| Number of points   | Time   | Area Improvement |
| :---               | :----: |         :---:    |
| 100                | 269 ms | -931868          |
| 500                | 621 ms | -4727694         |
| 1000               | 1443 ms| -7216280         |
| 2000               | 6015 ms| -28694320        |

For `Instance: uniform-1`:

`Instance: uniform-1 / max, local Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time  |    Area Improvement     |Total Time: |
| :---              | :----: |       :---:             |:---:       |
| 10                | 6 ms   | 0                       | 8 ms       |
| 100               | 36 ms  | 43602                   | 78 ms      |
| 200               | 70 ms  | 800426                  | 603 ms     |
| 300               | 108 ms | 714186                  | 2626 ms    |
| 500               | 219 ms | 2878424                 | 18535 ms   |

`Instance: uniform-1 / min, local Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time  |    Area Improvement     | Total Time: |
| :---              | :----: |       :---:             |:---:        |
| 10                | 5 ms   | 0                       | 6 ms        |
| 100               | 37 ms  | -474320                 | 56 ms       |
| 200               | 68 ms  | -3840782                | 419 ms      |
| 300               | 109 ms | -7586304                | 1774 ms     |
| 500               | 240 ms | -19696462               | 13209 ms    |

`Instance: uniform-1 / max, global Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time  |    Area Improvement     |Total Time:     |
| :---              | :----: |       :---:             |:---:           |
| 10                | 3 ms   | 0                       | 3 ms           |
| 100               | 20 ms  | 866224                  | 63 ms          |
| 200               | 33 ms  | 1333252                 | 560 ms         |
| 300               | 48 ms  | 1037444                 | 2538 ms        |
| 500               | 67 ms  | 1901768                 | 17581 ms       |

`Instance: uniform-1 / min, global Step (init_algorithm: convexHull), -L=5000:`

| Number of points  |  Time |    Area Improvement     | Total Time: |
| :---              | :----:|       :---:             |     :---:   |
| 10                | 4 ms  | -20748                  | 4 ms        |
| 100               | 25 ms | -2557984                | 47 ms       |
| 200               | 30 ms | -8629698                | 410 ms      |
| 300               | 50 ms | -6495518                | 2056 ms     |
| 500               | 67 ms | -6587680                | 12237 ms    |

Total time only for `simulated_annealing_algorithm_subdivision`:

`Instance: uniform-1 / max, subdivision(init_algorithm: convexHull), -L=5000:`

| Number of points   | Time   |    Area Improvement     |
| :---               | :----: |         :---:           |
| 100                | 280 ms | 3385188                 |
| 500                | 566 ms | 35184418                |
| 1000               | 1448 ms| 240249750               |
| 2000               | 459  ms| 1873695818              |

`Instance: uniform-1 / min, subdivision (init_algorithm: convexHull), -L=5000:`

| Number of points   | Time   | Area Improvement |
| :---               | :----: |         :---:    |
| 100                | 198 ms | -681708          |
| 500                | 575 ms | -86101538        |
| 1000               | 1425 ms| -519210578       |
| 2000               | 465 ms | -4235069196      |