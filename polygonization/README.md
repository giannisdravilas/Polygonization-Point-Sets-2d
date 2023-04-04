# Polygonization

Compilation:
```
cgal_create_CMakeLists -s to_polygon #(has already been executed)
cmake -DCGAL_DIR=/usr/include/CGAL -DCMAKE_BUILD_TYPE=Release
make
```
How to Run:
```
./to_polygon -i <point set input file> –ο <output file> –algorithm <incremental or convex_hull> -edge_selection <1 or 2 or 3> -initialization <1a or 1b or 2a or 2b | only on incremental>
```
e.g.
```
./to_polygon -i ./instances/data/images/euro-night-0000100.instance -ο output.txt -algorithm convex_hull -edge_selection 1
```
Clean Executables:
```
make clean
```

- Where:   
    - `-i [InputFile]`
    - `-ο [OutputFile]`.
    - `algorithm [convex_hull ΟR incremental]`
    - `edge_selection [1 OR 2 OR 3]` for random edge selection or edge selection so that minimum or maximum area is added to the polygon, respectively.
    - `initialization [1a OR 1b OR 2a OR 2b]` sort points in the beginning using 1. x: (a) decreasing / (b) increasing or 2. y: (a) decreasing / (b) increasing
    - `onion_initialization [1 - 5]` Not used here

## Results for point sets of different size (in ms):

### Incremental Algorithm:

`Instance: euro-night / initialization 1a:`

| Number of points  | Edge Selection 1 | Edge Selection 2 | Edge Selection 3 |
| :---              |      :----:      |      :---:       |      :---:       |
| 10                | 0                | 0                | 0                |
| 50                | 0                | 0                | 0                |
| 100               | 1                | 1                | 2                |
| 200               | 6                | 6                | 8                |
| 500               | 49               | 49               | 163              |
| 1000              | 286              | 230              | 303              |
| 2000              | 930              | 1102             | 1054             |
| 5000              | 7446             | 8574             | 22386            |
| 10000             | 29105            | 30656            | 34650            |
| 20000             | 129988           | 141268           | 377706           |
| 50000             | 1104734          | 1152760          | 3704387          |

`Instance: euro-night / edge_selection 1:`

| Number of points  | Initialization 1a | Initialization 1b | Initialization 2a | Initialization 2b |
| :---              |      :----:       |      :---:        |      :---:        |      :---:        |
| 10                | 0                 | 0                 | 0                 | 0                 |
| 50                | 0                 | 0                 | 0                 | 0                 |
| 100               | 1                 | 2                 | 3                 | 2                 |
| 200               | 6                 | 8                 | 10                | 21                |
| 500               | 49                | 54                | 80                | 67                |
| 1000              | 288               | 198               | 303               | 242               |
| 2000              | 930               | 1109              | 1220              | 1177              |
| 5000              | 7477              | 9533              | 10089             | 9219              | 
| 10000             | 29381             | 37682             | 42887             | 41696             |
| 20000             | 130055            | 147331            | 165484            | 140792            |
| 50000             | 1232006           | 1276391           | 1589224           | 1189562           |

`Instance: uniform-1 / initialization 1a:`

| Number of points  | Edge Selection 1 | Edge Selection 2 | Edge Selection 3 |
| :---              |      :----:      |      :---:       |      :---:       |
| 100               | 6                | 6                | 6                |
| 1000              | 179              | 169              | 224              |
| 10000             | 18481            | 17574            | 21705            |
| 20000             | 80205            | 68677            | 95650            |

`Instance: uniform-1 / edge_selection 1:`

| Number of points  | Initialization 1a | Initialization 1b | Initialization 2a | Initialization 2b |
| :---              |      :----:       |      :---:        |      :---:        |      :---:        |
| 100               | 5                 | 6                 | 6                 | 7                 |
| 1000              | 180               | 202               | 235               | 287               |
| 10000             | 18454             | 19897             | 29598             | 31800             |
| 20000             | 80568             | 90414             | 107997            | 100126            |

### Algorithm using Convex Hull:

`Instance: euro-night`

| Number of points  | Edge Selection 1 | Edge Selection 2 | Edge Selection 3 |
| :---              |    :----:        |      :----:      |      :----:      |
| 10                | 0 ms             | 0 ms             | 0 ms             |
| 50                | 0 ms             | 8 ms             | 3 ms             |
| 100               | 6 ms             | 16 ms            | 39 ms            |
| 200               | 19 ms            | 358 ms           | 568 ms           |
| 500               | 165 ms           | 11911 ms         | 18134 ms         |
| 1000              | 1003 ms          | 202886 ms        | 290098 ms        |
| 2000              | 6909 ms          | 3214851 ms       | 3752492 ms       |
| 5000              | 93224 ms         | > 5800000 ms     | > 5800000 ms     |

`Instance: uniform-1`

| Number of points  | Edge Selection 1 | Edge Selection 2 | Edge Selection 3 |
| :---              |      :----:      |      :---:       |      :---:       |
| 10                | 0 ms             | 0 ms             | 0 ms             |
| 50                | 0 ms             | 1 ms             | 3 ms             |
| 100               | 2 ms             | 20 ms            | 44 ms            |
| 500               | 186 ms           | 12899 ms         | 18806 ms         |