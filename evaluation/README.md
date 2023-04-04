# Evaluation

Compilation:
```
cgal_create_CMakeLists -s evaluate #(έχει εκτελεστεί από εμάς)
cmake -DCGAL_DIR=/usr/include/CGAL -DCMAKE_BUILD_TYPE=Release
make
```
How to Run:
```
./evaluate -i <point set path> –ο <output file> -preprocess <optional>
```
e.g.
```
./evaluate -i inputs -o output.txt -preprocess
```
Clean Executables:
```
make clean
```

## Results for point sets of different size (in ms):

score = (area from polygonization+optimization algorithms) / (convex hull area)

When no feasible solution is found or cut-off time is up, score=1 is used for minimization problems and score=0 is used for maximization problems.

| Number of points | Local Search (Incremental) | --/-- | --/-- | --/-- | Simulated Annealing (Incremental) | --/-- | --/-- | --/-- | Local Search (Convex Hull) | --/-- | --/-- | --/-- | Simulated Annealing (Convex Hull) | --/-- | --/-- | --/-- | Simulated Annealing Subdivision | --/-- | --/-- | --/-- |
| :--- | :----: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Size | min score | max score | min bound | max bound | min score | max score | min bound | max bound | min score | max score | min bound | max bound | min score | max score | min bound | max bound | min score | max score | min bound | max bound |
| 10   | 1.27246   | 3.54474   | 0.363867  | 0.869049  | 1.23807   | 3.54296   | 0.351727  | 0.867276  | 2.03064   | 3.30284   | 0.890685  | 0.743282  | 1.44596   | 3.47785   | 0.455102  | 0.82156   | 1.26691   | 2.58696   | 0.406871  | 0.529454  |
| 20   | 0.852228  | 3.58056   | 0.332764  | 0.840884  | 0.952228  | 3.55323   | 0.317462  | 0.845038  | 1.61605   | 3.48791   | 1         | 0.84078   | 1.63272   | 3.57582   | 1         | 0.842028  | 1.27682   | 3.25127   | 0.363719  | 0.726507  |
| 50   | 0.709723  | 3.56806   | 0.249745  | 0.876024  | 0.624006  | 3.529     | 0.183338  | 0.847197  | 1.46733   | 2.72495   | 1         | 0         | 1.59097   | 2.72761   | 1         | 0         | 0.929413  | 3.24364   | 0.314923  | 0.76774   |
| 100  | 0.705417  | 3.36984   | 0.191648  | 0.799245  | 0.778901  | 3.33277   | 0.216446  | 0.783517  | 1.5687    | 2.66316   | 1         | 0         | 1.60884   | 2.65122   | 1         | 0         | 0.979219  | 2.91539   | 0.35736   | 0.707532  |
| 200  | 0.801393  | 3.22373   | 0.265084  | 0.743775  | 0.854492  | 3.16571   | 0.262174  | 0.75001   | 2.4296    | 2.63491   | 1         | 0         | 1.73351   | 2.63748   | 1         | 0         | 0.82771   | 2.96282   | 0.310235  | 0.686687  |
| 500  | 3.12213   | 0.782273  | 1         | 0         | 0.751382  | 3.05309   | 0.254925  | 0.732747  | 4         | 0.918169  | 1         | 0         | 1.71951   | 2.62186   | 1         | 0         | 0.785592  | 3.1461    | 0.287392  | 0.730856  |
| 1000 | 3         | 0         | 1         | 0         | 0.555045  | 2.33161   | 0.258352  | 0.710222  | 3         | 0         | 1         | 0         | 1.49372   | 1.82242   | 1         | 0         | 0.564033  | 2.36196   | 0.284507  | 0.770118  |

The above table (also in `comparison.txt`) contains results for the following point sets:

`euro-night-0000010.instance`, `euro-night-0000020.instance`, `euro-night-0000050.instance`, `euro-night-0000100.instance`, `euro-night-0000200.instance`, `euro-night-0000500.instance`, `euro-night-0001000.instance`,<br>
`stars-0000010.instance`, `stars-0000020.instance`, `stars-0000050.instance`, `stars-0000100.instance`, `stars-0000200.instance`, `stars-0000500.instance`,<br>
`uniform-0000010-1.instance`, `uniform-0000020-1.instance`, `uniform-0000050-1.instance`, `uniform-0000100-1.instance`, `uniform-0000200-1.instance`, `uniform-0000500-1.instance`, `uniform-0001000-1.instance`,<br>
`us-night-0000010.instance`, `us-night-0000020.instance`, `us-night-0000050.instance`, `us-night-0000100.instance`, `us-night-0000200.instance`, `us-night-0000500.instance`, `us-night-0001000.instance`