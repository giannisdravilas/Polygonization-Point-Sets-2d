# Polygonization-Point-Sets-2d

📐 Implementation of Polygonization and Area Optimization algorithms for point sets in the plane. An app that evaluates the performance of these algorithms in terms of accuracy and speed is also included.

In collaboration with [Panagiotis Drivas](https://github.com/PanagiotisDrivas)

## Polygonization

`polygonization` directory contains the implementation of two polygonization algorithms for point sets containing 2d points.

- Incremental Algorithm
- Algorithm using Convex Hull

## Area Optimization

`optimization` directory contains the implementation of two area optimization algorithms for the polygons created using the polygonization algorithms.

- Local Search Algorithm
- Simulated Annealing Algorithm (with Local and Global Transition Steps)

## Evaluation

`evaluation` directory contains the implementation of an app that compares the accuracy and speed for each of the 9 polygonization-optimization combinations formed using the previously mentioned algorithms.

## Instances

`instances` directory contains example instances used for the development and evaluation of these programs.

<br>

`CGAL 5.5` along with `CMake` are used for this project.