# Parallel Gauss-Jordan KJI in MPI C

An implementation of  parallel Gauss-Jordan method in KJI form written in MPI C.

You could give your augmented matrix in a txt as an argument. If you don't then a random augmented matrix is generated. The "matrix_example.txt" in src folder is a txt file example.

We implement two variants of Gauss-Jordan method. The first seperates columns by zones and the second by card-shuffling.

In figures folder we give executions of both variation for p=2,4,8,16. As you can notice the executional times, the speedup and the efficiency for both variations follows the theoritical results.

We use MPI 3.1 and gcc 6.

In src folder run:

~$ mpicc suffle.c -o executable   (or ~$ mpicc zones_col.c -o executable)

~$ mpirun -np number_of_processes ./executable aug_matrix.txt

These implementations was developed for the course of Parallel Algorithms in University.
