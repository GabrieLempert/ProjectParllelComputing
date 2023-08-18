## Parallel Point Validation using MPI - README

This repository contains a parallelized program written in C that performs validation on a set of points using the Message Passing Interface (MPI). The program calculates coordinates for each point based on different time steps and checks if proximity criteria are met between points. The validated points are then recorded and saved to an output file.

### Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [How it Works](#how-it-works)
- [Parallelization Strategy](#parallelization-strategy)
- [Compiling and Running](#compiling-and-running)
- [Input File Format](#input-file-format)
- [Output](#output)
- [Sample Input and Output](#sample-input-and-output)
- [Authors](#authors)
- [License](#license)

### Overview

This program is designed to validate points in parallel by distributing the validation process across multiple MPI processes. It leverages MPI's message-passing capabilities to efficiently distribute and gather data among processes.

### Requirements

To compile and run the program, you will need:

- A C compiler (e.g., GCC)
- MPI library (e.g., OpenMPI or MPICH)

### Usage

1. Clone this repository to your local machine.

2. Compile the program using your C compiler and MPI flags. For example:

   ```sh
   make
   ```

3. Prepare your input data in a text file named `input.txt`. See [Input File Format](#input-file-format) for details.

4. Run the compiled program using MPI. For example:

   ```sh
   mpiexec -n <num_processes> ./build/program
   ```

   Replace `<num_processes>` with the number of MPI processes you want to use.

5. Check the generated `output.txt` file for the validation results.

### How it Works

The program works by performing the following steps:

1. **Initialization**: The program reads input parameters from the `input.txt` file, including the number of points, criteria parameters, and the number of time steps.

2. **Data Distribution**: In the root process (rank 0), the data is read from the file and then distributed to other processes using MPI. The time steps are divided equally among processes, and each process receives its share of the points.

3. **Parallel Validation**: Each process calculates the coordinates of points for its assigned time steps using the `CalculateCoordinate` function. It then validates the points using the `CheckCretirieaByT` function.

4. **Local Results**: Each process maintains a local matrix that stores IDs of points satisfying the criteria for its assigned time steps.

5. **Data Gathering**: The root process gathers the local matrices from all processes and combines them into a global matrix.

6. **Final Validation and Output**: The root process further validates the global matrix, and the results are written to the `output.txt` file.

### Parallelization Strategy

The program uses MPI to parallelize the validation process. It employs a master-worker model, where the root process (rank 0) reads the data, distributes it to worker processes, collects the results, and generates the final output.

### Compiling and Running

- Compile the program using your C compiler and MPI flags, as shown in the [Usage](#usage) section.
- Run the program using MPI, specifying the desired number of MPI processes.

### Input File Format

The input data is read from a text file named `input.txt`. The format of the file is as follows:

```
N K D TCount
id1 x1 x2 a b
id2 x1 x2 a b
...
idN x1 x2 a b
```

- `N`: Number of points.
- `K`: Number of points required to satisfy criteria.
- `D`: Maximum distance threshold.
- `TCount`: Total number of time steps.
- For each point, provide `id`, `x1`, `x2`, `a`, and `b`.

### Output

The validation results are written to the `output.txt` file. The output includes information about points that satisfy the proximity criteria at each time step. If no points satisfy the criteria for a time step, an appropriate message is displayed.

### Sample Input and Output

- Sample `input.txt`:

  ```
  4 2 1.5 3
  1 0 5 2 1
  2 1 6 1 0
  3 2 7 3 2
  4 3 8 0 -1
  ```

- Sample `output.txt`:

  ```
  Points 1, 2, 3 satisfy Proximity Criteria at t = t1
  Points 1, 3, 2 satisfy Proximity Criteria at t = t2
  Points 1, 3, 2 satisfy Proximity Criteria at t = t3
  No points satisfy the proximity criteria.
  ```

### Authors

- Gabriel Lempert

### License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.