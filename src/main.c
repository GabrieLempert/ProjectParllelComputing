// Created by: Gabriel Lempert

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#define MAX_POINTS 100000
#define FILE_NAME "input.txt"
#define OUTPUT_NAME "output.txt"
/*Data*/
typedef struct
{
    int id;
    double x1, x2, a, b;
} Point;

typedef struct
{
    int N, K, TCount;
    float D;
    Point *points;
} init_p;

typedef struct
{
    int id;
    double x, y;
} Coordinate;
/*Functions*/
void initData(init_p *data, const char *filename);
int checkDistanceSmallerThanD(double x1, double y1, double x2, double y2, float D);
void CalculateCoordinate(Point *p, int t, Coordinate *coord, int N, int index);
int CheckCretirieaByT(init_p *data, Coordinate *coord, int index);

/**
 * initData - Initialize the data structure and read input from a file.
 *
 * This function initializes the fields of the init_p structure with default values,
 * reads input parameters from the specified input file, and allocates memory for the
 * points array. It then reads the point data from the file and assigns it to the
 * corresponding fields in the structure. If any errors occur during this process,
 * appropriate error messages are displayed, and the program exits with an error code.
 *
 * @param data        Pointer to the init_p structure to be initialized.
 * @param filename    Name of the input file containing the data.
 */
void initData(init_p *data, const char *filename)
{
    // Initialize the struct members
    data->N = 0;
    data->K = 0;
    data->TCount = 0;
    data->D = 0.0;
    data->points = NULL;

    // Read the data from the file and assign it to the struct members
    FILE *inputFile = fopen(filename, "r");
    if (inputFile == NULL)
    {
        printf("Failed to open the input file.\n");
        exit(1);
    }

    int result = fscanf(inputFile, "%d %d %f %d", &(data->N), &(data->K), &(data->D), &(data->TCount));
    if (result != 4)
    {
        printf("Error reading input parameters from the file.\n");
        fclose(inputFile);
        exit(1);
    }

    if (data->N <= 0 || data->N > MAX_POINTS)
    {
        printf("Invalid value for N in the input file.\n");
        fclose(inputFile);
        exit(1);
    }

    // Allocate memory for the points array based on N
    data->points = (Point *)malloc(data->N * sizeof(Point));
    if (data->points == NULL)
    {
        printf("Memory allocation error for points array.\n");
        fclose(inputFile);
        exit(1);
    }

    // Read point data from the file and assign it to the points array
    for (int i = 0; i < data->N; i++)
    {
        result = fscanf(inputFile, "%d %lf %lf %lf %lf", &(data->points[i].id), &(data->points[i].x1), &(data->points[i].x2), &(data->points[i].a), &(data->points[i].b));
        if (result != 5)
        {
            printf("Error reading data for Point %d from the file.\n", i + 1);
            fclose(inputFile);
            exit(1);
        }
    }

    fclose(inputFile);
}


int checkDistanceSmallerThanD(double x1, double y1, double x2, double y2, float D)
{
    double distance = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    return (distance < D) ? 1 : 0;
}

/**
 * CalculateCoordinate - Calculate the coordinates of points based on time step.
 *
 * This function calculates the coordinates of points for a given time step 't'
 * using the given parameters and point data. The calculation is parallelized
 * using OpenMP to improve performance.
 *
 * @param p      Pointer to the array of points.
 * @param t      The time step value for which coordinates are calculated.
 * @param coord  Pointer to the array of Coordinate structures to store results.
 * @param N      The number of points in the array.
 * @param index  The index for identifying the range of time steps.
 */
void CalculateCoordinate(Point *p, int t, Coordinate *coord, int N, int index)
{   
    // Parallelize the loop using OpenMP
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        // Calculate x and y coordinates using given formulas
        double x = ((p[i].x2 - p[i].x1) / 2.0) * sin(t * M_PI / 2.0) + (p[i].x2 + p[i].x1) / 2.0;
        double y = p[i].a * x + p[i].b;

        // Assign calculated values to the Coordinate structure
        coord[i].id = p[i].id + N * index;
        coord[i].x = x;
        coord[i].y = y;
    }
}


/**
 * CheckCretirieaByT - Check if the proximity criteria are met for a specific time index.
 *
 * This function checks if the proximity criteria are met for a given time index 'index'.
 * It compares the coordinates of the point at the specified index with the coordinates
 * of other points and calculates the number of points that satisfy the distance condition.
 * If the count of satisfying points is equal to or greater than 'K_approximates', it indicates
 * that the criteria are met.
 *
 * @param data        Pointer to the init_p structure containing input data.
 * @param coord       Array of coordinates calculated for the specific time index.
 * @param index       The time index for which the criteria are being checked.
 * @return            Returns 1 if criteria are met, otherwise returns 0.
 */
int CheckCretirieaByT(init_p *data, Coordinate *coord, int index)
{
    // Create a temporary coordinate for the point at the specified index
    Coordinate temp;
    temp.x = coord[index].x;
    temp.y = coord[index].y;
    temp.id = coord[index].id;
    
    int K_approximates = data->K; // Number of required satisfying points
    float D_maximum = data->D;     // Maximum distance threshold
    int count = 0;                 // Counter for satisfying points

    // Loop through all points to check if the criteria are met
    for (int i = 0; i < data->N; i++)
    {
        if (index != i) // Avoid checking the point against itself
        {
            // Check if the distance between the current point and the temp point is smaller than D_maximum
            count += checkDistanceSmallerThanD(temp.x, temp.y, coord[i].x, coord[i].y, D_maximum);
        }

        if (count >= K_approximates)
        {
            // Criteria are met if the count of satisfying points is equal to or greater than K_approximates
            return 1;
        }
    }

    return 0; // Criteria are not met
}

/**
 * validate - Perform validation based on specified criteria and populate the idMatrix.
 *
 * This function validates the points based on specified criteria using the provided data and
 * coordinate array. It calculates the coordinates for each point using the CalculateCoordinate
 * function, then checks the criteria using the CheckCretirieaByT function. If the criteria
 * are satisfied, the point's ID is added to the idMatrix for the corresponding time step and
 * validation count. The process is repeated for each time step in the range determined by
 * startTCount and size.
 *
 * @param data           Pointer to the init_p structure containing input data.
 * @param idMatrix       2D array to store IDs of points satisfying criteria for each time step.
 * @param startTCount    Starting time step for validation range.
 * @param size           Number of time steps to validate.
 * @param rank           Rank of the MPI process.
 */
void validate(init_p *data, int idMatrix[][3], int startTCount, int size, int rank)
{
    Coordinate *coord = (Coordinate *)malloc(sizeof(Coordinate) * data->N);
    // Check if memory allocation for coord was successful
    if (coord == NULL)
    {
        printf("Memory allocation error for coord array.\n");
        exit(1);
    }
    int validationCount = 0;

    // For each cell in the idMatrix, initialize with -1 as the default value
    for (int i = 0; i < size; i++)
    {
        int calculateFrom = startTCount + i;
        int t = 2.0 * calculateFrom / data->TCount - 1;
        CalculateCoordinate(data->points, t, coord, data->N, calculateFrom);

        for (int j = 0; j < data->N; j++)
        {
            if (validationCount == 3)
            {
                break;
            }

            if (CheckCretirieaByT(data, coord, j))
            {
                idMatrix[i][validationCount] = coord[j].id;
                validationCount++;
            }
        }
        validationCount = 0;
    }
}

/**
 * WriteToOutputFile - Write validation results to an output file.
 *
 * This function writes the validation results to an output file with the specified filename.
 * It iterates through the idMatrix and writes information about points that satisfy the
 * proximity criteria at each time step. If no points satisfy the criteria, it writes a
 * corresponding message. The count of points that satisfy the criteria is also calculated.
 *
 * @param filename   Name of the output file to write the results to.
 * @param idMatrix   2D array containing IDs of points that satisfy criteria for each time step.
 * @param TCount     Total number of time steps.
 */
void WriteToOutputFile(const char *filename, int idMatrix[][3], int TCount)
{
    FILE *outputFile = fopen(filename, "w");
    if (outputFile == NULL)
    {
        printf("Error opening output file: %s\n", filename);
        return;
    }
    int count = 0;

    for (int i = 0; i < TCount; i++)
    {
        int id1 = idMatrix[i][0];
        int id2 = idMatrix[i][1];
        int id3 = idMatrix[i][2];

        if (id3 > 0)
        {
            count++;
            fprintf(outputFile, "Points %d, %d, %d satisfy Proximity Criteria at t = t%d\n", id1, id2, id3, i + 1);
        }
    }
    if (count == 0)
    {
        fprintf(outputFile, "No points satisfy the proximity criteria.\n");
    }

    fclose(outputFile);
}

int main(int argc, char *argv[])
{
    init_p *data =(init_p*)malloc(sizeof(init_p));
    int rank, size;
    int *tCountPerProccess;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
    {

        initData(data, FILE_NAME);
        // calculate each proccess start and end t count and save each of them in an array and send it to each proccess
        tCountPerProccess = (int *)malloc(sizeof(int) * size);
        int remainder = data->TCount % size;

        if (remainder > 0)
        {
            // assign each tCountPerProccess and add the to the last proccess
            for (int i = 0; i < size; i++)
            {
                tCountPerProccess[i] = data->TCount / size;
            }
            tCountPerProccess[size - 1] += remainder;
        }
        else
        {
            // distribute the tCountPerProccess equally
            for (int i = 0; i < size; i++)
            {
                tCountPerProccess[i] = data->TCount / size;
            }
        }
        //send data to each proccess
        for (int i = 1; i < size; i++)
        {
            MPI_Send((&data->N), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send((&data->K), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send((&data->D), 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            MPI_Send((&data->TCount), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send((data->points), sizeof(Point) * data->N, MPI_BYTE, i, 0, MPI_COMM_WORLD);
            MPI_Send((&tCountPerProccess[i]), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        // recive the data from the root
        int startCount, localTCount;
        MPI_Recv(&data->N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data->K, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data->D, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data->TCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Allocate memory for points array
        data->points = (Point *)malloc(sizeof(Point) * data->N);
        //check if points were allocated if not exit
        if (data->points == NULL)
        {
            printf("Memory allocation error for points array.\n");
            exit(1);
        }
        MPI_Recv(data->points, sizeof(Point) * data->N, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&localTCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank == size - 1)
            startCount = (data->TCount / size) * (size - 1);
        else
        {
            startCount = localTCount * rank;
        }
        int localMatrix[localTCount][3];
        //assing each cell with -1 in the local matrix
        for (int i = 0; i < localTCount; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                localMatrix[i][j] = -1;
            }
        }
        validate(data, localMatrix, startCount, localTCount, rank);
        // send the local matrix to the root 

        MPI_Send(&localMatrix, localTCount * 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        int idMatrix[data->TCount][3];
        for (int i = 0; i < data->TCount; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                idMatrix[i][j] = -1;
            }
        } 
       
        for (int i = 1; i < size; i++)
        { 
            // calculate the start and end count for each proccess the last proccess can have the remainder
            int startCount = (data->TCount / size) * i;
            int endCount = startCount + tCountPerProccess[i];
            int temp[tCountPerProccess[i] * 3];
            //Recive the local matrix from each proccess and save it in the global matrix in the correct position from the startCount to endCount
            MPI_Recv(&idMatrix[startCount], tCountPerProccess[i] * 3, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        validate(data, idMatrix, 0, tCountPerProccess[0], rank);
        WriteToOutputFile(OUTPUT_NAME, idMatrix, data->TCount);
        
    }

    
    MPI_Finalize();
    return 0;
}
