#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
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
void freeLocalIdMatrix(int **matrix, int rows);
void initData(init_p *data, const char *filename);
int checkDistanceSmallerThanD(double x1, double y1, double x2, double y2, float D);
void CalculateCoordinate(Point *p, int t, Coordinate *coord, int N, int index);
int CheckCretirieaByT(init_p *data, Coordinate *coord, int index);

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
void CalculateCoordinate(Point *p, int t, Coordinate *coord, int N, int index)
{
    double x = ((p[index].x2 - p[index].x1) / 2.0) * sin(t * M_PI / 2.0) + (p[index].x2 + p[index].x1) / 2.0;
    double y = p[index].a * x + p[index].b;

    coord[index].id = p[index].id + N * index;
    coord[index].x = x;
    coord[index].y = y;
}
void printInputValues(init_p *data, int rank)
{
    printf("My rank is,%d\n", rank);
    printf("N: %d\n", data->N);
    printf("K: %d\n", data->K);
    printf("D: %.2f\n", data->D);
    printf("TCount: %d\n", data->TCount);

    printf("Points:\n");
    for (int i = 0; i < data->N; i++)
    {
        printf("ID: %d, x1: %.2f, x2: %.2f, a: %.2f, b: %.2f\n",
               data->points[i].id, data->points[i].x1, data->points[i].x2, data->points[i].a, data->points[i].b);
    }
}

int CheckCretirieaByT(init_p *data, Coordinate *coord, int index)

{
    Coordinate temp = coord[index];
    int K_approximates = data->K;
    float D_maximum = data->D;
    int count = 0;

    for (int i = 0; i < data->N; i++)
    {
        if (index != i)
        {
            count += checkDistanceSmallerThanD(temp.x, temp.y, coord[i].x, coord[i].y, D_maximum);
        }

        if (count >= K_approximates)
        {
            return 1;
        }
    }

    return 0;
}
void validate(init_p *data, int **idMatrix, int startTCount, int endTcount, int size, int rank)
{
    printf("My rank is %d", rank);
    Coordinate *coord = (Coordinate *)malloc(sizeof(Coordinate) * data->N);
    int validationCount = 0;
    printf("Checking,%d", rank);
    for (int i = 0; i < size; i++)
    {
        int calculateFrom = startTCount + i;
        int t = 2.0 * calculateFrom / data->TCount - 1;
        // calculate coordinate for specif t each time
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
void initializeIdMatrix(int **localMatrix)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            localMatrix[i][j] = -1;
        }
    }
}
void WriteToOutputFile(const char *filename, int **idMatrix, int TCount)
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
        fprintf(outputFile, "Aren't any points");
    }

    fclose(outputFile);
}
void freeLocalIdMatrix(int **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}
MPI_Datatype createPointType()
{
    MPI_Datatype pointType;
    MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int blockLengths[5] = {1, 1, 1, 1, 1};
    MPI_Aint displacements[5];

    Point samplePoint;
    MPI_Get_address(&samplePoint.id, &displacements[0]);
    MPI_Get_address(&samplePoint.x1, &displacements[1]);
    MPI_Get_address(&samplePoint.x2, &displacements[2]);
    MPI_Get_address(&samplePoint.a, &displacements[3]);
    MPI_Get_address(&samplePoint.b, &displacements[4]);

    for (int i = 1; i < 5; i++)
    {
        displacements[i] -= displacements[0];
    }
    displacements[0] = 0;

    MPI_Type_create_struct(5, blockLengths, displacements, types, &pointType);
    MPI_Type_commit(&pointType);

    return pointType;
}

MPI_Datatype createInitPType()
{
    MPI_Datatype initPType;
    MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT};
    int blockLengths[4] = {1, 1, 1, 1};
    MPI_Aint displacements[4];
    init_p sampleInitP; // Create a sample instance of init_p struct
    MPI_Get_address(&sampleInitP.N, &displacements[0]);
    MPI_Get_address(&sampleInitP.K, &displacements[1]);
    MPI_Get_address(&sampleInitP.TCount, &displacements[2]);
    MPI_Get_address(&sampleInitP.D, &displacements[3]);

    for (int i = 1; i < 4; i++)
    {
        displacements[i] -= displacements[0];
    }
    displacements[0] = 0;

    MPI_Type_create_struct(4, blockLengths, displacements, types, &initPType);
    MPI_Type_commit(&initPType);

    return initPType;
}
void printPoints(Point *points, int N, int rank)
{
    printf("Rank %d Points:\n", rank);
    for (int i = 0; i < N; i++)
    {
        printf("Point %d: ID: %d, x1: %.2f, x2: %.2f, a: %.2f, b: %.2f\n",
               i + 1, points[i].id, points[i].x1, points[i].x2, points[i].a, points[i].b);
    }
    printf("\n");
}
void printIntMatrix(int **matrix, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    int rank, size;
    int **idMatrix;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
    {
        init_p data;
        // Allocate memory for the data struct

        initData(&data, FILE_NAME);
        // initializeIdMatrix(data->TCount);
        for (int i = 1; i < size; i++)
        {
            MPI_Send((&data.N), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send((&data.K), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send((&data.D), 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            MPI_Send((data.TCount), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send((data.points), sizeof(Point) * data.N, MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        init_p data;
        MPI_Recv(&data.N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data.K, 1, MPI_INT,0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data.D, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data.TCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Allocate memory for points array
        data.points = (Point *)malloc(sizeof(Point) * data.N);

        MPI_Recv(data.points, sizeof(Point) * data.N, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

       
    }

    // initializeIdMatrix(rank);

    // Broadcast the data structure to all processes
    MPI_Bcast(data, 1, MPIData, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    // Broadcast the points array from rank 0 to all other processes
    MPI_Bcast(data->points, data->N * sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);
    // Calculate the range of T values for each process
    int localTCount, startTCount, endTcount;
    int rangePerProcess = data->TCount / size;
    int remainder = data->TCount % size;
    if (rank < remainder)
    {
        startTCount = rank * (rangePerProcess + 1);
        endTcount = startTCount + rangePerProcess + 1;
    }
    else
    {
        startTCount = remainder * (rangePerProcess + 1) + (rank - remainder) * rangePerProcess;
        endTcount = startTCount + rangePerProcess;
    }
    localTCount = endTcount - startTCount;
    // Calculate the coordinates and validate for the local range
    printf("Test");
    int localMatrix[localTCount][3];
    printf("Test2");

    initializeIdMatrix(&localMatrix);
    printIntMatrix(&localMatrix, localTCount, 3);

    // validate(data,localMatrix, startTCount, endTcount,localTCount,rank);
    //  Gather idMatrix data from all processes to rank 0
    // MPI_Gather(localIdMatrix[0], 3 * localTCount, MPI_INT, idMatrix[0], 3 * localTCount, MPI_INT, 0, MPI_COMM_WORLD);

    // Free memory and finalize MPI
    if (rank == 0)
    {
        free(data->points);
        free(data);
    }
    else
    {
        free(data->points);
        free(data);
    }

    MPI_Finalize();

    return 0;
}
