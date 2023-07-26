#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_POINTS 100000

typedef struct
{
    int id;
    double x1, x2, a, b;
} Point;

typedef struct
{
    int id;
    double x, y;
    int t;
} Coordinate;

void calculateCoordinates(int N, int TCount, Point points[], Coordinate coordinates[][N])
{
    double t;
    int id = 1;
    for (int i = 0; i < TCount; i++)
    {
        t = 2.0 * i / TCount - 1.0;

        for (int j = 0; j < N; j++)
        {
            double x = ((points[j].x2 - points[j].x1) / 2.0) * sin(t * M_PI / 2.0) + (points[j].x2 + points[j].x1) / 2.0;
            double y = points[j].a * x + points[j].b;

            coordinates[i][j].id = id;
            coordinates[i][j].x = x;
            coordinates[i][j].y = y;
            id++;
        }
    }
}

void readInputFile(const char *filename, int *N, int *K, float *D, int *TCount, Point points[])
{
    FILE *inputFile = fopen(filename, "r");
    if (inputFile == NULL)
    {
        printf("Failed to open the input file.\n");
        exit(1);
    }

    // Read the input parameters
    int result = fscanf(inputFile, "%d %d %f %d", N, K, D, TCount);
    if (result != 4)
    {
        printf("Error reading input parameters from the file.\n");
        fclose(inputFile);
        exit(1);
    }

    int size = *N;
    if (size <= 0 || size > MAX_POINTS)
    {
        printf("Invalid value for N in the input file.\n");
        fclose(inputFile);
        exit(1);
    }

    // Read the points data
    for (int i = 0; i < size; i++)
    {
        result = fscanf(inputFile, "%d %lf %lf %lf %lf", &points[i].id, &points[i].x1, &points[i].x2, &points[i].a, &points[i].b);
        if (result != 5)
        {
            printf("Error reading data for Point %d from the file.\n", i + 1);
            fclose(inputFile);
            exit(1);
        }
    }

    fclose(inputFile);
}

void printCoordinates(int N, int TCount, Coordinate coordinates[][N])
{
    for (int i = 0; i < TCount; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("Coordinate ID: %d, x: %.2f, y: %.2f\n", coordinates[i][j].id, coordinates[i][j].x, coordinates[i][j].y);
        }
    }
}

int checkDistanceSmallerThanD(double x1, double y1, double x2, double y2, float D)
{
    double distance = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    return (distance < D) ? 1 : 0;
}

int **CheckRowsCoordinates(int N, float D, Coordinate coordinate[][N], int TCount, int K)
{
    int currentCol = 0;
    int currentRow = 0;
    int KCount = 0;
    int pointCount = 0;                                    
    int **idMatrix = (int **)malloc(TCount * sizeof(int *)); // Allocate memory for all rows at once

    // Initialize idMatrix with all zeroes using calloc
    for (int i = 0; i < TCount; i++)
    {
        idMatrix[i] = (int *)calloc(3, sizeof(int));
    }

    while (currentRow < TCount)
    {
        Coordinate temp = coordinate[currentRow][currentCol];

        // Loops through the point and checks if it's valid
        for (int i = 0; i < N; i++)
        {
            if (currentCol != i)
            {
                Coordinate currentCord = coordinate[currentRow][i];
                KCount += checkDistanceSmallerThanD(temp.x, temp.y, currentCord.x, currentCord.y, D);
            }
            if (KCount == K)
            {
                pointCount++;
                // Store the point ID in the array
                idMatrix[currentRow][pointCount] = temp.id;
                // Reset KCount for the next point
                KCount = 0;
            }
        }

        currentCol++;
        if (currentCol == N)
        {
            // Move to the next row
            currentRow++;
            currentCol = 0;
            pointCount = 0; // Reset pointCount for the new row
        }
    }

    return idMatrix;
}

void printInputValues(int N, int K, float D, int TCount, Point points[])
{
    printf("N: %d\n", N);
    printf("K: %d\n", K);
    printf("D: %.2f\n", D);
    printf("TCount: %d\n", TCount);

    printf("Points:\n");
    for (int i = 0; i < N; i++)
    {
        printf("ID: %d, x1: %.2f, x2: %.2f, a: %.2f, b: %.2f\n",
               points[i].id, points[i].x1, points[i].x2, points[i].a, points[i].b);
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

        if (id1 > 0 && id2 > 0 && id3 > 0)
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
void PrintMatrix(int **matrix, int rows, int cols)
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

int main()
{

    int N, K, TCount;
    float D;
    Point points[MAX_POINTS];

    readInputFile("input.txt", &N, &K, &D, &TCount, points);
    Coordinate coordinates[TCount][N];
    // Print the read values

    calculateCoordinates(N, TCount, points, coordinates);

    printCoordinates(N, TCount, coordinates);

    int **idMatrix = CheckRowsCoordinates(N, D, coordinates, TCount, K);
    PrintMatrix(idMatrix, TCount, 3);
    WriteToOutputFile("output.txt", idMatrix, TCount);

    // Don't forget to free the memory allocated for idMatrix
    for (int i = 0; i < TCount; i++)
    {
        free(idMatrix[i]);
    }
    free(idMatrix);

    return 0;
}
