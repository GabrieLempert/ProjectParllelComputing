#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_POINTS 100000
#define FILE_NAME "input.txt"
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
void CalculateCoordinate(Point *p, int t, Coordinate *coord, int N,int index);
int CheckCretirieaByT(init_p data, Coordinate *coord, int index, int t);

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

void CalculateCoordinate(Point *p, int t, Coordinate *coord, int N,int index)
{
    for (int j = 0; j < N; j++)
    {
        double x = ((p[j].x2 - p[j].x1) / 2.0) * sin(t * M_PI / 2.0) + (p[j].x2 + p[j].x1) / 2.0;
        double y = p[j].a * x + p[j].b;

        coord[j].id = p[j].id+N*index;
        coord[j].x = x;
        coord[j].y = y;
    }
}

int checkDistanceSmallerThanD(double x1, double y1, double x2, double y2, float D)
{
    double distance = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    return (distance < D) ? 1 : 0;
}

int CheckCretirieaByT(init_p data, Coordinate *coord, int index, int t)
{
    Coordinate temp = coord[index];
    int K_approximates = data.K;
    float D_maximum = data.D;
    int count = 0;
    for (int i = 0; i < data.N; i++)
    {
        if (index != i)
        {
            count += checkDistanceSmallerThanD(temp.x, temp.y, coord[i].x, coord[i].y, D_maximum);
        }
    }
    if (count >= K_approximates)
    {
        return 1;
    }

    return 0;
}

void validate(init_p data, int **idMatrix)
{
    Coordinate *coord = (Coordinate *)malloc(sizeof(Coordinate) * data.N);
    int validationCount = 0;

    for (int i = 0; i < data.TCount; i++)
    {
        int t = 2.0 * i / data.TCount - 1;
        // calculate coordinate for specif t each time
        CalculateCoordinate(data.points, t, coord, data.N,i);
        for(int j=0;j<data.N;j++){
            if (validationCount==3)
            {
                break;
            }
            
            if(CheckCretirieaByT(data,coord,j,t)){
                idMatrix[i][validationCount] = coord[j].id;
                validationCount++;
            }
        }
        
        validationCount = 0;
        
    }
}
// Section for debugging.
void printInputValues(init_p data)
{
    printf("N: %d\n", data.N);
    printf("K: %d\n", data.K);
    printf("D: %.2f\n", data.D);
    printf("TCount: %d\n", data.TCount);

    printf("Points:\n");
    for (int i = 0; i < data.N; i++)
    {
        printf("ID: %d, x1: %.2f, x2: %.2f, a: %.2f, b: %.2f\n",
               data.points[i].id, data.points[i].x1, data.points[i].x2, data.points[i].a, data.points[i].b);
    }
}
int **initializeIdMatrix(int rows)
{
    int **matrix = (int **)malloc(rows * sizeof(int *));
    if (matrix == NULL)
    {
        printf("Memory allocation error for idMatrix.\n");
        exit(1);
    }

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = (int *)malloc(3 * sizeof(int));
        if (matrix[i] == NULL)
        {
            printf("Memory allocation error for idMatrix row %d.\n", i);
            exit(1);
        }

        for (int j = 0; j < 3; j++)
        {
            matrix[i][j] = -1;
        }
    }

    return matrix;
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
       int id1 =idMatrix[i][0];
       int id2= idMatrix[i][1];
       int id3 = idMatrix[i][2];

        if (id3>0)
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

int main()
{

    init_p data;
    initData(&data, FILE_NAME);
    // printInputValues(data);
    int** idMatrix = initializeIdMatrix(data.TCount);
    validate(data,idMatrix);
    WriteToOutputFile("output.txt",idMatrix,data.TCount);




    // Don't forget to free the memory allocated for idMatrix

    return 0;
}