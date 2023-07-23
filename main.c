#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_POINTS 1000

typedef struct {
    int id;
    double x1, x2, a, b;
} Point;

typedef struct{
    int id;
    double x,y;
    int t;
}Coordinate;

void calculateCoordinates(int N,int TCount, Point points[], Coordinate coordinates[][N]) {
    double t;
    int id = 1;
    for (int i = 0; i < TCount; i++) {
        t = 2.0 * i / TCount - 1.0;
        
        
        for (int j = 0; j < N; j++) {
            double x = ((points[j].x2 - points[j].x1) / 2.0) * sin(t * M_PI / 2.0) + (points[j].x2 + points[j].x1) / 2.0;
            double y = points[j].a * x + points[j].b;
            
            coordinates[i][j].id = id;
            coordinates[i][j].x = x;
            coordinates[i][j].y = y;
            id++;
        }
       
    }
}

void readInputFile(const char* filename, int* N, int* K, float* D, int* TCount, Point points[]) {
    FILE* inputFile = fopen(filename, "r");
    if (inputFile == NULL) {
        printf("Failed to open the input file.\n");
        exit(1);
    }
    
    // Read the input parameters
    fscanf(inputFile, "%d %d %f %d", N, K, D, TCount);

    // Read the points data
    for (int i = 0; i < *N; i++) {
        fscanf(inputFile, "%d %lf %lf %lf %lf", &points[i].id, &points[i].x1, &points[i].x2, &points[i].a, &points[i].b);
    }
    
    fclose(inputFile);
}

void printCoordinates(int N, int TCount, Coordinate coordinates[][N]) {
    for (int i = 0; i < TCount; i++) {
        for (int j = 0; j < N; j++) {
            printf("Coordinate ID: %d, x: %.2f, y: %.2f\n", coordinates[i][j].id, coordinates[i][j].x, coordinates[i][j].y);
        }
    }
}

int checkDistanceSmallerThanD(double x1, double y1, double x2, double y2, float D) {
    double distance = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    return (distance < D) ? 1 : 0;
}


int** CheckRowsCoordinates(int N,float D,Coordinate* coordinate[],int TCount,int K){   
    int currentCol = 0;
    int currentRow = 0;
    int KCount = 0;
    int pointCount = 0;
  
    int** idMatrix=(int**)malloc(sizeof(int*)) ;
    Coordinate temp = coordinate[currentRow][currentCol];
    int array[3] = {};
    while(currentRow < TCount){
        
        for (int i = 0; i < N; i++)
        {
            /* code */
            if (currentCol != i)
            {
                Coordinate currentCord = coordinate[currentRow][i];
                KCount += checkDistanceSmallerThanD(temp.x,temp.y,currentCord.x,currentCord.y,D);
            }

            if(KCount == K){
                pointCount++;
                array[pointCount] = temp.id;
                break;
            }
        }
        
       
        

        

    }
   
    return idMatrix;
}


void printInputValues(int N, int K, float D, int TCount, Point points[]) {
    printf("N: %d\n", N);
    printf("K: %d\n", K);
    printf("D: %.2f\n", D);
    printf("TCount: %d\n", TCount);
    
    printf("Points:\n");
    for (int i = 0; i < N; i++) {
        printf("ID: %d, x1: %.2f, x2: %.2f, a: %.2f, b: %.2f\n",
               points[i].id, points[i].x1, points[i].x2, points[i].a, points[i].b);
    }
}

int main() {
    
    int N, K, TCount;
    float D;
    Point points[MAX_POINTS];
    
    readInputFile("input.txt", &N, &K, &D, &TCount, points);
    Coordinate coordinates[TCount][N];
    // Print the read values
    
    calculateCoordinates(N,TCount, points, coordinates);
    
    printCoordinates(N, TCount, coordinates);


    
    return 0;
}
