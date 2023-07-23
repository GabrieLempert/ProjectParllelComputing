#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int x1;
    int x2;
    int a;
    int b;
} PointParams;

typedef struct {
    int x;
    int y;
    int id;
} Point;

void readInputFile(const char* filename, int* N, int* K, double* D, int* TCount, PointParams** params) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Failed to open the input file.\n");
        return;
    }

    fscanf(file, "%d %d %lf %d\n", N, K, D, TCount);

    *params = (PointParams*)malloc((*N) * sizeof(PointParams));

    for (int i = 0; i < *N; i++) {
        fscanf(file, "%d %d %d %d\n", &((*params)[i].x1), &((*params)[i].x2), &((*params)[i].a), &((*params)[i].b));
    }

    fclose(file);
}

void writeOutputFile(const char* filename, double* tValues, int** pointIDs, int resultCount) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Failed to open the output file.\n");
        return;
    }

    if (resultCount > 0) {
        for (int i = 0; i < resultCount; i++) {
            fprintf(file, "Points  %d, %d, %d satisfy Proximity Criteria at t = %lf\n",
                pointIDs[i][0], pointIDs[i][1], pointIDs[i][2], tValues[i]);
        }
    } else {
        fprintf(file, "There were no 3 points found for any t.\n");
    }

    fclose(file);
}

double calculateX(int x1, int x2, double t) {
    return ((x2 - x1) / 2.0) * sin(t * M_PI / 2.0) + (x2 + x1) / 2.0;
}

int checkProximityCriteria(PointParams params, Point p, double D, int N, Point* points, int K) {
    int count = 0;
    for (int i = 0; i < N; i++) {
        if (i != p.id) {
            double distance = sqrt(pow(p.x - points[i].x, 2) + pow(p.y - points[i].y, 2));
            if (distance < D)
                count++;
            if (count >= K)
                return 1;
        }
    }
    return 0;
}

void findSatisfyingPoints(int N, PointParams* params, int K, double D, int TCount, Point* points, int* resultCount, double** tValues, int*** pointIDs) {
    *resultCount = 0;
    *pointIDs = NULL;
    *tValues = (double*)malloc(TCount * sizeof(double));

    for (int i = 0; i <= TCount; i++) {
        double t = 2.0 * i / TCount - 1.0;
        (*tValues)[i] = t;

        for (int j = 0; j < N - 2; j++) {
            for (int k = j + 1; k < N - 1; k++) {
                for (int l = k + 1; l < N; l++) {
                    Point p1, p2, p3;
                    p1.id = j;
                    p2.id = k;
                    p3.id = l;

                    p1.x = calculateX(params[j].x1, params[j].x2, t);
                    p2.x = calculateX(params[k].x1, params[k].x2, t);
                    p3.x = calculateX(params[l].x1, params[l].x2, t);

                    p1.y = params[j].a * p1.x + params[j].b;
                    p2.y = params[k].a * p2.x + params[k].b;
                    p3.y = params[l].a * p3.x + params[l].b;

                    if (checkProximityCriteria(params[j], p1, D, N, points, K) &&
                        checkProximityCriteria(params[k], p2, D, N, points, K) &&
                        checkProximityCriteria(params[l], p3, D, N, points, K)) {
                        (*resultCount)++;
                        *pointIDs = (int**)realloc(*pointIDs, (*resultCount) * sizeof(int*));
                        (*pointIDs)[(*resultCount) - 1] = (int*)malloc(3 * sizeof(int));
                        (*pointIDs)[(*resultCount) - 1][0] = p1.id;
                        (*pointIDs)[(*resultCount) - 1][1] = p2.id;
                        (*pointIDs)[(*resultCount) - 1][2] = p3.id;
                        break;
                    }
                }
                if (*resultCount > 0) break;
            }
            if (*resultCount > 0) break;
        }
    }
}

int main() {
    int N, K, TCount;
    double D;
    PointParams* params;
    Point* points;
    int resultCount;
    double* tValues;
    int** pointIDs;

    readInputFile("input.txt", &N, &K, &D, &TCount, &params);

    points = (Point*)malloc(N * sizeof(Point));
    for (int i = 0; i < N; i++) {
        points[i].id = i;
    }

    findSatisfyingPoints(N, params, K, D, TCount, points, &resultCount, &tValues, &pointIDs);

    writeOutputFile("output.txt", tValues, pointIDs, resultCount);

    free(params);
    free(points);
    free(tValues);
    for (int i = 0; i < resultCount; i++) {
        free(pointIDs[i]);
    }
    free(pointIDs);

    return 0;
}
