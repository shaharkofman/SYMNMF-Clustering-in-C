#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"

/*
 * ========================================STRING_COMPARE==========================================
 * a function for comparing strings
 * returns 1 if the strings are equal, 0 otherwise
 */
int string_compare(const char *str1, const char *str2) {
    while (*str1 && *str2) {
        if (*str1 != *str2) {
            return 0; 
        }
        str1++;
        str2++;
    }
    return (*str1 == *str2); 
}
/*
 * ========================================READ_DATA_POINTS========================================
 * reads data points from file into a dynamically allocated 2D array
 * returns a pointer to the array and sets N (number of points) and d (dimension)
 */
double** read_data_points(char *file_name, int *N, int *d) {
    FILE *file;
    double value;
    char c;
    int point_count = 0;
    int dim = 0;
    double **matrix;
    int i, j;

    file = fopen(file_name, "r");
    if (file == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    /* count the number of points and dimension */
    while (fscanf(file, "%lf%c", &value, &c) != EOF) { /* %lf is for floating point number and %c is for delimeter (comma or newline) */
        if (point_count == 0) {
            dim++; /* first point, count dimension */
        }
        if (c == '\n') {
            point_count++; /* end of point */
        }
    }
    /* in case no newline is in last line */
    if (c != '\n' && ftell(file) > 0) {
        point_count++;
    }
    *N = point_count; /* used to access point_count globally */
    *d = dim; /* used to access dim globally */

    rewind(file); /* reset file pointer to the beginning */

    /* allocate memory for the matrix */
    matrix = (double **)malloc(point_count * sizeof(double *));
    if (matrix == NULL) {
        printf("An Error Has Occurred\n");
        fclose(file);
        exit(1);
    }
    for (i = 0; i < point_count; i++) {
        matrix[i] = (double *)malloc(dim * sizeof(double));
        if (matrix[i] == NULL) {
            printf("An Error Has Occured\n");
            for (j = 0; j < i; j++) {
                free(matrix[j]); /* free previously allocated rows */
            }
            free(matrix);
            fclose(file);
            exit(1);
        }
    }
    /* populate matrix*/
    for (i = 0; i < point_count; i++) {
        for (j = 0; j < dim; j++) {
            fscanf(file, "%lf%c", &matrix[i][j], &c);
        }
    }

    fclose(file);
    return matrix;
}
/*
 * ========================================EUCLIDEAN_DISTANCE=====================================
 * calculate squared euclidean distance between two vectors
*/
double squared_euclidean_distance(double *vec1, double *vec2, int d) {
    double dist = 0.0;
    int i;
    for (i = 0; i < d; i++) {
        dist += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return dist;
}
