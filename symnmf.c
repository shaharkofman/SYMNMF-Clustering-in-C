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
/*
 * ===================================CALCULATE_SYM_MATRIX=========================================
 * allocates memofy for matrix A of size NxN and populates it using formula 1.1
*/
double** calculate_sym_matrix(double **data_points, int N, int d) {
    int i, j;
    double dist;
    double **sym_matrix = (double **)malloc(N * sizeof(double *));
    if (sym_matrix == NULL) {
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < N; i++) {
        sym_matrix[i] = (double *)malloc(N * sizeof(double));
        if (sym_matrix[i] == NULL) { /* mid-allocation error, free previous allocations */
            for (j = 0; j < i; j++) {
                free(sym_matrix[j]);
            }
            free(sym_matrix);
            return NULL; /* caller is the handler */
        }
        for (j = 0; j < N; j++) {
            if (i == j) { /* same point */
                sym_matrix[i][j] = 0;
            }
            else {
                dist = squared_euclidean_distance(data_points[i], data_points[j], d);
                sym_matrix[i][j] = exp(-dist / 2.0);
            }
        }
    }
    return sym_matrix;
}
/*
 * ========================================FREE_MATRIX=============================================
 * free allocated memory for a 2D matrix of size NxN
*/
void free_matrix(double **matrix, int N) {
    int i;
    if (matrix != NULL) {
        for (i = 0; i < N; i++) {
            free(matrix[i]); 
        }
        free(matrix); 
    }
}
/*
 * ========================================PRINT_MATRIX============================================
 * print matrix with 4 decimal places format and comma seperation
*/
void print_matrix(double **matrix, int rows, int cols) {
    int i;
    int j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1) {
                printf(", ");
            }
        }
        printf("\n");
    }
}
/*
 * ========================================CALCULATE_DDG_MATRIX====================================
 * calculate the diagonal degree matrix D from the symmetric matrix A
 * this method reuses calculate_sym_matrix() to calculate the symmetric matrix A
*/
double** calculate_ddg_matrix(double **data_points, int N, int d) {
    int i;
    int j;
    double row_sum;
    double **sym_matrix;
    double **ddg_matrix;

    /* get similarity matrix A */
    sym_matrix = calculate_sym_matrix(data_points, N, d);
    if (sym_matrix == NULL) {
        return NULL; /* caller is the handler */
    }

    /* allocate memort for ddg matrix D */
    ddg_matrix = (double **)malloc(N * sizeof(double *));
    if (ddg_matrix == NULL) {
        free_matrix(sym_matrix, N); /* free intermediate matrix */
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < N; i++) {
        ddg_matrix[i] = (double *)calloc(N, sizeof(double)); /* we use calloc for the 0 initialization */
        if (ddg_matrix[i] == NULL) { 
            free_matrix(sym_matrix, N);
            free_matrix(ddg_matrix, i);
            return NULL; /* caller is the handler */
        }
        /* calculate row sum for each row, place the sum on the diagonal of that row */
        row_sum = 0.0;
        for (j = 0; j < N; j++) {
            row_sum += sym_matrix[i][j];
        }
        ddg_matrix[i][i] = row_sum;
    }
    free_matrix(sym_matrix, N); /* we dont need matrix A */

    return ddg_matrix;
}
/*
 * ========================================CALCULATE_NORM_MATRIX===================================
 * calculates normalized similarity matrix W
 * this method reuses calculate_sym_matrix() to calculate the symmetric matrix A
*/
double** calculate_norm_matrix(double **data_points, int N, int d) {
    int i;
    int j;
    double* degrees;
    double **sym_matrix;
    double **norm_matrix;

    /* get similarity matrix A */
    sym_matrix = calculate_sym_matrix(data_points, N, d);
    if (sym_matrix == NULL) {
        return NULL; /* caller is the handler */
    }
    /* calculate array degrees, in which every entry is a row sum */
    degrees = (double *)calloc(N, sizeof(double)); /* again for 0 initialization */
    if (degrees == NULL) {
        free_matrix(sym_matrix, N);
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            degrees[i] += sym_matrix[i][j];
        }
    }
    /* allocate memory for W */
    norm_matrix = (double **)malloc(N * sizeof(double *));
    if (norm_matrix == NULL) {
        free_matrix(sym_matrix, N);
        free(degrees);
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < N; i++) {
        norm_matrix[i] = (double *)malloc(N * sizeof(double));
        if (norm_matrix[i] == NULL) {
            free_matrix(sym_matrix, N);
            free(degrees);
            free_matrix(norm_matrix, i);
            return NULL; /* caller is the handler */
        }
        /* to calculate entry i,j of W we simply calculate A_ij/(sqrt(degree_i) * sqrt(degree_j)) */
        for (j = 0; j < N; j++) {
            if (degrees[i] == 0 || degrees[j] == 0) {
                norm_matrix[i][j] = 0; /* avoid division by zero */
            }
            else {
                norm_matrix[i][j] = sym_matrix[i][j] / (sqrt(degrees[i]) * sqrt(degrees[j]));
            }
        }
    }
    /* free un needed data */
    free_matrix(sym_matrix, N);
    free(degrees);

    return norm_matrix;
}
