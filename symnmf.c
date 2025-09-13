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
    int point_count = 0, dim = 0;
    double **matrix;
    int i, j;

    file = fopen(file_name, "r");
    if (file == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    /* count the number of points and dimension */
    while (fscanf(file, "%lf%c", &value, &c) != EOF) { /* %lf is for floating point number and %c is for delimeter (comma or newline) */
        if (point_count == 0) { dim++; } /* first point, count dimension */
        if (c == '\n') { point_count++; } /* end of point */
    }
    if (c != '\n' && ftell(file) > 0) { point_count++; } /* in case no newline is in last line */
    *N = point_count; *d = dim; /* used to access point_count and dim globally */
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
            for (j = 0; j < i; j++) { free(matrix[j]); } /* free previously allocated rows */
            free(matrix); fclose(file); exit(1);
        }
    }
    for (i = 0; i < point_count; i++) { /* populate matrix*/
        for (j = 0; j < dim; j++) { fscanf(file, "%lf%c", &matrix[i][j], &c); }
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
                printf(",");
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
        return NULL; 
    }
    for (i = 0; i < N; i++) {
        ddg_matrix[i] = (double *)calloc(N, sizeof(double)); /* we use calloc for the 0 initialization */
        if (ddg_matrix[i] == NULL) { 
            free_matrix(sym_matrix, N);
            free_matrix(ddg_matrix, i);
            return NULL; 
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
    int i, j;
    double* degrees;
    double **sym_matrix, **norm_matrix;

    /* get similarity matrix A */
    sym_matrix = calculate_sym_matrix(data_points, N, d);
    if (sym_matrix == NULL) { return NULL; } /* caller is the handler */
    /* calculate array degrees, in which every entry is a row sum */
    degrees = (double *)calloc(N, sizeof(double)); /* again for 0 initialization */
    if (degrees == NULL) {
        free_matrix(sym_matrix, N);
        return NULL; 
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) { degrees[i] += sym_matrix[i][j]; }
    }
    /* allocate memory for W */
    norm_matrix = (double **)malloc(N * sizeof(double *));
    if (norm_matrix == NULL) {
        free_matrix(sym_matrix, N);
        free(degrees);
        return NULL; 
    }
    for (i = 0; i < N; i++) {
        norm_matrix[i] = (double *)malloc(N * sizeof(double));
        if (norm_matrix[i] == NULL) {
            free_matrix(sym_matrix, N);
            free(degrees);
            free_matrix(norm_matrix, i);
            return NULL; 
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
    free_matrix(sym_matrix, N); free(degrees); /* free un needed data */
    return norm_matrix;
}
/*
 * ========================================MAT_MULTIPLY========================================
 * helper method to multiply two matrixes A and B
*/
double** mat_multiply(double** A, double** B, int rowsA, int colsA, int colsB) {
    double **return_matrix;
    int i;
    int j;
    int k;

    /* allocation */
    return_matrix = (double **)calloc(rowsA, sizeof(double *));
    if (return_matrix == NULL) {
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < rowsA; i++) {
        return_matrix[i] = (double *)calloc(colsB, sizeof(double));
        if (return_matrix[i] == NULL) {
            for (j = 0; j < i; j++) {
                free(return_matrix[j]);
            }
            free(return_matrix);
            return NULL; 
        }
    }
    /* multiplication */
    for (i = 0; i < rowsA; i++) {
        for (j = 0; j < colsB; j++) {
            for (k = 0; k < colsA; k++) {
                return_matrix[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return return_matrix;
}
/*
 * ========================================MAT_TRANSPOSE===========================================
 * helper method to transpose a matrix A
*/
double** mat_transpose(double** A, int rows, int cols) {
    double **return_matrix;
    int i;
    int j;

    /* allocation */
    return_matrix = (double **)malloc(cols * sizeof(double *)); /* rows and cols switch places */
    if (return_matrix == NULL) {
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < cols; i++) {
        return_matrix[i] = (double *)malloc(rows * sizeof(double));
        if (return_matrix[i] == NULL) {
            for (j = 0; j < i; j++) {
                free(return_matrix[j]);
            }
            free(return_matrix);
            return NULL;
        }
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            return_matrix[j][i] = A[i][j]; /* switch places */
        }
    }
    return return_matrix;
}
/*
 * =======================================FROBENIUS_NORM_SQUARED===================================
 * helper method to calculate the frobenius norm squared of two matrixes A and B
*/
double frobenius_norm_squared(double** A, double** B, int rows, int cols) {
    double norm = 0.0;
    int i;
    int j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            norm += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return norm;
}
/*
 * ===========================================OPTIMIZE_H===========================================
 * this method does the core optimization of the algorithm, iteratively.
 * the matrix operations are left for the helper methods
 */
double** optimize_h(double** W, double** init_H, int N, int k) {
    const int max_iter = 300;
    const double eps = 1e-4;
    const double beta = 0.5;
    double **curr_H, **next_H, **curr_H_transpose, **denominator, **numerator, **temp_matrix;
    int i, j, iter;
    double curr_frobenius_norm;

    /* allocate memory for curr_H, next_H */
    curr_H = (double **)malloc(N * sizeof(double *));
    if (curr_H == NULL) { return NULL; }
    for (i = 0; i < N; i++) {
        curr_H[i] = (double *)malloc(k * sizeof(double));
        if (curr_H[i] == NULL) {
            for (j = 0; j < i; j++) { free(curr_H[j]); }
            free(curr_H); return NULL; } /* caller is the handler */
        for (j = 0; j < k; j++) { curr_H[i][j] = init_H[i][j]; } /* copy initial H from the argument */
    }
    next_H = (double **)malloc(N * sizeof(double *));
    if (next_H == NULL) { free_matrix(curr_H, N); return NULL; }
    for (i = 0; i < N; i++) {
        next_H[i] = (double *)malloc(k * sizeof(double));
        if (next_H[i] == NULL) {
            for (j = 0; j < i; j++) { free(next_H[j]); }
            free(next_H); free_matrix(curr_H, N); return NULL; } 
    }
    /* optimization loop */
    for (iter = 0; iter < max_iter; iter++) {
        curr_H_transpose = mat_transpose(curr_H, N, k); /* get H^T, which is allocated in the helper method */
        if (curr_H_transpose == NULL) { free_matrix(curr_H, N); free_matrix(next_H, N); return NULL; } 
        temp_matrix = mat_multiply(curr_H, curr_H_transpose, N, k, N); /* calculate denominator */
        if (temp_matrix == NULL) { free_matrix(curr_H, N); free_matrix(next_H, N); free_matrix(curr_H_transpose, k); return NULL; } 
        denominator = mat_multiply(temp_matrix, curr_H, N, N, k);
        if (denominator == NULL) { free_matrix(curr_H, N); free_matrix(next_H, N); free_matrix(curr_H_transpose, k); free_matrix(temp_matrix, N); return NULL; } 
        free_matrix(temp_matrix, N); /* we dont need temp_matrix anymore */
        numerator = mat_multiply(W, curr_H, N, N, k); /* calculate numerator */
        if (numerator == NULL) { free_matrix(curr_H, N); free_matrix(next_H, N); free_matrix(curr_H_transpose, k); free_matrix(denominator, N); return NULL; } 
        for (i = 0; i < N; i++) { /* update next_H */
            for (j = 0; j < k; j++) {
                if (denominator[i][j] == 0) { next_H[i][j] = curr_H[i][j]; } /* avoid division by zero */
                else { next_H[i][j] = curr_H[i][j] * (1 - beta + beta * (numerator[i][j] / denominator[i][j])); }
            }
        }
        curr_frobenius_norm = frobenius_norm_squared(next_H, curr_H, N, k); /* check convergence */
        if (curr_frobenius_norm < eps) { free_matrix(curr_H, N); free_matrix(curr_H_transpose, k); free_matrix(denominator, N); free_matrix(numerator, N); return next_H; } /* converged */
        for (i = 0; i < N; i++) { /* swap curr_H and next_H if no convergence */
            for (j = 0; j < k; j++) { curr_H[i][j] = next_H[i][j]; }
        }
        for (i = 0; i < N; i++) { /* reset next_H for the next iteration */
            for (j = 0; j < k; j++) { next_H[i][j] = 0.0; }
        }
        free_matrix(curr_H_transpose, k); free_matrix(denominator, N); free_matrix(numerator, N); /* free intermediate matrices */
    }
    free_matrix(curr_H, N); free_matrix(next_H, N); /* if we reach here, we did not converge - then we should free curr_H, next_H and return*/
    return NULL; 
}
/*
 * ================================================================================================
 * ==============================================MAIN==============================================
 * ================================================================================================
*/
int main(int argc, char *argv[]) {
    char *goal; char *filename; int N = 0; int d = 0; double **data_points = NULL; double **sym_matrix = NULL; double **ddg_matrix = NULL; double **norm_matrix = NULL;
    if (argc != 3) { printf("An Error Has Occurred\n"); return 1; }
    goal = argv[1]; filename = argv[2];
    data_points = read_data_points(filename, &N, &d); if (data_points == NULL) { printf("An Error Has Occurred\n"); return 1; }
    if (string_compare(goal, "sym") == 1) {
        sym_matrix = calculate_sym_matrix(data_points, N, d);
        if (sym_matrix == NULL) {
            printf("An Error Has Occurred\n");
            free_matrix(data_points, N);
            return 1;
        }
        print_matrix(sym_matrix, N, N);
        free_matrix(sym_matrix, N);
    }
    else if (string_compare(goal, "ddg") == 1) {
        ddg_matrix = calculate_ddg_matrix(data_points, N, d);
        if (ddg_matrix == NULL) {
            printf("An Error Has Occurred\n");
            free_matrix(data_points, N);
            return 1;
        }
        print_matrix(ddg_matrix, N, N);
        free_matrix(ddg_matrix, N);
    }
    else if (string_compare(goal, "norm") == 1) {
        norm_matrix = calculate_norm_matrix(data_points, N, d);
        if (norm_matrix == NULL) {
            printf("An Error Has Occurred\n");
            free_matrix(data_points, N);
            return 1;
        }
        print_matrix(norm_matrix, N, N);
        free_matrix(norm_matrix, N);
    }
    else {
        printf("An Error Has Occurred\n");
        return 1; /* invalid goal */
    }
    /* free allocated memory */
    free_matrix(data_points, N); /* free original matrix */
    return 0;
}
