#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

static PyMethodDef symnmf_methods[] = {
    {"sym", (PyCFunction)sym_capi, METH_VARARGS, "calculate similarity matrix A"},
    {"ddg", (PyCFunction)ddg_capi, METH_VARARGS, "calculate diagonal degree matrix D"},
    {"norm", (PyCFunction)norm_capi, METH_VARARGS, "calculates normalized similarity matrix W"},
    {"symnmf", (PyCFunction)symnmf_capi, METH_VARARGS, "execute the full symnmf algorithm"},
    {NULL, NULL, 0, NULL} 
};
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",                                    
    NULL,
    -1,
    symnmf_methods
};
PyMODINIT_FUNC PyInit_symnmf(void) {
    return PyModule_Create(&symnmfmodule);
}

/*
===============================================PY_TO_C_MATRIX=========================================
 * similarly to how we read input, this function takes python objects and builds the c matrix 
 * to reduce code repetition since a lot of functions in the module need this conversion
*/
static double** py_to_c_matrix(PyObject* python_points_list, int *N, int *d) {
    double **data_points;
    int i;
    int j;
    double c_coord;

    if (!PyList_Check(python_points_list)) {
        PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
        return NULL;
    }
    /* parse number of points and their dimention */
    *N = PyList_Size(python_points_list);
    *d = PyList_Size(PyList_GetItem(python_points_list, 0)); /* length of first vector */

    /* now we build a C 2D array with the python object */
    data_points = (double **)malloc(*N * sizeof(double *));
    if (data_points == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; /* caller is the handler */
    }
    for (i = 0; i < *N; i++) {
        data_points[i] = (double *)malloc(*d * sizeof(double));
        if (data_points[i] == NULL) {
            for (j = 0; j < i; j++) {
                free(data_points[j]);
            }
            free(data_points);
            PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
            return NULL; 
        }
        PyObject *point = PyList_GetItem(python_points_list, i);
        if (point == NULL) {
            free_matrix(data_points, *N);
            return NULL; /* error is raised by parsing function */
        }
        for (j = 0; j < *d; j++) {
            PyObject *py_coord = PyList_GetItem(point, j);
            if (py_coord == NULL) {
                free_matrix(data_points, *N);
                return NULL; /* error is raised by parsing function */
            }

            c_coord = PyFloat_AsDouble(py_coord); /* actual conversion */
            /* if PyFloat() returns -1.0 then parsing has failed but we also need to consider that the coordinate itself can be -1.0 */
            if (c_coord == -1.0 && PyErr_Occurred()) {
                free_matrix(data_points, *N);
                return NULL; /* error is raised by parsing function */
            }
            data_points[i][j] = c_coord;
        }
    }
    return data_points;
}
/*
 * =======================================C_TO_PY_MATRIX============================================
 * this function converts a C matrix back to a python object
*/
static PyObject* c_to_py_matrix(double **matrix, int N, int d) {
    int i;
    int j;
    PyObject *py_matrix;
    PyObject *py_row;
    PyObject *py_coord;

    py_matrix = PyList_New(N);
    if (py_matrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        /* no need to dereference anything here */
        return NULL; /* error is raised by parsing function */
    }
    for (i = 0; i < N; i++) {
        py_row = PyList_New(d);
        if (py_row == NULL) {
            PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
            Py_DECREF(py_matrix); /* dereference the matrix */
            return NULL; 
        }
        for (j = 0; j < d; j++) {
            py_coord = PyFloat_FromDouble(matrix[i][j]);
            if (py_coord == NULL) {
                PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
                Py_DECREF(py_row); /* dereference the row we created */
                Py_DECREF(py_matrix);
                return NULL;
            }
            PyList_SetItem(py_row, j, py_coord); /* put the i,j coordinate in the current python row */
        }
        PyList_SetItem(py_matrix, i, py_row); /* put the completed row in the python matrix */
    }
    return py_matrix;
}
/*
 * ========================================SYM_CAPI=================================================
 * this function is the C API for calling calculate_sym_matrix from python
*/
PyObject* sym_capi(PyObject *self, PyObject *args) {
    PyObject *python_points_list;
    int N;
    int d;
    double **data_points;
    double **sym_matrix;
    PyObject *py_return_matrix;

    /* parse the python object */
    if (!PyArg_ParseTuple(args, "O", &python_points_list)) {
        return NULL; /* error is raised by parsing function */
    }
    /* convert python object to C matrix using the helper method */
    data_points = py_to_c_matrix(python_points_list, &N, &d);
    if (data_points == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    /* calculate the symmetric matrix */
    sym_matrix = calculate_sym_matrix(data_points, N, d);
    free_matrix(data_points, N); /* we dont need the data points anymore */
    if (sym_matrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; 
    }
    /* convert the C matrix to a python object */
    py_return_matrix = c_to_py_matrix(sym_matrix, N, N);
    free_matrix(sym_matrix, N); /* we dont need the symetric matrix anymore */
    if (py_return_matrix == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    return py_return_matrix; 
}
/*
 * ========================================DDG_CAPI=================================================
 * this function is the C API for calling calculate_ddg_matrix from python
 * it is almost identical to the sym_capi function but calls calculate_ddg_matrix instead of calculate_sym_matrix
*/
PyObject* ddg_capi(PyObject *self, PyObject *args) {
    PyObject *python_points_list;
    int N;
    int d;
    double **data_points;
    double **ddg_matrix;
    PyObject *py_return_matrix;

    /* parse the python object */
    if (!PyArg_ParseTuple(args, "O", &python_points_list)) {
        return NULL; /* error is raised by parsing function */
    }
    /* convert python object to C matrix using the helper method */
    data_points = py_to_c_matrix(python_points_list, &N, &d);
    if (data_points == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    /* calculate the degree matrix */
    ddg_matrix = calculate_ddg_matrix(data_points, N, d);
    free_matrix(data_points, N); /* we dont need the data points anymore */
    if (ddg_matrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; 
    }
    /* convert the C matrix to a python object */
    py_return_matrix = c_to_py_matrix(ddg_matrix, N, N);
    free_matrix(ddg_matrix, N); /* we dont need the degree matrix anymore */
    if (py_return_matrix == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    return py_return_matrix; 
}
/*
 * ========================================NORM_CAPI=================================================
 * this function is the C API for calling calculate_norm_matrix from python
 * it is almost identical to the sym_capi function but calls calculate_norm_matrix instead of calculate_sym_matrix
*/
PyObject* norm_capi(PyObject *self, PyObject *args) {
    PyObject *python_points_list;
    int N;
    int d;
    double **data_points;
    double **norm_matrix;
    PyObject *py_return_matrix;

    /* parse the python object */
    if (!PyArg_ParseTuple(args, "O", &python_points_list)) {
        return NULL; /* error is raised by parsing function */
    }
    /* convert python object to C matrix using the helper method */
    data_points = py_to_c_matrix(python_points_list, &N, &d);
    if (data_points == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    /* calculate the degree matrix */
    norm_matrix = calculate_norm_matrix(data_points, N, d);
    free_matrix(data_points, N); /* we dont need the data points anymore */
    if (norm_matrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; 
    }
    /* convert the C matrix to a python object */
    py_return_matrix = c_to_py_matrix(norm_matrix, N, N);
    free_matrix(norm_matrix, N); /* we dont need the normalized matrix anymore */
    if (py_return_matrix == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    return py_return_matrix; 
}
/* 
 * ========================================SYMNMF_CAPI=============================================
 * this function is the C API for calling symnmf from python
 * in additio to the data points received from python, it also receives the initial randomized matrix H
*/
PyObject* symnmf_capi(PyObject *self, PyObject *args) {
    PyObject* python_points_list;
    PyObject* python_init_H;
    int N;
    int d;
    int N_H; /* place holder for py_to_c output */
    int k; 
    double **data_points;
    double **norm_matrix;
    double **init_H;
    double **optimized_H;
    PyObject* py_return_matrix;

    /* parse the python object */
    if (!PyArg_ParseTuple(args, "OO", &python_points_list, &python_init_H)) {
        return NULL; /* error is raised by parsing function */
    }
    /* get matrix W by converting python_points_list to C and calling C methods */
    data_points = py_to_c_matrix(python_points_list, &N, &d);
    if (data_points == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    norm_matrix = calculate_norm_matrix(data_points, N, d);
    free_matrix(data_points, N); /* we dont need the data points anymore */
    if (norm_matrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; 
    }

    /* convert initial H matrix to C */
    init_H = py_to_c_matrix(python_init_H, &N_H, &k);
    if (init_H == NULL) {
        free_matrix(norm_matrix, N);
        return NULL;
    }

    /* use optimize_h to execute the algorithm */
    optimized_H = optimize_h(norm_matrix, init_H, N, k);
    free_matrix(norm_matrix, N);
    free_matrix(init_H, N_H); /* no problem with dimensions even though H is Nxk since free_matrix works by rows */
    if (optimized_H == NULL) {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }

    /* convert H back to python and return it */
    py_return_matrix = c_to_py_matrix(optimized_H, N, k);
    free_matrix(optimized_H, N);
    if (py_return_matrix == NULL) {
        return NULL; /* error is raised by parsing function */
    }
    return py_return_matrix; 
}
