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
