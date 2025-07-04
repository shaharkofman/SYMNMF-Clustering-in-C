#include <Python.h>

double** calculate_sym_matrix(double **data_points, int N, int d);
double** calculate_ddg_matrix(double **data_points, int N, int d);
double** calculate_norm_matrix(double **data_points, int N, int d);
double** optimize_h(double** W, double** init_H, int N, int k);
void free_matrix(double **matrix, int N);
static PyObject* sym_capi(PyObject *self, PyObject *args);
static PyObject* ddg_capi(PyObject *self, PyObject *args);
static PyObject* norm_capi(PyObject *self, PyObject *args);
static PyObject* symnmf_capi(PyObject *self, PyObject *args);
