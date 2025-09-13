double** calculate_sym_matrix(double **data_points, int N, int d);
double** calculate_ddg_matrix(double **data_points, int N, int d);
double** calculate_norm_matrix(double **data_points, int N, int d);
double** optimize_h(double** W, double** init_H, int N, int k);
void free_matrix(double **matrix, int N);

/* C API methods, ifdef block in order to expose only if <Python.h> is included */
#ifdef PY_SSIZE_T_CLEAN
#include <Python.h>
PyObject* sym_capi(PyObject *self, PyObject *args);
PyObject* ddg_capi(PyObject *self, PyObject *args);
PyObject* norm_capi(PyObject *self, PyObject *args);
PyObject* symnmf_capi(PyObject *self, PyObject *args);
#endif
