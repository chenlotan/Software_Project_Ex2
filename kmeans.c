#include <math.h>
#include <Python.h>

int dimension, k, N;


double **initialize_2d_double_array(double *arr[], int d1, int d2);

int calc_argmin(double *mu[], double *vector);

double compute_distance(double vec1[], double vec2[]);

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[]);

double calculating_epsilon(double *mu[], double *new_mu[]);

int check_allocation(const double *p);

void free_memory(double **array, int len);

double **transform_PyObject_to_2dArray(PyObject *mat, int rows, int columns);

PyObject *transform_2dArray_to_PyObject(double **mat, int rows, int columns);

static PyObject *fit(PyObject *self, PyObject *args);


double compute_distance(double vec1[], double vec2[]) {
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        sum += (double) pow(vec1[i] - vec2[i], 2);
    }
    return sum;
}


double **initialize_2d_double_array(double **arr, int d1, int d2) {
    int i;
    arr = (double **) malloc(d1 * sizeof(double *));
    for (i = 0; i < d1; ++i) {
        arr[i] = (double *) calloc(d2, sizeof(double));
    }
    return arr;
}

int calc_argmin(double *mu[], double *vector) {
    double min_val = HUGE_VAL, sum_p;
    int min_mu = 0, i;
    for (i = 0; i < k; ++i) {
        sum_p = compute_distance(mu[i], vector);
        if (sum_p < min_val) {
            min_val = sum_p;
            min_mu = i;
        }
    }
    return min_mu;
}

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[]) {
    int *count;
    double *vec;
    int t, j, r, s, q, min_mu;
    count = calloc(k, sizeof(int));

    for (t = 0; t < N; t++) {
        min_mu = calc_argmin(mu, vectors_list[t]);
        count[min_mu]++;
        vec = vectors_list[t];
        for (j = 0; j < dimension; ++j) {
            new_sum[min_mu][j] += *vec;
            vec++;
        }
    }
    for (r = 0; r < k; ++r) {
        if (count[r] == 0) {
            for (s = 0; s < dimension; ++s) {
                new_sum[r][s] = mu[r][s];
            }
        } else {
            for (q = 0; q < dimension; ++q) {
                new_sum[r][q] = new_sum[r][q] / (double) count[r];
            }
        }
    }
    free(count);
}

double calculating_epsilon(double *mu[], double *new_mu[]) {
    double eps = 0, dist;
    int i;
    for (i = 0; i < k; ++i) {
        dist = compute_distance(mu[i], new_mu[i]);
        if (eps < dist) {
            eps = dist;
        }
    }
    return eps;
}


int check_allocation(const double *p) {
    if (p == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    return 0;
}


void free_memory(double **array, int len) {
    int q;
    for (q = 0; q < len; q++) {
        free(array[q]);
    }
    free(array);
}


static PyMethodDef k_means_func[] = {
        {"fit", (PyCFunction) fit, METH_VARARGS, NULL},
        {NULL, NULL, 0,                          NULL}
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        k_means_func
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m)
        return NULL;
    return m;
}


static PyObject *fit(PyObject *self, PyObject *args) {
    int max_iter, i, j, q;
    double eps, curr_eps, **centroids, **data_points, **final_centroids = NULL;
    PyObject *centroids_copy, *data_points_copy;
    if (!PyArg_ParseTuple(args, "iiiidOO", &k, &dimension, &N, &max_iter, &eps, &centroids_copy, &data_points_copy)) {
        return NULL;
    }
    if (!(PyList_Check(centroids_copy)) || !(PyList_Check(data_points_copy))) {
        printf("Error in Types of Arguments\n");
    }
    final_centroids = initialize_2d_double_array(final_centroids, k, dimension);
    centroids = transform_PyObject_to_2dArray(centroids_copy, k, dimension);
    data_points = transform_PyObject_to_2dArray(data_points_copy, N, dimension);

    for (i = 0; i < max_iter; ++i) {
        reset_clusters(data_points, centroids, final_centroids);
        curr_eps = calculating_epsilon(centroids, final_centroids);
        for (j = 0; j < k; ++j) {
            for (q = 0; q < dimension; ++q) {
                centroids[j][q] = final_centroids[j][q];
                final_centroids[j][q] = 0;
            }
        }
        if (curr_eps < 0.000001) {
            break;
        }
    }
    centroids_copy = transform_2dArray_to_PyObject(centroids, k, dimension);
    free_memory(final_centroids, k);
    free_memory(centroids, k);
    free_memory(data_points, N);
    return centroids_copy;
    }

    double **transform_PyObject_to_2dArray(PyObject *mat, int rows, int columns) {
        double **new_mat = NULL;
        PyObject *row, *column;
        int i, j;
        new_mat = initialize_2d_double_array(new_mat, rows, columns);
        for (i = 0; i < rows; ++i) {
            row = PyList_GetItem(mat, i);
            for (j = 0; j < columns; ++j) {
                column = PyList_GetItem(row, j);
                new_mat[i][j] = PyFloat_AsDouble(column);
            }
        }
        return new_mat;
    }

    PyObject *transform_2dArray_to_PyObject(double **mat, int rows, int columns) {
        PyObject *new_mat, *row;
        int i, j;
        new_mat = PyList_New(rows);
        for (i = 0; i < rows; ++i) {
            row = PyList_New(columns);
            for (j = 0; j < columns; ++j) {
                PyList_SetItem(row, j, Py_BuildValue("f", mat[i][j]));
            }
            PyList_SetItem(new_mat, i, row);
        }
        return new_mat;
    }


    static PyMethodDef Kmeans_Methods[] = {
            {"fit",
                    (PyCFunction) fit,
                    METH_VARARGS,
                            NULL},
            {NULL, NULL, 0, NULL}
    };


    static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "mykmeanssp",
            NULL,
            -1,
            Kmeans_Methods
    };

    PyMODINIT_FUNC
    PyInit_kmeans(void) {
        PyObject *m;
        m = PyModule_Create(&moduledef);
        if (!m) {
            return NULL;
        }
        return m;
    }