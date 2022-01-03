#include <math.h>
#include <Python.h>

int dimension, k, N;


// double **read_file(char fileName[]);

// void number_of_vectors(char fileName[]);

// void find_dimension(char line[]);

// void initialize(double **vectors_list, double *mu[], double *new_mu[]);

double **initialize_2d_double_array(double *arr[], int d1, int d2);

int calc_argmin(double *mu[], double *vector);

double compute_distance(double vec1[], double vec2[]);

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[]);

double calculating_epsilon(double *mu[], double *new_mu[]);

void create_output(double *mu[], char op_filename[]);

int check_allocation(const double *p);

void free_memory(double **array, int len);

double **transform_PyObject_to_2dArray(PyObject *mat, int rows, int columns);

PyObject *transform_2dArray_to_PyObject(double **mat, int rows, int columns);

static PyObject *fit(PyObject *self, PyObject *args);

// args: k, max_iter, epsilon, N, d, array of vectors, mu, output_file
int main(int argc, char *argv[]) {
    // char *input_file = argv[argc - 2];
//    char *output_file = argv[argc - 1];
//    int max_iter, i, j, q;
//    double eps;
//    double **new_mu, **mu = argv[argc - 2]; // check how to get 2d array from python
//    double **vectors_list = argv[argc - 3];
    // number_of_vectors(input_file);
    // vectors_list = read_file(input_file);
//    k = (int) strtol(argv[1], NULL, 10);
//    max_iter = (int) strtol(argv[2], NULL, 10);
//    eps = (double) strtol(argv[3], NULL, 10);
//    N = (int) strtol(argv[4], NULL, 10);
//    dimension = (int) strtol(argv[5], NULL, 10);
    // if (argc == 5) {
    //     max_iter = (int) strtol(argv[2], NULL, 10);
    // }
    // if (k < 1 || k > N) {
    //     printf("Invalid Input!");
    //     exit(1);
    // }

//    mu = (double **) malloc(k * sizeof(double *));
//    new_mu = (double **) malloc(k * sizeof(double *));
//    initialize_2d_double_array(new_mu, k, dimension);
//    // initialize(vectors_list, mu, new_mu);
//    for (i = 0; i < max_iter; ++i) {
//        reset_clusters(vectors_list, mu, new_mu);
//        eps = calculating_epsilon(mu, new_mu);
//        for (j = 0; j < k; ++j) {
//            for (q = 0; q < dimension; ++q) {
//                mu[j][q] = new_mu[j][q];
//                new_mu[j][q] = 0;
//            }
//        }
//        if (eps < 0.000001) {
//            break;
//        }
//    }
//    create_output(mu, output_file);
//    free_memory(vectors_list, N);
//    free_memory(mu, k);
//    free_memory(new_mu, k);
//    return 0;
}


// double **read_file(char fileName[]) {
//     FILE *file = fopen(fileName, "r");
//     char buff[1024], *ptr;
//     int ch, i = 0;
//     double *vector, *place;
//     double **all_vectors;
//     if (file) {
//         all_vectors = (double **) malloc(N * sizeof(double *));
//         ch = fscanf(file, "%s", buff);
//         while ((ch != '\n') && (ch != EOF)) {
//             vector = (double *) (malloc(dimension * sizeof(double)));
//             check_allocation(vector);
//             place = vector;
//             ptr = strtok(buff, ",");
//             while (ptr != NULL) {
//                 *(vector) = strtod(ptr, NULL);
//                 ptr = strtok(NULL, ",");
//                 vector++;
//             }
//             all_vectors[i] = place;
//             i++;
//             ch = fscanf(file, "%s", buff);
//         }
//         fclose(file);
//         return all_vectors;
//     } else {
//         printf("Invalid Input!");
//         exit(1);
//     }

// }

double compute_distance(double vec1[], double vec2[]) {
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        sum += (double) pow(vec1[i] - vec2[i], 2);
    }
    return sum;
}

// void find_dimension(char line[]) {
//     int i = 0;
//     char *ptr = strtok(line, ",");
//     while (ptr != NULL) {
//         ptr = strtok(NULL, ",");
//         i++;
//     }
//     dimension = i;
// }


// void number_of_vectors(char fileName[]){
//     FILE *file = fopen(fileName, "r");
//     char buff[1024];
//     int ch, n = 0;
//     ch = fscanf(file, "%s", buff);
//     while((ch != 'n') && (ch != EOF)){
//         if (n == 0){
//             find_dimension(buff);
//         }
//         n++;
//         ch = fscanf(file, "%s", buff);
//     }
//     N = n;
//     fclose(file);
// }

// void initialize(double **vectors_list, double *mu[], double *new_mu[]) {
//     int i, j;
//     double *vec;
//     for (i = 0; i < k; ++i) {
//         mu[i] = (double *) calloc(dimension, sizeof(double));
//         new_mu[i] = (double *) calloc(dimension, sizeof(double));
//         check_allocation(mu[i]);
//         check_allocation(new_mu[i]);
//         vec = vectors_list[i];
//         for (j = 0; j < dimension; j++)
//         {
//             mu[i][j] = vec[j];
//         }
//     }
// }

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

//void create_output(double *mu[], char op_filename[]) {
//    FILE *f;
//    int i, j;
//    f = fopen(op_filename, "w");
//    for (i = 0; i < k; ++i) {
//        for (j = 0; j < dimension; ++j) {
//            fprintf(f, "%0.4f", mu[i][j]);
//            if (j != dimension - 1) {
//                fprintf(f, ",");
//            }
//        }
//        fprintf(f, "\n");
//    }
//    fclose(f);
//}


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
    double eps, curr_eps, **centroids, **data_points, **final_centroids;
    PyObject *centroids_copy, *data_points_copy;
    if (!PyArg_ParseTuple(args, "iiiidOO", &k, &dimension, &N, &max_iter, &eps, &centroids_copy, &data_points_copy)) {
        return NULL;
    }
    if (!(PyList_Check(centroids_copy)) || !(PyList_Check(data_points_copy))) {
        printf("Error in Types of Arguments");
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
        centroids_copy = transform_2dArray_to_PyObject(final_centroids, k, dimension);
        free_memory(final_centroids, k);
        free_memory(centroids, k);
        free_memory(data_points, N);
        return centroids_copy;
    }

    double **transform_PyObject_to_2dArray(PyObject *mat, int rows, int columns) {
        double **new_mat;
        PyObject *row, *column;
        int i, j;
        new_mat = initialize_2d_double_array(new_mat, rows, columns);
        for (i = 0; i < rows; ++i) {
//        printf("i = %zd\n", i);
            row = PyList_GetItem(mat, i);
//        printf("Py_List_GetItem for row succeeded\n");
            for (j = 0; j < columns; ++j) {
//            printf("j = %zd\n", j);
                column = PyList_GetItem(row, j);
//            printf("PyList_GetItem for col");
                new_mat[i][j] = PyFloat_AsDouble(column);
//            printf("located value in new_mat\n");
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