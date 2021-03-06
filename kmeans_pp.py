import numpy as np
import pandas as pd
import sys
from mykmeanssp import fit

def validate(condition):
    if not condition:
        print('Invalid Input!')
        exit(1)


def check_if_float(num):
    try:
        float(num)
    except ValueError:
        print('Invalid Input!')
        exit(1)


def read_file_to_df(file):
    if file.split(".")[-1] == "txt":
        return pd.read_csv(file, sep=",", header=None)
    else:
        return pd.read_csv(file)


def initialize_centroids(vectors, k):
    mu = np.zeros((k, vectors.shape[1]))
    mu_index = np.zeros(k)
    np.random.seed(0)
    N = vectors.shape[0]
    # print(vectors.shape[1])
    mu_index[0] = np.random.choice(N)
    mu[0] = vectors[(int(mu_index[0]))]
    for i in range(1, k):
        D = np.apply_along_axis(lambda vector: min(np.linalg.norm(vector-mu[:i], axis=1)**2), axis=1, arr=vectors)
        p = D/np.sum(D)
        mu_index[i] = np.random.choice(N, p=p)
        mu[i] = vectors[int(mu_index[i])]
    return mu_index.astype(int), mu.astype(float)


args = sys.argv[1:]
validate(4 <= len(args) <= 5)
k, max_iter, epsilon, input1_filename, input2_filename = args[0], args[1] if len(args) == 5 else 300, args[-3], args[-2], args[-1]
validate(args[0].isdigit())
check_if_float(args[-3])
validate(float(args[-3]) >= 0)
if len(args) == 5:
    validate(max_iter.isdigit())
    max_iter = int(max_iter)
k = int(k)
epsilon = float(epsilon)
df1 = read_file_to_df(input1_filename)
df2 = read_file_to_df(input2_filename)
vectors = pd.merge(df1, df2, on=0).sort_values(by=0).iloc[:, 1:].to_numpy()
centroids_index, centroids = initialize_centroids(vectors, k)
final_centroids = fit(k, vectors.shape[1], vectors.shape[0], max_iter, epsilon, centroids.tolist(), vectors.tolist())


for i in range(len(centroids_index)):
    if i < len(centroids_index) - 1:
        print(str(centroids_index[i]) + ",", end="")
    else:
        print(str(centroids_index[i]))
for i in range(len(final_centroids)):
    for j in range(len(final_centroids[0])):
            if j < len(final_centroids[0]) - 1:
                print(str("{0:.4f}".format(final_centroids[i][j])) + ",", end="")
            else:
                print(str("{0:.4f}".format(final_centroids[i][j])))

