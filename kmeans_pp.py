import numpy as np
import pandas as pd
import sys


def validate(condition):
    if not condition:
        print('Invalid Input!')
        exit(1)


def join_inputs(file1, file2):  #TODO
    pass


# Open the file and create array of vectors
def read_file(file_name):
    f = open(file_name, "r")
    vectors = []
    for line in f:
        if line != "\n":
            numbers = line.split(',')
            for i in range(len(numbers)):
                numbers[i] = float("%0.4f" % float(numbers[i]))
            vectors.append(numbers)
    return vectors


def initialize_centroids(vectors, k):  # TODO
    mu = []
    np.random.seed(0)
    mu.append(np.random.choice(vectors))
    for i in range(k-1):
        pass
    return mu


args = sys.argv[1:]
validate(3 <= len(args) <= 4)
k, max_iter, epsilon, input1_filename, input2_filename = args[0], args[1] if len(args) == 4 else 300, args[-3], args[-2], args[-1]
validate(args[0].isdigit())
if len(args) == 4:
    validate(max_iter.isdigit())
    max_iter = int(max_iter)
k = int(k)
input_filename = join_inputs(input1_filename, input2_filename)  #TODO
vectors = read_file(input_filename)
d = len(vectors[0])
mu = initialize_centroids(vectors, k)
