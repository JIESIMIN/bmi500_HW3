import math
import numbers
import random
import time


# Guard: v must be a list/tuple of finite real numbers
def check_vector(v, name):
    if not isinstance(v, (list, tuple)):
        raise TypeError(f"{name} must be a list or tuple, got {type(v).__name__}")
    for i, val in enumerate(v):
        if not isinstance(val, numbers.Real):
            raise TypeError(f"{name}[{i}] must be a real number, got {type(val).__name__}")
        if not math.isfinite(val):
            raise ValueError(f"{name}[{i}] must be finite, got {val}")


# create a function to compute the dot product of two vectors using a for loop
# add comments for the selected function
def dot_product(a, b, check=True):
    # Validate inputs; matvec_multiply passes check=False after validating once
    if check:
        check_vector(a, "a")
        check_vector(b, "b")
        if len(a) != len(b):
            raise ValueError(f"vectors must have the same length, got {len(a)} and {len(b)}")

    # Sum of element-wise products
    result = 0
    for i in range(len(a)):
        result += a[i] * b[i]

    # Finite floats can still overflow to inf
    if isinstance(result, float) and not math.isfinite(result):
        raise OverflowError("dot product overflowed")
    return result


# create a function to compute the matrix-vector product using the dot_product function
# add comments for the selected function
def matvec_multiply(A, x):
    # Validate A (list of rows) and x once, so dot_product can skip re-checking x per row
    if not isinstance(A, (list, tuple)) or len(A) == 0:
        raise ValueError("A must be a non-empty list of rows")
    check_vector(x, "x")
    for i, row in enumerate(A):
        check_vector(row, f"A[{i}]")
        if len(row) != len(x):
            raise ValueError(f"A[{i}] has {len(row)} columns but x has length {len(x)}")

    # Each entry of the result is the dot product of one row of A with x
    return [dot_product(row, x, check=False) for row in A]


# create a main function to test the matrix-vector product function using randomly generated data of size 1000x1000
# add comments for the selected function
def main():
    # Random 1000x1000 matrix and length-1000 vector
    n = 1000
    A = [[random.random() for _ in range(n)] for _ in range(n)]
    x = [random.random() for _ in range(n)]

    # Time the product and show a sample of the result
    start = time.time()
    y = matvec_multiply(A, x)
    print(f"Output length: {len(y)}, time: {time.time() - start:.3f}s")
    print(f"First 5 entries: {y[:5]}")


if __name__ == "__main__":
    main()
