import math
import random

from matvec_multiply import dot_product, matvec_multiply


# Check that calling fn(*args) raises the given exception type
def raises(exc, fn, *args):
    try:
        fn(*args)
    except exc:
        return True
    return False


# generate tests for matrix-vector-product function in matvec_multiply.py
def test_matvec_known_values():
    assert matvec_multiply([[1, 2], [3, 4]], [5, 6]) == [17, 39]


def test_matvec_identity():
    # Identity matrix returns the vector unchanged
    x = [3, -1, 7]
    I = [[int(i == j) for j in range(3)] for i in range(3)]
    assert matvec_multiply(I, x) == x


def test_matvec_non_square():
    # 2x3 matrix times length-3 vector gives length-2 vector
    assert matvec_multiply([[1, 0, 2], [0, 1, -1]], [1, 2, 3]) == [7, -1]


def test_matvec_zero_vector():
    assert matvec_multiply([[1, 2], [3, 4]], [0, 0]) == [0, 0]


def test_matvec_dimension_mismatch():
    # Number of columns must equal vector length
    assert raises(ValueError, matvec_multiply, [[1, 2], [3, 4]], [1, 2, 3])


def test_matvec_invalid_inputs():
    # Ragged rows, empty/None/1-D matrix, and a column-vector x are rejected
    assert raises(ValueError, matvec_multiply, [[1, 2], [3]], [1, 1])
    assert raises(ValueError, matvec_multiply, [], [1, 2])
    assert raises(ValueError, matvec_multiply, None, [1, 2])
    assert raises(TypeError, matvec_multiply, [1, 2], [1, 2])
    assert raises(TypeError, matvec_multiply, [[1, 2]], [[1], [2]])


def test_matvec_random_linearity():
    # A(x + z) should equal Ax + Az
    n = 50
    A = [[random.random() for _ in range(n)] for _ in range(n)]
    x = [random.random() for _ in range(n)]
    z = [random.random() for _ in range(n)]
    lhs = matvec_multiply(A, [xi + zi for xi, zi in zip(x, z)])
    rhs = [a + b for a, b in zip(matvec_multiply(A, x), matvec_multiply(A, z))]
    assert all(math.isclose(l, r) for l, r in zip(lhs, rhs))


# generate tests for dot-product function in matvec_multiply.py
def test_dot_known_values():
    assert dot_product([1, 2, 3], [4, 5, 6]) == 32


def test_dot_orthogonal():
    assert dot_product([1, 0], [0, 1]) == 0


def test_dot_negative_and_float():
    assert math.isclose(dot_product([-1.5, 2.0], [2.0, -0.5]), -4.0)


def test_dot_commutative():
    a = [random.random() for _ in range(100)]
    b = [random.random() for _ in range(100)]
    assert math.isclose(dot_product(a, b), dot_product(b, a))


def test_dot_empty():
    assert dot_product([], []) == 0


def test_dot_length_mismatch():
    assert raises(ValueError, dot_product, [1, 2], [1, 2, 3])


def test_dot_invalid_inputs():
    # Non-list inputs, non-numeric entries, NaN/inf, and overflow are rejected
    assert raises(TypeError, dot_product, None, [1])
    assert raises(TypeError, dot_product, "ab", [2, 3])
    assert raises(TypeError, dot_product, ["1", 2], [2, 3])
    assert raises(ValueError, dot_product, [float("nan"), 1], [1, 1])
    assert raises(ValueError, dot_product, [float("inf")], [1])
    assert raises(OverflowError, dot_product, [1e200], [1e200])


# Run all tests without pytest
if __name__ == "__main__":
    tests = [f for name, f in dict(globals()).items() if name.startswith("test_")]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"All {len(tests)} tests passed")
