"""
Numerical verification of the BCH implementation in Interferometry.py.

Strategy: Represent the quadratic operator algebra {p², p, (px+xp), x, 1, x²}
as matrices using a truncated Fock space (quantum harmonic oscillator basis).
Then verify that exp(BCHN(X,Y)) ≈ exp(X)·exp(Y) by comparing matrix exponentials.

We compare only the interior block of the matrices to avoid Fock space truncation artifacts.
"""

import numpy as np
from scipy.linalg import expm
import sympy as sy
import sys
sys.path.insert(0, '/home/user/InterferometryAnalytic')
from Interferometry import OpEx, BCHN, C

# --- Build matrix representation of the operator algebra ---
N = 50   # Fock space truncation
M = 20   # Interior block size for comparison (avoids edge effects)

# Creation and annihilation operators
a_op = np.zeros((N, N), dtype=complex)
for i in range(N - 1):
    a_op[i, i + 1] = np.sqrt(i + 1)
adag_op = a_op.T.copy()

# Position and momentum with hbar=1
x_mat = (a_op + adag_op) / np.sqrt(2)
p_mat = -1j * (a_op - adag_op) / np.sqrt(2)
I_mat = np.eye(N, dtype=complex)

# Basis operator matrices
p2_mat = p_mat @ p_mat
x2_mat = x_mat @ x_mat
xp_px_mat = x_mat @ p_mat + p_mat @ x_mat


def opex_to_matrix(a, b, c, d, e, f):
    return a * p2_mat + b * p_mat + c * xp_px_mat + d * x_mat + e * I_mat + f * x2_mat


def opex_to_params(opex_obj):
    """Extract complex numeric parameters from an OpEx (after substituting hbar=1)."""
    hbar = sy.symbols('hbar')
    params = []
    for attr in ['a', 'b', 'c', 'd', 'e', 'f']:
        val = getattr(opex_obj, attr)
        val = complex(sy.simplify(val).subs(hbar, 1))
        params.append(val)
    return params


def interior_norm(A, B):
    """Frobenius norm of (A-B) restricted to interior MxM block, relative to B."""
    diff = np.linalg.norm((A - B)[:M, :M])
    ref = max(np.linalg.norm(B[:M, :M]), 1e-15)
    return diff / ref


def verify_commutator():
    """Verify C() against matrix commutator for random inputs."""
    print("=== Verifying commutator C() ===")

    np.random.seed(42)
    all_pass = True
    for trial in range(5):
        xp = [round(v, 4) for v in np.random.uniform(-0.1, 0.1, 6)]
        yp = [round(v, 4) for v in np.random.uniform(-0.1, 0.1, 6)]

        # Symbolic commutator via C()
        X_sym = OpEx([sy.Rational(str(v)) for v in xp])
        Y_sym = OpEx([sy.Rational(str(v)) for v in yp])
        CXY_sym = C(X_sym, Y_sym)
        CXY_params = opex_to_params(CXY_sym)

        # Matrix commutator
        X_m = opex_to_matrix(*xp)
        Y_m = opex_to_matrix(*yp)
        CXY_direct = X_m @ Y_m - Y_m @ X_m
        CXY_code = opex_to_matrix(*CXY_params)

        err = interior_norm(CXY_direct, CXY_code)
        ok = err < 1e-10
        all_pass = all_pass and ok
        print(f"  Trial {trial+1}: relative error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")

    return all_pass


def verify_bch(X_params, Y_params, order, label="", tol=1e-6):
    """Verify that exp(BCHN(X,Y)) ≈ exp(X)·exp(Y)."""
    X_sym = OpEx([sy.Rational(str(v)) for v in X_params])
    Y_sym = OpEx([sy.Rational(str(v)) for v in Y_params])

    Z_sym = BCHN(X_sym, Y_sym, order)
    Z_params = opex_to_params(Z_sym)

    X_m = opex_to_matrix(*X_params)
    Y_m = opex_to_matrix(*Y_params)
    Z_m = opex_to_matrix(*Z_params)

    expZ = expm(Z_m)
    expX_expY = expm(X_m) @ expm(Y_m)

    err = interior_norm(expZ, expX_expY)
    ok = err < tol
    print(f"  {label}BCH order {order}: relative error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")
    return err, ok


def verify_bch_orders():
    """Verify BCH at each order with small parameters."""
    print("\n=== Verifying BCHN accuracy at each order ===")

    # Use small parameters so series converges well
    scale = 0.03
    np.random.seed(123)
    all_pass = True

    for trial in range(3):
        xp = [round(v, 6) for v in np.random.uniform(-scale, scale, 6)]
        yp = [round(v, 6) for v in np.random.uniform(-scale, scale, 6)]
        print(f"\n  Trial {trial+1}:")
        for order in [2, 4, 6, 8]:
            # Tolerance scales with truncation error: O(scale^(order+1))
            tol = 100 * (scale ** (order + 1))
            _, ok = verify_bch(xp, yp, order, label="", tol=tol)
            all_pass = all_pass and ok

    return all_pass


def verify_bch_convergence():
    """Verify that higher BCH orders give strictly better approximations."""
    print("\n=== Verifying BCH convergence (error should decrease with order) ===")

    scale = 0.05
    np.random.seed(456)
    xp = [round(v, 6) for v in np.random.uniform(-scale, scale, 6)]
    yp = [round(v, 6) for v in np.random.uniform(-scale, scale, 6)]

    X_m = opex_to_matrix(*xp)
    Y_m = opex_to_matrix(*yp)
    target = expm(X_m) @ expm(Y_m)

    prev_err = float('inf')
    all_pass = True
    errors = []

    for order in range(1, 9):
        X_sym = OpEx([sy.Rational(str(v)) for v in xp])
        Y_sym = OpEx([sy.Rational(str(v)) for v in yp])
        Z_sym = BCHN(X_sym, Y_sym, order)
        Z_params = opex_to_params(Z_sym)
        Z_m = opex_to_matrix(*Z_params)

        expZ = expm(Z_m)
        err = interior_norm(expZ, target)
        errors.append(err)

        converging = err <= prev_err * 1.1  # allow small noise
        status = "PASS" if converging else "FAIL"
        if not converging:
            all_pass = False
        print(f"  Order {order}: error = {err:.2e}  [{status}]")
        prev_err = err

    # Also check that order 8 is much better than order 1
    if errors[-1] > errors[0] * 0.5:
        print("  WARNING: Order 8 not significantly better than order 1")
        all_pass = False

    return all_pass


def verify_known_identities():
    """Verify known algebraic identities as sanity checks."""
    print("\n=== Verifying known algebraic identities ===")
    all_pass = True

    # BCH(0, Y) = Y
    zero = OpEx([0, 0, 0, 0, 0, 0])
    Y = OpEx([sy.Rational('1/7'), sy.Rational('1/3'), 0, sy.Rational('1/5'), sy.Rational('1/11'), 0])
    Z = BCHN(zero, Y, 8)
    hbar = sy.symbols('hbar')
    err = sum(abs(complex(sy.simplify(getattr(Z, a) - getattr(Y, a)).subs(hbar, 1)))
              for a in 'abcdef')
    ok = err < 1e-15
    all_pass = all_pass and ok
    print(f"  BCH(0, Y) = Y: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")

    # BCH(X, 0) = X
    X = OpEx([sy.Rational('1/7'), sy.Rational('1/3'), 0, sy.Rational('1/5'), sy.Rational('1/11'), 0])
    Z = BCHN(X, zero, 8)
    err = sum(abs(complex(sy.simplify(getattr(Z, a) - getattr(X, a)).subs(hbar, 1)))
              for a in 'abcdef')
    ok = err < 1e-15
    all_pass = all_pass and ok
    print(f"  BCH(X, 0) = X: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")

    # Grade-4 identity: [Y,[X,[X,Y]]] = [X,[Y,[X,Y]]]
    X = OpEx([sy.Rational('1/3'), sy.Rational('2/5'), sy.Rational('1/7'),
              sy.Rational('3/11'), sy.Rational('1/13'), sy.Rational('2/9')])
    Y = OpEx([sy.Rational('2/7'), sy.Rational('1/9'), sy.Rational('3/11'),
              sy.Rational('1/5'), sy.Rational('2/13'), sy.Rational('1/3')])
    E2 = C(X, Y)
    lhs = C(Y, C(X, E2))  # [Y,[X,[X,Y]]]
    rhs = C(X, C(Y, E2))  # [X,[Y,[X,Y]]]
    err = sum(abs(complex(sy.simplify(getattr(lhs, a) - getattr(rhs, a)).subs(hbar, 1)))
              for a in 'abcdef')
    ok = err < 1e-10
    all_pass = all_pass and ok
    print(f"  Grade-4 identity [Y,[X,[X,Y]]]=[X,[Y,[X,Y]]]: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")

    return all_pass


if __name__ == '__main__':
    print("=" * 60)
    print("BCH Implementation Verification Test Suite")
    print("=" * 60)

    results = []
    results.append(("Known identities", verify_known_identities()))
    results.append(("Commutator C()", verify_commutator()))
    results.append(("BCH orders", verify_bch_orders()))
    results.append(("BCH convergence", verify_bch_convergence()))

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    all_ok = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        all_ok = all_ok and passed

    print(f"\nOverall: {'ALL TESTS PASSED' if all_ok else 'SOME TESTS FAILED'}")
    sys.exit(0 if all_ok else 1)
