"""Tests for the generalized polynomial operator algebra."""

import sys
sys.path.insert(0, '/home/user/InterferometryAnalytic')

import numpy as np
from scipy.linalg import expm
import sympy as sy
from sympy import Rational, I

from Interferometry import OpEx, C, BCHN
from poly_operator import (PolyOpEx, moyal_commutator, bchn, BCH_WORDS,
                           interaction_picture, hbar)

# --- Fock space matrix representation (reused from test_bch.py) ---
N = 50; M_BLOCK = 20

a_op = np.zeros((N, N), dtype=complex)
for i in range(N - 1):
    a_op[i, i + 1] = np.sqrt(i + 1)
adag_op = a_op.T.copy()

x_mat = (a_op + adag_op) / np.sqrt(2)
p_mat = -1j * (a_op - adag_op) / np.sqrt(2)
I_mat = np.eye(N, dtype=complex)
p2_mat = p_mat @ p_mat
x2_mat = x_mat @ x_mat
xp_px_mat = x_mat @ p_mat + p_mat @ x_mat


def poly_to_matrix(poly):
    """Convert a PolyOpEx to a Fock-space matrix (hbar=1, numeric coefficients)."""
    mat = np.zeros((N, N), dtype=complex)
    for (a, b), c in poly.coeffs.items():
        coeff = complex(sy.simplify(c).subs(hbar, 1))
        # Build Sym(x^a p^b) matrix: Weyl-ordered = average over orderings
        # For small a, b this is straightforward
        mat += coeff * _weyl_ordered_matrix(a, b)
    return mat


def _weyl_ordered_matrix(a, b):
    """Matrix for the Weyl-ordered operator Sym(x^a p^b).

    Uses the identity: the Weyl-ordered matrix equals the matrix obtained
    by evaluating the Weyl symbol x^a p^b at the operator matrices,
    using symmetric (Weyl) ordering.

    For practical computation, we use the fact that for the Weyl symbol x^a p^b,
    the corresponding operator is (1/2^n) Σ_perm x^a_1 p^b_1 x^a_2 p^b_2 ...
    averaged over all orderings.  For small a+b, we compute this directly.

    Actually, a cleaner approach: use the relationship between Weyl and
    standard ordering via the Bopp shift.
    """
    # For monomials up to reasonable degree, compute via explicit symmetrization
    # of the operator product.  This is exact in the Fock space.
    from itertools import permutations

    # Build the word: a copies of 'x' and b copies of 'p'
    word = ['x'] * a + ['p'] * b
    n = len(word)
    if n == 0:
        return I_mat.copy()

    op_map = {'x': x_mat, 'p': p_mat}

    # Average over all distinct permutations
    seen = set()
    total = np.zeros((N, N), dtype=complex)
    count = 0
    for perm in permutations(word):
        if perm in seen:
            continue
        seen.add(perm)
        mat = I_mat.copy()
        for letter in perm:
            mat = mat @ op_map[letter]
        total += mat
        count += 1

    return total / count


def interior_norm(A, B):
    """Relative Frobenius norm of (A-B) on interior block."""
    diff = np.linalg.norm((A - B)[:M_BLOCK, :M_BLOCK])
    ref = max(np.linalg.norm(B[:M_BLOCK, :M_BLOCK]), 1e-15)
    return diff / ref


# ===== Tests ================================================================

def test_conversion_roundtrip():
    """Test OpEx -> PolyOpEx -> OpEx round-trip."""
    print("=== Conversion round-trip ===")
    all_pass = True

    for trial in range(5):
        params = [Rational(i + trial, 7 + trial) for i in range(6)]
        original = OpEx(params)
        poly = PolyOpEx.from_opex(original)
        restored = poly.to_opex()

        ok = all(sy.simplify(getattr(original, a) - getattr(restored, a)) == 0
                 for a in 'abcdef')
        all_pass = all_pass and ok
        print(f"  Trial {trial+1}: [{'PASS' if ok else 'FAIL'}]")

    return all_pass


def test_moyal_vs_C():
    """Test that moyal_commutator matches C() for quadratic operators."""
    print("\n=== Moyal bracket vs C() for quadratic operators ===")
    all_pass = True
    np.random.seed(42)

    for trial in range(10):
        xp = [Rational(int(v * 100), 100)
              for v in np.random.uniform(-0.5, 0.5, 6)]
        yp = [Rational(int(v * 100), 100)
              for v in np.random.uniform(-0.5, 0.5, 6)]

        X_opex = OpEx(xp)
        Y_opex = OpEx(yp)

        # C() result
        C_result = C(X_opex, Y_opex)

        # Moyal bracket result
        X_poly = PolyOpEx.from_opex(X_opex)
        Y_poly = PolyOpEx.from_opex(Y_opex)
        M_result = moyal_commutator(X_poly, Y_poly)
        M_as_opex = M_result.to_opex()

        # Compare
        err = sum(abs(complex(sy.simplify(
            getattr(C_result, a) - getattr(M_as_opex, a)).subs(hbar, 1)))
            for a in 'abcdef')
        ok = err < 1e-12
        all_pass = all_pass and ok
        if not ok:
            print(f"  Trial {trial+1}: error = {err:.2e}  [FAIL]")
    if all_pass:
        print(f"  All 10 trials: [PASS]")
    return all_pass


def test_moyal_cubic():
    """Test Moyal bracket for cubic operators against matrix commutator."""
    print("\n=== Moyal bracket for cubic operators (matrix verification) ===")
    all_pass = True

    # [x^3, p^2] should give 6*i*hbar * x^2*p  (Moyal bracket)
    A = PolyOpEx({(3, 0): 1})
    B = PolyOpEx({(0, 2): 1})
    C_AB = moyal_commutator(A, B)

    expected_coeff = complex((6 * I * hbar).subs(hbar, 1))
    actual_coeff = complex(C_AB.coeffs.get((2, 1), 0).subs(hbar, 1) if (2, 1) in C_AB.coeffs else 0)
    ok1 = abs(actual_coeff - expected_coeff) < 1e-12
    print(f"  [x^3, p^2] = 6iℏ·Sym(x^2·p): [{'PASS' if ok1 else 'FAIL'}]")
    all_pass = all_pass and ok1

    # Numerical check: random cubic vs quadratic
    np.random.seed(77)
    for trial in range(3):
        # Random cubic and quadratic
        c_coeffs = {(a, b): Rational(int(v * 100), 100)
                    for (a, b), v in zip(
                        [(3,0), (2,1), (1,2), (0,3), (2,0), (1,1), (0,2), (1,0), (0,1), (0,0)],
                        np.random.uniform(-0.1, 0.1, 10))}
        d_coeffs = {(a, b): Rational(int(v * 100), 100)
                    for (a, b), v in zip(
                        [(2,0), (1,1), (0,2), (1,0), (0,1), (0,0)],
                        np.random.uniform(-0.1, 0.1, 6))}

        C_poly = PolyOpEx(c_coeffs)
        D_poly = PolyOpEx(d_coeffs)

        # Moyal bracket
        CD_moyal = moyal_commutator(C_poly, D_poly)

        # Matrix commutator
        C_mat = poly_to_matrix(C_poly)
        D_mat = poly_to_matrix(D_poly)
        CD_mat_direct = C_mat @ D_mat - D_mat @ C_mat
        CD_mat_moyal = poly_to_matrix(CD_moyal)

        err = interior_norm(CD_mat_direct, CD_mat_moyal)
        ok = err < 1e-10
        all_pass = all_pass and ok
        print(f"  Cubic trial {trial+1}: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")

    return all_pass


def test_bchn_vs_BCHN():
    """Test that the word-based bchn() matches the original BCHN() for quadratic operators."""
    print("\n=== Word-based bchn() vs original BCHN() ===")
    all_pass = True
    np.random.seed(321)

    for trial in range(3):
        xp = [Rational(int(v * 1000), 1000)
              for v in np.random.uniform(-0.05, 0.05, 6)]
        yp = [Rational(int(v * 1000), 1000)
              for v in np.random.uniform(-0.05, 0.05, 6)]

        X_opex = OpEx(xp)
        Y_opex = OpEx(yp)

        X_poly = PolyOpEx.from_opex(X_opex)
        Y_poly = PolyOpEx.from_opex(Y_opex)

        for order in [2, 4, 6, 8]:
            # Original BCHN
            Z_opex = BCHN(X_opex, Y_opex, order)

            # Word-based bchn
            Z_poly = bchn(X_poly, Y_poly, order)
            Z_poly_as_opex = Z_poly.to_opex()

            # Compare
            err = sum(abs(complex(sy.simplify(
                getattr(Z_opex, a) - getattr(Z_poly_as_opex, a)).subs(hbar, 1)))
                for a in 'abcdef')
            ok = err < 1e-10
            all_pass = all_pass and ok
            if not ok:
                print(f"  Trial {trial+1}, order {order}: error = {err:.2e}  [FAIL]")

        print(f"  Trial {trial+1}, all orders: [PASS]" if all_pass else "")

    return all_pass


def test_bchn_cubic_convergence():
    """Test BCH for cubic operators converges against matrix exponential."""
    print("\n=== BCH convergence for cubic operators ===")

    scale = 0.02
    np.random.seed(555)

    # Random cubic operators
    keys = [(3,0), (2,1), (1,2), (0,3), (2,0), (1,1), (0,2), (1,0), (0,1), (0,0)]
    xc = {k: Rational(int(v * 1000), 1000) for k, v in zip(keys, np.random.uniform(-scale, scale, len(keys)))}
    yc = {k: Rational(int(v * 1000), 1000) for k, v in zip(keys, np.random.uniform(-scale, scale, len(keys)))}

    X_poly = PolyOpEx(xc)
    Y_poly = PolyOpEx(yc)

    X_mat = poly_to_matrix(X_poly)
    Y_mat = poly_to_matrix(Y_poly)
    target = expm(X_mat) @ expm(Y_mat)

    prev_err = float('inf')
    all_pass = True

    # Cap at order 7 for cubic operators: order 8 is too slow because
    # polynomial degree grows with each commutator (algebra not closed)
    for order in range(1, 8):
        Z_poly = bchn(X_poly, Y_poly, order)
        Z_mat = poly_to_matrix(Z_poly)
        expZ = expm(Z_mat)
        err = interior_norm(expZ, target)

        converging = err <= prev_err * 1.1
        ok = converging
        all_pass = all_pass and ok
        print(f"  Order {order}: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")
        prev_err = err

    return all_pass


def test_interaction_picture():
    """Test interaction picture for a simple case."""
    print("\n=== Interaction picture ===")
    all_pass = True

    # H0 = p^2/(2m) + m*omega^2*x^2/2 (harmonic oscillator)
    # For simplicity, set m=1, omega=1: H0 = p^2/2 + x^2/2
    H0 = PolyOpEx({(0, 2): Rational(1, 2), (2, 0): Rational(1, 2)})

    # Perturbation: V = x^3
    V = PolyOpEx({(3, 0): 1})

    t = sy.Symbol('t')
    V_I = interaction_picture(H0, V, t)

    # Check numerically at t=0: V_I(0) should equal V
    V_I_at_0 = PolyOpEx(
        {k: sy.simplify(v.subs(t, 0)) for k, v in V_I.coeffs.items()})

    err = 0
    for k in set(list(V.coeffs.keys()) + list(V_I_at_0.coeffs.keys())):
        err += abs(complex(
            (V.coeffs.get(k, 0) - V_I_at_0.coeffs.get(k, 0))))
    ok = err < 1e-12
    print(f"  V_I(t=0) == V: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")
    all_pass = all_pass and ok

    # Numerical check at t=0.1: compare with matrix e^{iH0t} V e^{-iH0t}
    t_val = 0.1
    H0_mat = poly_to_matrix(H0)
    V_mat = poly_to_matrix(V)
    V_I_expected = expm(1j * H0_mat * t_val) @ V_mat @ expm(-1j * H0_mat * t_val)

    V_I_num = PolyOpEx(
        {k: v.subs([(t, t_val), (hbar, 1)]) for k, v in V_I.coeffs.items()})
    V_I_mat = poly_to_matrix(V_I_num)

    err2 = interior_norm(V_I_mat, V_I_expected)
    ok2 = err2 < 1e-8
    print(f"  V_I(t=0.1) vs matrix: error = {err2:.2e}  [{'PASS' if ok2 else 'FAIL'}]")
    all_pass = all_pass and ok2

    return all_pass


def test_pert_order_tracking():
    """Test that perturbative order is correctly tracked and truncated."""
    print("\n=== Perturbative order tracking ===")
    all_pass = True

    # Create a quadratic H0 (pert_order=0) and cubic perturbation (pert_order=1)
    H0 = PolyOpEx({(0, 2): Rational(1, 2), (2, 0): Rational(1, 2)},
                  pert_order=0, max_pert_order=3)
    V = PolyOpEx({(3, 0): Rational(1, 10)},
                 pert_order=1, max_pert_order=3)

    # Commutator [H0, V]: should have pert_order = 0 + 1 = 1
    comm = moyal_commutator(H0, V)
    ok = all(comm.get_pert_order(k) == 1 for k in comm.coeffs)
    print(f"  [H0, V] pert_order = 1: [{'PASS' if ok else 'FAIL'}]")
    all_pass = all_pass and ok

    # [V, [H0, V]]: pert_order = 1 + 1 = 2
    comm2 = moyal_commutator(V, comm)
    ok2 = all(comm2.get_pert_order(k) == 2 for k in comm2.coeffs)
    print(f"  [V, [H0, V]] pert_order = 2: [{'PASS' if ok2 else 'FAIL'}]")
    all_pass = all_pass and ok2

    # [V, [V, [H0, V]]]: pert_order = 1 + 2 = 3
    comm3 = moyal_commutator(V, comm2)
    ok3 = all(comm3.get_pert_order(k) == 3 for k in comm3.coeffs)
    print(f"  [V, [V, [H0, V]]] pert_order = 3: [{'PASS' if ok3 else 'FAIL'}]")
    all_pass = all_pass and ok3

    # [V, comm3] would be pert_order 4 -> should be truncated (max=3)
    comm4 = moyal_commutator(V, comm3)
    ok4 = len(comm4.coeffs) == 0
    print(f"  pert_order 4 truncated to 0 terms: [{'PASS' if ok4 else 'FAIL'}]")
    all_pass = all_pass and ok4

    # Addition preserves per-term orders correctly
    mixed = H0 + V
    ok5 = (mixed.get_pert_order((0, 2)) == 0 and
           mixed.get_pert_order((2, 0)) == 0 and
           mixed.get_pert_order((3, 0)) == 1)
    print(f"  H0 + V preserves per-term orders: [{'PASS' if ok5 else 'FAIL'}]")
    all_pass = all_pass and ok5

    return all_pass


def test_max_degree_truncation():
    """Test that max_degree truncation makes cubic BCH tractable."""
    print("\n=== max_degree truncation ===")
    all_pass = True

    scale = 0.02
    np.random.seed(555)

    keys = [(3,0), (2,1), (1,2), (0,3), (2,0), (1,1), (0,2), (1,0), (0,1), (0,0)]
    xc = {k: Rational(int(v * 1000), 1000)
          for k, v in zip(keys, np.random.uniform(-scale, scale, len(keys)))}
    yc = {k: Rational(int(v * 1000), 1000)
          for k, v in zip(keys, np.random.uniform(-scale, scale, len(keys)))}

    # Without truncation: order 7 gives degree-15 intermediates (slow)
    # With max_degree=6: intermediates capped, much faster, still accurate
    X_trunc = PolyOpEx(xc, max_degree=6)
    Y_trunc = PolyOpEx(yc, max_degree=6)

    import time

    # Time the truncated version at order 6 (order 8 with symbolics is too slow)
    t0 = time.time()
    Z_trunc = bchn(X_trunc, Y_trunc, 6)
    t_trunc = time.time() - t0

    # Verify against matrix exponential
    X_mat = poly_to_matrix(PolyOpEx(xc))
    Y_mat = poly_to_matrix(PolyOpEx(yc))
    target = expm(X_mat) @ expm(Y_mat)

    Z_mat = poly_to_matrix(Z_trunc)
    expZ = expm(Z_mat)
    err = interior_norm(expZ, target)

    ok1 = err < 1e-2  # degree truncation introduces controllable error
    print(f"  BCH order 6 (max_degree=6): error = {err:.2e}, time = {t_trunc:.1f}s  [{'PASS' if ok1 else 'FAIL'}]")
    all_pass = all_pass and ok1

    # Verify truncation actually limits degree
    ok2 = Z_trunc.degree <= 6
    print(f"  Result degree ≤ 6: degree={Z_trunc.degree}  [{'PASS' if ok2 else 'FAIL'}]")
    all_pass = all_pass and ok2

    return all_pass


def test_pert_truncated_bch():
    """Test BCH with perturbative truncation: H0 + εV composed with itself."""
    print("\n=== Perturbative BCH truncation ===")
    all_pass = True

    # H = p²/2 + x²/2 + ε*x³  where ε = 0.01
    eps = Rational(1, 100)

    # Create H*(-i*t/hbar) as the BCH input, split into O(1) + O(ε)
    t_val = Rational(1, 10)
    quad_coeffs = {(0, 2): -I * t_val / (2 * hbar),
                   (2, 0): -I * t_val / (2 * hbar)}
    pert_coeffs = {(3, 0): -I * t_val * eps / hbar}

    X = PolyOpEx(quad_coeffs, max_degree=8, pert_order=0, max_pert_order=2)
    V_part = PolyOpEx(pert_coeffs, max_degree=8, pert_order=1, max_pert_order=2)
    full_X = X + V_part

    # Y = same (two equal time steps)
    Y = PolyOpEx(dict(full_X.coeffs), full_X.max_degree,
                 pert_order=dict(full_X._pert_order),
                 max_pert_order=full_X.max_pert_order)

    import time
    t0 = time.time()
    Z = bchn(full_X, Y, 4)
    elapsed = time.time() - t0

    # Check: result should contain O(ε^0), O(ε^1), O(ε^2) terms
    orders_present = set(Z._pert_order.values())
    ok1 = 0 in orders_present and 1 in orders_present
    print(f"  BCH has O(1) and O(ε) terms: {orders_present}  [{'PASS' if ok1 else 'FAIL'}]")
    all_pass = all_pass and ok1

    # No O(ε³) or higher terms should be present
    ok2 = all(po <= 2 for po in Z._pert_order.values())
    print(f"  No terms above O(ε²): [{'PASS' if ok2 else 'FAIL'}]")
    all_pass = all_pass and ok2

    # Performance: should be fast with truncation
    ok3 = elapsed < 30.0  # generous bound
    print(f"  Completed in {elapsed:.1f}s (< 30s): [{'PASS' if ok3 else 'FAIL'}]")
    all_pass = all_pass and ok3

    # Numerical check: compare with matrix exp at ε=0.01
    Z_num = PolyOpEx(
        {k: complex(sy.simplify(v).subs(hbar, 1))
         for k, v in Z.coeffs.items()})
    Z_mat = poly_to_matrix(Z_num)
    expZ = expm(Z_mat)

    # Build full X matrix
    full_X_num = PolyOpEx(
        {k: complex(sy.simplify(v).subs(hbar, 1))
         for k, v in full_X.coeffs.items()})
    X_mat = poly_to_matrix(full_X_num)
    target = expm(X_mat) @ expm(X_mat)

    err = interior_norm(expZ, target)
    ok4 = err < 1e-3
    print(f"  Numerical accuracy: error = {err:.2e}  [{'PASS' if ok4 else 'FAIL'}]")
    all_pass = all_pass and ok4

    return all_pass

if __name__ == '__main__':
    print("=" * 60)
    print("Polynomial Operator Algebra Test Suite")
    print("=" * 60)

    results = []
    results.append(("Conversion round-trip", test_conversion_roundtrip()))
    results.append(("Moyal vs C()", test_moyal_vs_C()))
    results.append(("Moyal cubic", test_moyal_cubic()))
    results.append(("bchn vs BCHN", test_bchn_vs_BCHN()))
    results.append(("BCH cubic convergence", test_bchn_cubic_convergence()))
    results.append(("Interaction picture", test_interaction_picture()))
    results.append(("Pert order tracking", test_pert_order_tracking()))
    results.append(("max_degree truncation", test_max_degree_truncation()))
    results.append(("Pert truncated BCH", test_pert_truncated_bch()))

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
