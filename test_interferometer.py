"""Tests for the Interferometer class and the unified OpEx/PolyOpEx pipeline.

Includes:
  - Mach-Zehnder interferometer with free-fall Hamiltonian: known phase = k*g*T^2
  - Mach-Zehnder with harmonic trap: perturbative expansion in ω
  - Symbolic perturbation parameter ε with collect_pert_orders()
  - PolyOpEx-based Interferometer via use_poly=True
"""

import sys
sys.path.insert(0, '/home/user/InterferometryAnalytic')

import numpy as np
from scipy.linalg import expm
import sympy as sy
from sympy import Rational, I, symbols, simplify, expand, series

from Interferometry import (OpEx, Hamiltonian, Pulse, U, Interferometer,
                            C, BCHN)
from poly_operator import PolyOpEx, moyal_commutator, bchn, hbar

# --- Fock space infrastructure (from test_bch.py) ---
N_FOCK = 50; M_BLOCK = 20

a_op = np.zeros((N_FOCK, N_FOCK), dtype=complex)
for i in range(N_FOCK - 1):
    a_op[i, i + 1] = np.sqrt(i + 1)
adag_op = a_op.T.copy()
x_mat = (a_op + adag_op) / np.sqrt(2)
p_mat = -1j * (a_op - adag_op) / np.sqrt(2)
I_mat = np.eye(N_FOCK, dtype=complex)


def opex_to_matrix(opex):
    """Convert OpEx to Fock-space matrix (hbar=1)."""
    a = complex(sy.simplify(opex.a).subs(hbar, 1))
    b = complex(sy.simplify(opex.b).subs(hbar, 1))
    c = complex(sy.simplify(opex.c).subs(hbar, 1))
    d = complex(sy.simplify(opex.d).subs(hbar, 1))
    e = complex(sy.simplify(opex.e).subs(hbar, 1))
    f = complex(sy.simplify(opex.f).subs(hbar, 1))
    return (a * p_mat @ p_mat + b * p_mat +
            c * (x_mat @ p_mat + p_mat @ x_mat) +
            d * x_mat + e * I_mat + f * x_mat @ x_mat)


def interior_norm(A, B):
    diff = np.linalg.norm((A - B)[:M_BLOCK, :M_BLOCK])
    ref = max(np.linalg.norm(B[:M_BLOCK, :M_BLOCK]), 1e-15)
    return diff / ref


# ===== Tests ================================================================

def test_mach_zehnder_free_fall():
    """Mach-Zehnder interferometer in free fall: phase = k*g*T^2.

    Sequence:
      Upper: Pulse(k) → evolve(T) → Pulse(-k) → evolve(T)
      Lower: evolve(T) → Pulse(k) → evolve(T) → Pulse(-k)

    For H = p²/(2m) + m*g*x, the well-known result is φ = k*g*T².
    """
    print("=== Mach-Zehnder free fall (phase = k*g*T²) ===")

    m, g_acc, k_val, T_val = symbols('m g k T', positive=True)
    hbar_sym = symbols('hbar')

    H = Hamiltonian([Rational(1, 2) / m, 0, 0, m * g_acc, 0, 0])

    upper = [Pulse(k_val), U(H, T_val), Pulse(-k_val), U(H, T_val)]
    lower = [U(H, T_val), Pulse(k_val), U(H, T_val), Pulse(-k_val)]

    interf = Interferometer(upper, lower, BCHOrder=4)
    res_dic, res_opex = interf.overlap()

    # The overlap operator should be exp(i*phase) times identity
    # phase = coefficient of the identity operator, divided by i
    # In the BCH result, the constant term (e) gives the scalar part
    phase_raw = simplify(res_dic['const'])

    # For a Mach-Zehnder with uniform gravity, the result should be:
    # const = i * k * g * T^2  (with appropriate sign convention)
    # We check that const / (I * k * g * T^2) simplifies to ±1
    # and all other quadratic terms vanish at leading order

    # Extract the leading gravity phase by substituting simple values
    phase_numeric = phase_raw.subs([(m, 1), (hbar_sym, 1)])
    expected = I * k_val * g_acc * T_val**2
    ratio = simplify(phase_numeric / expected)

    ok = ratio == 1 or ratio == -1
    print(f"  Phase = {'+' if ratio == 1 else '-'}k*g*T²: [{'PASS' if ok else 'FAIL'}]")
    if not ok:
        print(f"    phase = {phase_numeric}")
        print(f"    ratio = {ratio}")

    # Also check that p, x, p², x², px+xp terms vanish
    operator_terms_ok = True
    for name in ['p2', 'p', 'px_xp', 'x', 'x2']:
        val = simplify(res_dic[name].subs([(m, 1), (hbar_sym, 1)]))
        if val != 0:
            print(f"    {name} = {val} (should be 0) [FAIL]")
            operator_terms_ok = False

    ok2 = operator_terms_ok
    print(f"  All operator terms vanish: [{'PASS' if ok2 else 'FAIL'}]")

    return ok and ok2


def test_mach_zehnder_numerical():
    """Numerical verification of Mach-Zehnder overlap against matrix exponentiation."""
    print("\n=== Mach-Zehnder numerical verification ===")

    # Use small numeric parameters so BCH converges well
    m_val, g_val, k_num, T_num = Rational(1), Rational(1, 100), Rational(1, 10), Rational(1, 50)
    hbar_sym = symbols('hbar')

    H = Hamiltonian([Rational(1, 2) / m_val, 0, 0, m_val * g_val, 0, 0])

    upper = [Pulse(k_num), U(H, T_num), Pulse(-k_num), U(H, T_num)]
    lower = [U(H, T_num), Pulse(k_num), U(H, T_num), Pulse(-k_num)]

    interf = Interferometer(upper, lower, BCHOrder=8)
    _, res_opex = interf.overlap()

    # Build the overlap matrix from the analytic result
    Z_mat = opex_to_matrix(res_opex)
    expZ = expm(Z_mat)

    # Build the overlap matrix directly from the sequence
    # Upper: exp(-iH*T/hbar) produces unitary evolution
    # Each element in the sequence is exp(-i*H*t/hbar)
    H_mat = opex_to_matrix(H)

    def pulse_mat(k):
        """exp(i*k*x) matrix."""
        return expm(1j * float(k) * x_mat)

    def evolve_mat(t):
        """exp(-i*H*t/hbar) matrix with hbar=1."""
        return expm(-1j * H_mat * float(t))

    # Upper arm: Pulse(k) → U(T) → Pulse(-k) → U(T)
    U_upper = evolve_mat(T_num) @ pulse_mat(-k_num) @ evolve_mat(T_num) @ pulse_mat(k_num)

    # Lower arm: U(T) → Pulse(k) → U(T) → Pulse(-k)
    U_lower = pulse_mat(-k_num) @ evolve_mat(T_num) @ pulse_mat(k_num) @ evolve_mat(T_num)

    # Overlap = U_lower^† @ U_upper
    overlap_mat = U_lower.conj().T @ U_upper

    err = interior_norm(expZ, overlap_mat)
    ok = err < 1e-8
    print(f"  Analytic vs matrix overlap: error = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")

    return ok


def test_interferometer_poly_mode():
    """Test Interferometer with use_poly=True matches OpEx result for quadratic H."""
    print("\n=== Interferometer poly mode vs OpEx mode ===")

    m_val = Rational(1)
    g_val = Rational(1, 10)
    k_num = Rational(2)
    T_num = Rational(1, 10)

    H = Hamiltonian([Rational(1, 2) / m_val, 0, 0, m_val * g_val, 0, 0])

    upper = [Pulse(k_num), U(H, T_num), Pulse(-k_num), U(H, T_num)]
    lower = [U(H, T_num), Pulse(k_num), U(H, T_num), Pulse(-k_num)]

    # OpEx mode (original)
    interf_opex = Interferometer(upper, lower, BCHOrder=4)
    dic_opex, _ = interf_opex.overlap()

    # PolyOpEx mode
    interf_poly = Interferometer(upper, lower, BCHOrder=4, use_poly=True)
    dic_poly, _ = interf_poly.overlap()

    # Compare all terms
    all_pass = True
    for name in ['p2', 'p', 'px_xp', 'x', 'const', 'x2']:
        v1 = dic_opex[name]
        v2 = dic_poly[name]
        diff = abs(complex(simplify(v1 - v2).subs(hbar, 1)))
        ok = diff < 1e-12
        if not ok:
            print(f"  {name}: diff = {diff:.2e}  [FAIL]")
            all_pass = False

    print(f"  All terms match (BCH order 4): [{'PASS' if all_pass else 'FAIL'}]")
    return all_pass


def test_symbolic_epsilon():
    """Test symbolic perturbation parameter ε with collect_pert_orders()."""
    print("\n=== Symbolic perturbation parameter ε ===")
    all_pass = True

    eps = symbols('epsilon')

    # H0 = p²/2 + x²/2 (quadratic, order 0)
    # V = x³ (cubic perturbation, order 1)
    H0 = PolyOpEx({(0, 2): Rational(1, 2), (2, 0): Rational(1, 2)},
                  pert_order=0, max_pert_order=2, pert_symbol=eps)
    V = PolyOpEx({(3, 0): Rational(1, 10)},
                 pert_order=1, max_pert_order=2, pert_symbol=eps)

    # Build X = -i*t*(H0 + V)/hbar with t small
    t_val = Rational(1, 50)
    X_H0 = H0 * (-I * t_val / hbar)
    X_V = V * (-I * t_val / hbar)
    X = X_H0 + X_V
    Y = X_H0 + X_V  # two identical time steps

    Z = bchn(X, Y, 4)

    # collect_pert_orders should give {0: ..., 1: ..., 2: ...}
    orders = Z.collect_pert_orders()

    ok1 = set(orders.keys()) == {0, 1, 2} or set(orders.keys()) == {0, 1}
    print(f"  Orders present: {sorted(orders.keys())}  [{'PASS' if ok1 else 'FAIL'}]")
    all_pass = all_pass and ok1

    # Order-0 terms should NOT contain ε
    if 0 in orders:
        coeffs_0 = orders[0].coeffs
        has_eps = any(eps in sy.sympify(c).free_symbols for c in coeffs_0.values())
        ok2 = not has_eps
        print(f"  O(1) terms free of ε: [{'PASS' if ok2 else 'FAIL'}]")
        all_pass = all_pass and ok2

    # Order-1 terms should contain ε^1 (after collect_pert_orders multiplies by ε)
    if 1 in orders:
        coeffs_1 = orders[1].coeffs
        # Each coefficient should be divisible by ε
        all_have_eps = all(eps in sy.sympify(c).free_symbols for c in coeffs_1.values())
        ok3 = all_have_eps
        print(f"  O(ε) terms contain ε: [{'PASS' if ok3 else 'FAIL'}]")
        all_pass = all_pass and ok3

    # pert_symbol propagates through operations
    ok4 = Z.pert_symbol == eps
    print(f"  pert_symbol propagated: [{'PASS' if ok4 else 'FAIL'}]")
    all_pass = all_pass and ok4

    return all_pass


def test_harmonic_trap_perturbative():
    """Mach-Zehnder in a harmonic trap, treating ω² perturbatively.

    H = p²/(2m) + m*g*x + m*ω²*x²/2

    The trap term is perturbative. At O(ω⁰), the phase is k*g*T².
    At O(ω²), there's a correction proportional to ω²*T⁴.
    """
    print("\n=== Harmonic trap perturbative (ε = ω²) ===")

    # Use small k and T so BCH converges well across the full sequence
    m_val = Rational(1)
    g_val = Rational(1)
    k_val = Rational(1, 50)
    T_val = Rational(1, 50)

    # Quadratic H0 = p²/2m + m*g*x
    H0_coeffs = {(0, 2): Rational(1, 2) / m_val,
                 (1, 0): m_val * g_val}

    # Perturbation: V = m*ω²*x²/2 → coefficient 1/2, ω² absorbed into ε
    V_coeffs = {(2, 0): m_val * Rational(1, 2)}

    eps = symbols('epsilon')

    H0 = PolyOpEx(H0_coeffs, pert_order=0, max_pert_order=2, pert_symbol=eps)
    V = PolyOpEx(V_coeffs, pert_order=1, max_pert_order=2, pert_symbol=eps)
    full_H = H0 + V

    # Build MZ interferometer using PolyOpEx BCH
    # Pulse(k) exponent: i*k*x (from H*(-i/hbar*1) with H={x: -hbar*k})
    pulse_plus = PolyOpEx({(1, 0): I * k_val}, pert_order=0,
                          max_pert_order=2, pert_symbol=eps)
    pulse_minus = PolyOpEx({(1, 0): -I * k_val}, pert_order=0,
                           max_pert_order=2, pert_symbol=eps)

    evolve_exp = full_H * (-I * T_val / hbar)

    # Upper arm
    Oup = PolyOpEx.zero(max_pert_order=2, pert_symbol=eps)
    for arg in [pulse_plus, evolve_exp, pulse_minus, evolve_exp]:
        Oup = bchn(Oup, arg, 8).simplify()

    # Lower arm
    Olow = PolyOpEx.zero(max_pert_order=2, pert_symbol=eps)
    for arg in [evolve_exp, pulse_plus, evolve_exp, pulse_minus]:
        Olow = bchn(arg, Olow, 8).simplify()

    Z = bchn(Olow, Oup, 8).simplify()
    orders = Z.collect_pert_orders()

    all_pass = True

    # Check O(1) phase ≈ k*g*T² (in the constant term)
    # With BCH order 8 and small parameters, this should be very close
    if 0 in orders:
        const_0 = orders[0].coeffs.get((0, 0), 0)
        const_0_val = simplify(const_0.subs(hbar, 1))
        expected_phase = I * k_val * g_val * T_val**2
        # Use numerical comparison since BCH truncation gives small corrections
        diff = abs(complex(const_0_val - expected_phase))
        magnitude = abs(complex(expected_phase))
        rel_err = diff / magnitude if magnitude > 0 else diff
        # 9 composed BCH operations at order 8 accumulate truncation error
        ok = rel_err < 1e-3
        print(f"  O(1) phase ≈ k*g*T²: rel_err = {rel_err:.2e}  [{'PASS' if ok else 'FAIL'}]")
        all_pass = all_pass and ok
    else:
        print(f"  O(1) terms missing  [FAIL]")
        all_pass = False

    # Check that O(ε) correction exists
    if 1 in orders:
        ok2 = len(orders[1].coeffs) > 0
        print(f"  O(ε) correction present ({len(orders[1].coeffs)} terms): [{'PASS' if ok2 else 'FAIL'}]")
        all_pass = all_pass and ok2
    else:
        print(f"  O(ε) correction missing  [FAIL]")
        all_pass = False

    return all_pass


if __name__ == '__main__':
    print("=" * 60)
    print("Interferometer Integration Test Suite")
    print("=" * 60)

    results = []
    results.append(("MZ free fall (analytic)", test_mach_zehnder_free_fall()))
    results.append(("MZ free fall (numerical)", test_mach_zehnder_numerical()))
    results.append(("Interferometer poly mode", test_interferometer_poly_mode()))
    results.append(("Symbolic epsilon", test_symbolic_epsilon()))
    results.append(("Harmonic trap perturbative", test_harmonic_trap_perturbative()))

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
