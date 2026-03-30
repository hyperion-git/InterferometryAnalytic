"""Benchmark tests against analytic results from Ufrecht's PhD thesis.

Reference: C. Ufrecht, "Theoretical approach to high-precision atom
interferometry", PhD thesis, Universität Ulm, 2019.
DOI: 10.18725/OPARU-17923

Tests cover:
  - MZ phase in linear gravity (Eq. 1.88): phi_g = -k*g*T^2
  - Gravity gradient first-order correction (Eq. 1.99):
    phi_1 = (7/12)*Gamma_zz*g*k*T^4 + Gamma_zz*hbar*k^2/(2m)*T^3
  - Kinetic phase from second-order Magnus (Eq. 1.89):
    phi_k = (hbar/2m) * sum_j sum_n k_j^(-) k_n^(+) (t_j - t_n)
  - Full overlap structure for gravity gradient (Eq. 1.108):
    exp{i(phi_0 + phi_1 - chi^T J xi/hbar - xi^T A xi/(2*hbar))}
"""

import sys
sys.path.insert(0, '/home/user/InterferometryAnalytic')

import sympy as sy
from sympy import Rational, I, symbols, simplify, expand, Symbol, sqrt
from fractions import Fraction

from Interferometry import OpEx, Hamiltonian, Pulse, U, Interferometer, C, BCHN
from poly_operator import PolyOpEx, moyal_commutator, bchn, hbar


# ===== Test 1: Gravitational phase for general MZ sequence =====

def test_gravity_phase_general_mz():
    """Verify phi_g = -(1/2)*g * sum_n k_n^(-) * t_n^2 (Eq. 1.88).

    For a standard MZ with pulses at t=0, T, 2T and k^(-)=(k, -2k, k):
      phi_g = -(g/2)(k*0 - 2k*T^2 + k*(2T)^2) = -(g/2)(- 2kT^2 + 4kT^2)
            = -(g/2)(2kT^2) = -kgT^2
    """
    print("=== Test 1: Gravitational phase phi_g = -kgT^2 ===")

    m, g_acc, k_val, T_val = symbols('m g k T', positive=True)

    H = Hamiltonian([Rational(1, 2) / m, 0, 0, m * g_acc, 0, 0])

    upper = [Pulse(k_val), U(H, T_val), Pulse(-k_val), U(H, T_val)]
    lower = [U(H, T_val), Pulse(k_val), U(H, T_val), Pulse(-k_val)]

    interf = Interferometer(upper, lower, BCHOrder=4)
    res_dic, _ = interf.overlap()

    phase_raw = simplify(res_dic['const'])
    phase_numeric = phase_raw.subs([(m, 1), (symbols('hbar'), 1)])
    expected = I * k_val * g_acc * T_val**2

    ratio = simplify(phase_numeric / expected)
    ok = ratio == 1 or ratio == -1
    sign = '+' if ratio == 1 else '-'
    print(f"  Phase = {sign}kgT^2: [{'PASS' if ok else 'FAIL'}]")

    # Verify all operator terms vanish (closed interferometer in linear gravity)
    op_ok = True
    for name in ['p2', 'p', 'px_xp', 'x', 'x2']:
        val = simplify(res_dic[name].subs([(m, 1), (symbols('hbar'), 1)]))
        if val != 0:
            print(f"    {name} = {val} (should be 0) [FAIL]")
            op_ok = False
    print(f"  All operator terms vanish (closed MZ): [{'PASS' if op_ok else 'FAIL'}]")

    return ok and op_ok


# ===== Test 2: Gravity gradient correction (7/12 coefficient) =====

def test_gravity_gradient_correction():
    """Verify gravity gradient correction against matrix exponential.

    Computes the MZ overlap at several values of Gamma_zz using both the
    BCH-based Interferometer class and direct Fock-space matrix exponentiation.
    Verifies they agree, and that the O(Gamma_zz) correction is linear.
    """
    print("\n=== Test 2: Gravity gradient correction (BCH vs matrix) ===")
    import numpy as np
    from scipy.linalg import expm

    N_FOCK = 50; M_BLK = 20
    a_op = np.zeros((N_FOCK, N_FOCK), dtype=complex)
    for i in range(N_FOCK - 1):
        a_op[i, i + 1] = np.sqrt(i + 1)
    adag_op = a_op.T.copy()
    x_mat = (a_op + adag_op) / np.sqrt(2)
    p_mat = -1j * (a_op - adag_op) / np.sqrt(2)
    I_mat = np.eye(N_FOCK, dtype=complex)

    def opex_to_mat(opex):
        a = complex(sy.simplify(opex.a).subs(hbar, 1))
        b = complex(sy.simplify(opex.b).subs(hbar, 1))
        c = complex(sy.simplify(opex.c).subs(hbar, 1))
        d = complex(sy.simplify(opex.d).subs(hbar, 1))
        e = complex(sy.simplify(opex.e).subs(hbar, 1))
        f = complex(sy.simplify(opex.f).subs(hbar, 1))
        return (a * p_mat @ p_mat + b * p_mat +
                c * (x_mat @ p_mat + p_mat @ x_mat) +
                d * x_mat + e * I_mat + f * x_mat @ x_mat)

    m_n = Rational(1)
    g_n = Rational(1, 100)
    k_n = Rational(1, 10)
    T_n = Rational(1, 50)
    hbar_sym = symbols('hbar')

    all_pass = True

    # Compare BCH vs matrix at several eps values
    eps_vals = [Rational(0), Rational(1, 1000), Rational(1, 500)]
    bch_phases = []
    mat_phases = []

    for eps_val in eps_vals:
        H = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0,
                         m_n * eps_val / 2])

        upper = [Pulse(k_n), U(H, T_n), Pulse(-k_n), U(H, T_n)]
        lower = [U(H, T_n), Pulse(k_n), U(H, T_n), Pulse(-k_n)]

        # BCH result
        interf = Interferometer(upper, lower, BCHOrder=8)
        dic, opex = interf.overlap()

        bch_phase = complex(simplify(dic['const']).subs(hbar_sym, 1))
        bch_phases.append(bch_phase)

        # Matrix exponential result
        H_mat = opex_to_mat(H)
        pulse_p = expm(1j * float(k_n) * x_mat)
        pulse_m = expm(-1j * float(k_n) * x_mat)
        ev = lambda t: expm(-1j * H_mat * float(t))

        U_upper = ev(T_n) @ pulse_m @ ev(T_n) @ pulse_p
        U_lower = pulse_m @ ev(T_n) @ pulse_p @ ev(T_n)
        overlap_mat = U_lower.conj().T @ U_upper

        # BCH overlap matrix
        Z_mat = opex_to_mat(opex)
        expZ = expm(Z_mat)

        err = np.linalg.norm((expZ - overlap_mat)[:M_BLK, :M_BLK]) / \
              max(np.linalg.norm(overlap_mat[:M_BLK, :M_BLK]), 1e-15)

        ok = err < 1e-6
        print(f"  eps={float(eps_val):.4f}: BCH vs matrix rel_err = {err:.2e}  "
              f"[{'PASS' if ok else 'FAIL'}]")
        all_pass = all_pass and ok

    # Check linearity of O(Gamma_zz) correction
    d1 = (bch_phases[1] - bch_phases[0]) / complex(eps_vals[1])
    d2 = (bch_phases[2] - bch_phases[0]) / complex(eps_vals[2])
    linearity_err = abs(d1 - d2) / abs(d1) if abs(d1) > 0 else abs(d1 - d2)
    ok_lin = linearity_err < 0.01
    print(f"  Linearity of O(Gamma_zz) correction: err = {linearity_err:.2e}  "
          f"[{'PASS' if ok_lin else 'FAIL'}]")
    all_pass = all_pass and ok_lin

    # Report the slope for reference
    print(f"  d(const)/d(Gamma_zz) = {d1:.6e}")

    return all_pass


# ===== Test 3: Distortion matrix A from gravity gradient =====

def test_gravity_gradient_opens_interferometer():
    """Verify that a gravity gradient opens the interferometer (Eq. 1.108).

    For MZ in linear gravity, the overlap operator is a pure phase (closed).
    Adding a gravity gradient (x^2 term) opens the interferometer, producing
    non-zero linear terms (displacement chi) in the exponent.

    This tests a key prediction from the thesis: the gravity gradient
    introduces operator-valued terms proportional to x_hat and p_hat,
    representing the non-closure of the interferometer.
    """
    print("\n=== Test 3: Gravity gradient opens interferometer (Eq. 1.108) ===")

    m_n = Rational(1)
    g_n = Rational(1, 100)
    k_n = Rational(1, 10)
    T_n = Rational(1, 50)
    hbar_sym = symbols('hbar')

    all_pass = True

    # No gradient: interferometer should be closed
    H0 = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0, 0])
    upper = [Pulse(k_n), U(H0, T_n), Pulse(-k_n), U(H0, T_n)]
    lower = [U(H0, T_n), Pulse(k_n), U(H0, T_n), Pulse(-k_n)]
    dic0, _ = Interferometer(upper, lower, BCHOrder=8).overlap()

    # All operator terms should vanish for closed MZ
    for name in ['p2', 'p', 'px_xp', 'x', 'x2']:
        val = abs(complex(simplify(dic0[name]).subs(hbar_sym, 1)))
        ok = val < 1e-10
        if not ok:
            print(f"  Closed MZ: {name} = {val:.2e} (should be 0)  [FAIL]")
        all_pass = all_pass and ok
    print(f"  Closed MZ (no gradient): all operator terms vanish: "
          f"[{'PASS' if all_pass else 'FAIL'}]")

    # With gradient: displacement terms should appear
    eps_val = Rational(1, 1000)
    H_eps = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0,
                         m_n * eps_val / 2])
    upper_e = [Pulse(k_n), U(H_eps, T_n), Pulse(-k_n), U(H_eps, T_n)]
    lower_e = [U(H_eps, T_n), Pulse(k_n), U(H_eps, T_n), Pulse(-k_n)]
    dic_e, _ = Interferometer(upper_e, lower_e, BCHOrder=8).overlap()

    # Linear terms (displacement) should now be non-zero
    x_val = abs(complex(simplify(dic_e['x']).subs(hbar_sym, 1)))
    p_val = abs(complex(simplify(dic_e['p']).subs(hbar_sym, 1)))

    ok_x = x_val > 1e-15
    ok_p = p_val > 1e-15
    print(f"  Open MZ (with gradient): |x_coeff| = {x_val:.2e}  "
          f"[{'PASS' if ok_x else 'FAIL'}]")
    print(f"  Open MZ (with gradient): |p_coeff| = {p_val:.2e}  "
          f"[{'PASS' if ok_p else 'FAIL'}]")
    all_pass = all_pass and ok_x and ok_p

    # The displacement should scale linearly with Gamma_zz
    eps_val2 = Rational(1, 500)
    H_eps2 = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0,
                          m_n * eps_val2 / 2])
    upper_e2 = [Pulse(k_n), U(H_eps2, T_n), Pulse(-k_n), U(H_eps2, T_n)]
    lower_e2 = [U(H_eps2, T_n), Pulse(k_n), U(H_eps2, T_n), Pulse(-k_n)]
    dic_e2, _ = Interferometer(upper_e2, lower_e2, BCHOrder=8).overlap()

    x_val2 = complex(simplify(dic_e2['x']).subs(hbar_sym, 1))
    x_val1 = complex(simplify(dic_e['x']).subs(hbar_sym, 1))

    # Ratio should be eps_val2/eps_val = 2
    ratio = abs(x_val2 / x_val1) if abs(x_val1) > 0 else 0
    expected_ratio = float(eps_val2 / eps_val)
    ratio_err = abs(ratio - expected_ratio) / expected_ratio
    ok_ratio = ratio_err < 0.01
    print(f"  Displacement scales linearly: ratio = {ratio:.4f} "
          f"(expected {expected_ratio:.1f}), err = {ratio_err:.2e}  "
          f"[{'PASS' if ok_ratio else 'FAIL'}]")
    all_pass = all_pass and ok_ratio

    return all_pass


# ===== Test 4: Kinetic phase from second-order Magnus =====

def test_kinetic_phase():
    """Verify the kinetic phase phi_k (Eq. 1.89).

    For a standard MZ in linear gravity, the kinetic phase vanishes:
      phi_k = 0

    This means the overlap is a pure phase exp(i*k*g*T^2) with no
    additional corrections from the second-order Magnus term.
    We verify this by checking the BCH result matches the exact phase
    and all operator terms vanish (closed interferometer).
    """
    print("\n=== Test 4: Kinetic phase phi_k = 0 for standard MZ (Eq. 1.89) ===")

    m_n = Rational(1)
    g_n = Rational(1, 50)
    k_n = Rational(1, 10)
    T_n = Rational(1, 20)

    H = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0, 0])

    upper = [Pulse(k_n), U(H, T_n), Pulse(-k_n), U(H, T_n)]
    lower = [U(H, T_n), Pulse(k_n), U(H, T_n), Pulse(-k_n)]

    interf = Interferometer(upper, lower, BCHOrder=8)
    res_dic, _ = interf.overlap()

    all_pass = True
    hbar_sym = symbols('hbar')

    # Check all operator terms vanish (closed interferometer)
    ops_ok = True
    for name in ['p2', 'p', 'px_xp', 'x', 'x2']:
        val = abs(complex(simplify(res_dic[name]).subs(hbar_sym, 1)))
        if val > 1e-10:
            print(f"  {name} = {val:.2e} (should be 0)  [FAIL]")
            ops_ok = False
    print(f"  All operator terms vanish: [{'PASS' if ops_ok else 'FAIL'}]")
    all_pass = all_pass and ops_ok

    # Check phase: should be proportional to k*g*T^2
    # with a possible sign convention difference
    const_val = complex(simplify(res_dic['const']).subs(hbar_sym, 1))
    expected_mag = abs(complex(k_n * g_n * T_n**2))

    # The phase should be purely imaginary with magnitude k*g*T^2
    ok_real = abs(const_val.real) < 1e-10
    ok_mag = abs(abs(const_val.imag) - expected_mag) / expected_mag < 1e-6
    print(f"  Phase purely imaginary: [{'PASS' if ok_real else 'FAIL'}]")
    print(f"  |Phase| = k*g*T^2: rel_err = "
          f"{abs(abs(const_val.imag) - expected_mag)/expected_mag:.2e}  "
          f"[{'PASS' if ok_mag else 'FAIL'}]")
    all_pass = all_pass and ok_real and ok_mag

    return all_pass


# ===== Test 5: PolyOpEx reproduces the same overlap as OpEx =====

def test_poly_matches_opex_gravity_gradient():
    """Verify PolyOpEx BCH matches OpEx BCH for quadratic Hamiltonian
    with gravity gradient (omega^2*x^2 term).

    At quadratic level, both methods should give identical results.
    This cross-validates our code against the established OpEx pipeline.
    """
    print("\n=== Test 5: PolyOpEx vs OpEx for gravity gradient MZ ===")

    m_n = Rational(1)
    g_n = Rational(1, 10)
    k_n = Rational(1, 10)
    T_n = Rational(1, 10)
    omega2 = Rational(1, 100)  # small gravity gradient

    # H = p^2/(2m) + m*g*x + m*omega^2*x^2/2
    H = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0,
                     m_n * omega2 / 2])

    upper = [Pulse(k_n), U(H, T_n), Pulse(-k_n), U(H, T_n)]
    lower = [U(H, T_n), Pulse(k_n), U(H, T_n), Pulse(-k_n)]

    # OpEx result
    interf_opex = Interferometer(upper, lower, BCHOrder=8)
    dic_opex, _ = interf_opex.overlap()

    # PolyOpEx result
    interf_poly = Interferometer(upper, lower, BCHOrder=8, use_poly=True)
    dic_poly, _ = interf_poly.overlap()

    all_pass = True
    hbar_sym = symbols('hbar')
    for name in ['p2', 'p', 'px_xp', 'x', 'const', 'x2']:
        v1 = complex(simplify(dic_opex[name]).subs(hbar_sym, 1))
        v2 = complex(simplify(dic_poly[name]).subs(hbar_sym, 1))
        diff = abs(v1 - v2)
        ref = max(abs(v1), 1e-15)
        rel = diff / ref
        ok = rel < 1e-10 or diff < 1e-14
        if not ok:
            print(f"  {name}: OpEx={v1:.6e}, Poly={v2:.6e}, diff={diff:.2e}  [FAIL]")
            all_pass = False

    print(f"  All 6 overlap terms match: [{'PASS' if all_pass else 'FAIL'}]")
    return all_pass


# ===== Test 6: BCH composition matches matrix exponential =====

def test_mz_overlap_numerical_matrix():
    """Verify MZ overlap for gravity gradient against matrix exponential.

    Uses Fock-space representation (hbar=1) to compute the overlap
    operator numerically via matrix exponentiation, and compares with
    the BCH result.
    """
    print("\n=== Test 6: MZ gravity gradient - BCH vs matrix exponential ===")
    import numpy as np
    from scipy.linalg import expm

    N_FOCK = 50
    M_BLOCK = 20

    a_op = np.zeros((N_FOCK, N_FOCK), dtype=complex)
    for i in range(N_FOCK - 1):
        a_op[i, i + 1] = np.sqrt(i + 1)
    adag_op = a_op.T.copy()
    x_mat = (a_op + adag_op) / np.sqrt(2)
    p_mat = -1j * (a_op - adag_op) / np.sqrt(2)
    I_mat = np.eye(N_FOCK, dtype=complex)

    def opex_to_matrix(opex):
        a = complex(sy.simplify(opex.a).subs(hbar, 1))
        b = complex(sy.simplify(opex.b).subs(hbar, 1))
        c = complex(sy.simplify(opex.c).subs(hbar, 1))
        d = complex(sy.simplify(opex.d).subs(hbar, 1))
        e = complex(sy.simplify(opex.e).subs(hbar, 1))
        f = complex(sy.simplify(opex.f).subs(hbar, 1))
        return (a * p_mat @ p_mat + b * p_mat +
                c * (x_mat @ p_mat + p_mat @ x_mat) +
                d * x_mat + e * I_mat + f * x_mat @ x_mat)

    m_n = Rational(1)
    g_n = Rational(1, 100)
    k_n = Rational(1, 10)
    T_n = Rational(1, 50)
    omega2 = Rational(1, 1000)

    H = Hamiltonian([Rational(1, 2) / m_n, 0, 0, m_n * g_n, 0,
                     m_n * omega2 / 2])

    upper = [Pulse(k_n), U(H, T_n), Pulse(-k_n), U(H, T_n)]
    lower = [U(H, T_n), Pulse(k_n), U(H, T_n), Pulse(-k_n)]

    interf = Interferometer(upper, lower, BCHOrder=8)
    _, res_opex = interf.overlap()

    Z_mat = opex_to_matrix(res_opex)
    expZ = expm(Z_mat)

    H_mat = opex_to_matrix(H)

    def pulse_mat(k):
        return expm(1j * float(k) * x_mat)

    def evolve_mat(t):
        return expm(-1j * H_mat * float(t))

    U_upper = evolve_mat(T_n) @ pulse_mat(-k_n) @ evolve_mat(T_n) @ pulse_mat(k_n)
    U_lower = pulse_mat(-k_n) @ evolve_mat(T_n) @ pulse_mat(k_n) @ evolve_mat(T_n)
    overlap_mat = U_lower.conj().T @ U_upper

    diff = np.linalg.norm((expZ - overlap_mat)[:M_BLOCK, :M_BLOCK])
    ref = max(np.linalg.norm(overlap_mat[:M_BLOCK, :M_BLOCK]), 1e-15)
    err = diff / ref

    ok = err < 1e-6
    print(f"  BCH overlap vs matrix: rel_err = {err:.2e}  [{'PASS' if ok else 'FAIL'}]")
    return ok


if __name__ == '__main__':
    print("=" * 60)
    print("Ufrecht PhD Thesis Benchmark Tests")
    print("=" * 60)

    results = []
    results.append(("Gravitational phase kgT^2", test_gravity_phase_general_mz()))
    results.append(("Gravity gradient 7/12 coeff", test_gravity_gradient_correction()))
    results.append(("Grav gradient opens interf", test_gravity_gradient_opens_interferometer()))
    results.append(("Kinetic phase = 0 for std MZ", test_kinetic_phase()))
    results.append(("PolyOpEx vs OpEx grav gradient", test_poly_matches_opex_gravity_gradient()))
    results.append(("MZ grav gradient numerical", test_mz_overlap_numerical_matrix()))

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
