"""
Generalized polynomial operator algebra for atom interferometry.

Extends the quadratic operator framework in Interferometry.py to arbitrary
polynomial degree using the Weyl symbol / Moyal bracket formalism.

A PolyOpEx represents an operator as a polynomial in x and p via its Weyl
symbol: W(x,p) = Σ c_{a,b} x^a p^b, which corresponds to the Weyl-ordered
operator Σ c_{a,b} Sym(x^a p^b).

The commutator is the Moyal bracket, which is exact for polynomials and
reduces to the Poisson bracket (times iℏ) at leading order.

Perturbative order tracking:
    Each monomial can carry an integer perturbative weight.  When a
    max_pert_order is set, terms exceeding that weight are automatically
    dropped in all operations (commutators, addition, BCH).  This makes
    high-order BCH with non-quadratic Hamiltonians tractable.

References:
    Arnal, Casas, Chiralt, Mediterr. J. Math. 18, 53 (2021) [arXiv:2006.15869]
"""

import math
from fractions import Fraction
import sympy as sy
from sympy import I, Rational, factorial
from collections import defaultdict

hbar = sy.symbols('hbar')

# ---------------------------------------------------------------------------
# Falling factorial (integer arguments only)
# ---------------------------------------------------------------------------

def _falling_factorial(n, k):
    """n · (n-1) · ... · (n-k+1).  Returns 0 if n < k or k < 0."""
    if k < 0 or n < k:
        return 0
    r = 1
    for i in range(k):
        r *= (n - i)
    return r


def _binom(n, k):
    """Binomial coefficient C(n, k) for non-negative integers."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


# ---------------------------------------------------------------------------
# Fast representation helpers: decompose sympy coeffs into Fraction × (iℏ)^h
# ---------------------------------------------------------------------------

def _decompose_ihbar(expr):
    """Decompose a sympy expression into {h: Fraction} where expr = Σ frac_h * (i*hbar)^h.

    Returns None if the expression contains symbols other than hbar or
    cannot be decomposed into rational × (i*hbar)^h form.
    """
    expr = sy.expand(expr)
    if expr == 0:
        return {0: Fraction(0)}

    # Check for unexpected free symbols
    free = expr.free_symbols - {hbar}
    if free:
        return None

    # Collect as polynomial in hbar
    try:
        poly = sy.Poly(expr, hbar, domain='QQ_I')  # Gaussian rationals
    except (sy.polys.polyerrors.PolynomialError, sy.polys.polyerrors.GeneratorsError):
        return None  # not a polynomial in hbar (e.g., contains 1/hbar)
    result = {}
    for monom, coeff in poly.as_dict().items():
        h = monom[0]  # power of hbar
        # coeff is a Gaussian rational: a + b*I
        # We need coeff * (i*hbar)^h to equal coeff_original * hbar^h
        # Since (i*hbar)^h = i^h * hbar^h, we need frac = coeff / i^h
        # i^0=1, i^1=i, i^2=-1, i^3=-i, i^4=1, ...
        i_power = h % 4
        # Divide coeff by i^h
        re_c = sy.re(coeff)
        im_c = sy.im(coeff)
        if i_power == 0:
            frac_re, frac_im = re_c, im_c
        elif i_power == 1:
            # coeff / i = coeff * (-i) = im_c - i*re_c
            frac_re, frac_im = im_c, -re_c
        elif i_power == 2:
            # coeff / i^2 = coeff / (-1) = -re_c - i*im_c
            frac_re, frac_im = -re_c, -im_c
        elif i_power == 3:
            # coeff / i^3 = coeff * i = -im_c + i*re_c
            frac_re, frac_im = -im_c, re_c

        if frac_im != 0:
            return None  # not a real rational × (i*hbar)^h
        result[h] = Fraction(int(sy.numer(frac_re)), int(sy.denom(frac_re)))

    return result


def _fast_to_sympy(fast):
    """Convert fast {(a, b, h): Fraction} back to {(a, b): sympy_expr}."""
    result = defaultdict(lambda: sy.Integer(0))
    for (a, b, h), frac in fast.items():
        if frac != 0:
            result[(a, b)] += sy.Rational(frac.numerator, frac.denominator) * (I * hbar)**h
    return {k: v for k, v in result.items() if v != 0}


# ---------------------------------------------------------------------------
# PolyOpEx class
# ---------------------------------------------------------------------------

class PolyOpEx:
    """Polynomial operator expression in Weyl-ordered basis.

    Stores Weyl symbol coefficients: {(a, b): c_ab} representing
    the operator  Σ c_{a,b} Sym(x^a p^b).

    Internally uses a dual representation for performance:
    - ``coeffs``: {(a, b): sympy_expr} — general symbolic coefficients
    - ``_fast``: {(a, b, h): Fraction} — hbar-separated rational coefficients

    When all coefficients are rational (no free symbols other than hbar),
    the fast path is used for commutator arithmetic, giving ~10-50x speedup
    by avoiding SymPy's expression tree overhead.

    Parameters
    ----------
    coeffs : dict
        Mapping (a, b) -> sympy expression.  Zero entries are dropped.
    max_degree : int or None
        If set, monomials with a+b > max_degree are silently truncated.
    pert_order : dict or int or None
        Perturbative order for each monomial.
        - dict {(a, b): int} — per-monomial weights
        - int — uniform weight for all monomials
        - None — no perturbative tracking (weight 0 for all)
    max_pert_order : int or None
        If set, monomials with perturbative order > max_pert_order are
        silently dropped in all operations.
    pert_symbol : sympy Symbol or None
        If set, a symbolic perturbation parameter (e.g., ε).  The method
        ``collect_pert_orders()`` returns the operator grouped by powers of
        this symbol, making perturbative structure explicit in the output.
    """

    def __init__(self, coeffs, max_degree=None, pert_order=None, max_pert_order=None,
                 pert_symbol=None):
        self.max_degree = max_degree
        self.max_pert_order = max_pert_order
        self.pert_symbol = pert_symbol

        # Normalize pert_order to a dict
        if pert_order is None:
            self._pert_order = {}    # default: 0 for everything
        elif isinstance(pert_order, int):
            self._pert_order = {k: pert_order for k in coeffs if coeffs[k] != 0}
        else:
            self._pert_order = dict(pert_order)

        self.coeffs = {}
        for (a, b), c in coeffs.items():
            if c == 0:
                continue
            if max_degree is not None and a + b > max_degree:
                continue
            po = self._pert_order.get((a, b), 0)
            if max_pert_order is not None and po > max_pert_order:
                continue
            self.coeffs[(a, b)] = c

        # Clean up pert_order to match actual coeffs
        self._pert_order = {k: self._pert_order.get(k, 0) for k in self.coeffs}

        # Try to build fast representation
        self._fast = None
        self._try_build_fast()

    def _try_build_fast(self):
        """Try to decompose coefficients into {(a, b, h): Fraction}.

        Each coefficient is expected to be a polynomial in (i*hbar) with
        rational coefficients.  If any coefficient contains other symbols,
        falls back to the slow SymPy path.
        """
        fast = {}
        for (a, b), c in self.coeffs.items():
            decomposed = _decompose_ihbar(c)
            if decomposed is None:
                self._fast = None
                return
            for h, frac in decomposed.items():
                if frac != 0:
                    fast[(a, b, h)] = frac
        self._fast = fast

    @classmethod
    def _from_fast(cls, fast, pert_order, max_degree=None, max_pert_order=None,
                   pert_symbol=None):
        """Construct from fast representation, bypassing SymPy."""
        obj = cls.__new__(cls)
        obj.max_degree = max_degree
        obj.max_pert_order = max_pert_order
        obj.pert_symbol = pert_symbol

        # Filter and store fast coeffs
        filtered = {}
        for (a, b, h), frac in fast.items():
            if frac == 0:
                continue
            if max_degree is not None and a + b > max_degree:
                continue
            po = pert_order.get((a, b), 0)
            if max_pert_order is not None and po > max_pert_order:
                continue
            filtered[(a, b, h)] = frac

        obj._fast = filtered
        obj._pert_order = {k: pert_order.get(k, 0)
                           for k in {(a, b) for (a, b, h) in filtered}}

        # Build sympy coeffs lazily (only when accessed)
        obj._coeffs_cache = None
        return obj

    @property
    def coeffs(self):
        if hasattr(self, '_coeffs_cache'):
            if self._coeffs_cache is None:
                self._coeffs_cache = _fast_to_sympy(self._fast)
            return self._coeffs_cache
        return self._coeffs_dict

    @coeffs.setter
    def coeffs(self, value):
        self._coeffs_dict = value

    def get_pert_order(self, key):
        """Get perturbative order for monomial (a, b)."""
        return self._pert_order.get(key, 0)

    # --- constructors -------------------------------------------------------

    @classmethod
    def zero(cls, max_degree=None, max_pert_order=None, pert_symbol=None):
        return cls({}, max_degree, max_pert_order=max_pert_order,
                   pert_symbol=pert_symbol)

    @classmethod
    def from_opex(cls, opex, max_degree=None, max_pert_order=None, pert_symbol=None):
        """Convert from OpEx [a,b,c,d,e,f] = a·p² + b·p + c·(xp+px) + d·x + e + f·x²."""
        coeffs = {}
        if opex.a != 0: coeffs[(0, 2)] = opex.a
        if opex.b != 0: coeffs[(0, 1)] = opex.b
        if opex.c != 0: coeffs[(1, 1)] = 2 * opex.c   # (xp+px) = 2·Sym(xp)
        if opex.d != 0: coeffs[(1, 0)] = opex.d
        if opex.e != 0: coeffs[(0, 0)] = opex.e
        if opex.f != 0: coeffs[(2, 0)] = opex.f
        md = max_degree if max_degree is not None else 2
        return cls(coeffs, md, pert_order=0, max_pert_order=max_pert_order,
                   pert_symbol=pert_symbol)

    def to_opex(self):
        """Convert back to OpEx.  Raises ValueError if degree > 2."""
        from Interferometry import OpEx
        for (a, b) in self.coeffs:
            if a + b > 2 and self.coeffs[(a, b)] != 0:
                raise ValueError(
                    f"Cannot convert degree-{a+b} term x^{a}p^{b} to OpEx (max degree 2)")
        return OpEx([
            self.coeffs.get((0, 2), 0),
            self.coeffs.get((0, 1), 0),
            self.coeffs.get((1, 1), 0) / 2,   # Sym(xp) → (xp+px)/2
            self.coeffs.get((1, 0), 0),
            self.coeffs.get((0, 0), 0),
            self.coeffs.get((2, 0), 0),
        ])

    # --- properties ---------------------------------------------------------

    @property
    def degree(self):
        if not self.coeffs:
            return 0
        return max(a + b for (a, b) in self.coeffs)

    @property
    def pert_degree(self):
        """Maximum perturbative order across all monomials."""
        if not self._pert_order:
            return 0
        return max(self._pert_order.values()) if self._pert_order else 0

    # --- truncation helpers -------------------------------------------------

    def _merge_truncation(self, other):
        """Determine combined max_degree, max_pert_order, pert_symbol from two operands."""
        md = self.max_degree
        if md is None:
            md = other.max_degree
        elif other.max_degree is not None:
            md = max(md, other.max_degree)

        mpo = self.max_pert_order
        if mpo is None:
            mpo = other.max_pert_order
        elif other.max_pert_order is not None:
            mpo = max(mpo, other.max_pert_order)

        ps = self.pert_symbol if self.pert_symbol is not None else other.pert_symbol

        return md, mpo, ps

    # --- arithmetic ---------------------------------------------------------

    @property
    def is_fast(self):
        return self._fast is not None

    def __add__(self, other):
        md, mpo, ps = self._merge_truncation(other)

        # Fast path: both operands in fast representation
        if self._fast is not None and other._fast is not None:
            result_fast = dict(self._fast)
            result_po = dict(self._pert_order)
            for k, v in other._fast.items():
                result_fast[k] = result_fast.get(k, Fraction(0)) + v
            for k in {(a, b) for (a, b, h) in other._fast}:
                if k in result_po:
                    result_po[k] = min(result_po[k], other.get_pert_order(k))
                else:
                    result_po[k] = other.get_pert_order(k)
            return PolyOpEx._from_fast(result_fast, result_po, md, mpo,
                                       pert_symbol=ps)

        # Slow path: general sympy coefficients
        result_coeffs = dict(self.coeffs)
        result_po = dict(self._pert_order)

        for k, v in other.coeffs.items():
            result_coeffs[k] = result_coeffs.get(k, 0) + v
            if k in result_po:
                result_po[k] = min(result_po[k], other.get_pert_order(k))
            else:
                result_po[k] = other.get_pert_order(k)

        return PolyOpEx(result_coeffs, md, pert_order=result_po,
                        max_pert_order=mpo, pert_symbol=ps)

    def __sub__(self, other):
        return self + (other * (-1))

    def __neg__(self):
        return self * (-1)

    def __mul__(self, scalar):
        # Fast path: scalar is a rational number
        if self._fast is not None:
            try:
                frac_scalar = Fraction(scalar)
                return PolyOpEx._from_fast(
                    {k: v * frac_scalar for k, v in self._fast.items()},
                    dict(self._pert_order), self.max_degree, self.max_pert_order,
                    pert_symbol=self.pert_symbol)
            except (TypeError, ValueError, ZeroDivisionError):
                pass

        return PolyOpEx(
            {k: v * scalar for k, v in self.coeffs.items()},
            self.max_degree,
            pert_order=dict(self._pert_order),
            max_pert_order=self.max_pert_order,
            pert_symbol=self.pert_symbol)

    def __rmul__(self, scalar):
        return self * scalar

    def __eq__(self, other):
        if not isinstance(other, PolyOpEx):
            return NotImplemented
        keys = set(self.coeffs) | set(other.coeffs)
        for k in keys:
            diff = sy.simplify(self.coeffs.get(k, 0) - other.coeffs.get(k, 0))
            if diff != 0:
                return False
        return True

    def simplify(self):
        if self._fast is not None:
            return self  # fast coeffs are already simplified (exact rationals)
        return PolyOpEx(
            {k: sy.simplify(v) for k, v in self.coeffs.items()},
            self.max_degree,
            pert_order=dict(self._pert_order),
            max_pert_order=self.max_pert_order,
            pert_symbol=self.pert_symbol)

    def expand(self):
        return PolyOpEx(
            {k: sy.expand(v) for k, v in self.coeffs.items()},
            self.max_degree,
            pert_order=dict(self._pert_order),
            max_pert_order=self.max_pert_order,
            pert_symbol=self.pert_symbol)

    def truncate(self, max_degree=None, max_pert_order=None):
        """Return a new PolyOpEx with tighter truncation bounds."""
        md = max_degree if max_degree is not None else self.max_degree
        mpo = max_pert_order if max_pert_order is not None else self.max_pert_order
        return PolyOpEx(dict(self.coeffs), md,
                        pert_order=dict(self._pert_order), max_pert_order=mpo,
                        pert_symbol=self.pert_symbol)

    def collect_pert_orders(self):
        """Return a dict {order: PolyOpEx} grouping terms by perturbative order.

        Useful for inspecting the perturbative structure of the result.
        If ``pert_symbol`` is set, the returned PolyOpExs have their
        coefficients multiplied by pert_symbol^order.
        """
        groups = defaultdict(dict)
        po_groups = defaultdict(dict)
        for (a, b), c in self.coeffs.items():
            po = self.get_pert_order((a, b))
            groups[po][(a, b)] = c
            po_groups[po][(a, b)] = po

        result = {}
        eps = self.pert_symbol
        for order, coeffs in sorted(groups.items()):
            p = PolyOpEx(coeffs, self.max_degree,
                         pert_order=po_groups[order],
                         max_pert_order=self.max_pert_order,
                         pert_symbol=self.pert_symbol)
            if eps is not None and order > 0:
                p = PolyOpEx(
                    {k: v * eps**order for k, v in p.coeffs.items()},
                    p.max_degree,
                    pert_order=dict(p._pert_order),
                    max_pert_order=p.max_pert_order,
                    pert_symbol=self.pert_symbol)
            result[order] = p
        return result

    # --- display ------------------------------------------------------------

    def __repr__(self):
        if not self.coeffs:
            return 'PolyOpEx(0)'
        terms = []
        for (a, b) in sorted(self.coeffs):
            c = self.coeffs[(a, b)]
            mon = ''
            if a == 1: mon += 'x'
            elif a > 1: mon += f'x^{a}'
            if b == 1: mon += 'p'
            elif b > 1: mon += f'p^{b}'
            if not mon:
                mon = '1'
            po = self._pert_order.get((a, b), 0)
            po_str = f'[O({po})]' if po > 0 else ''
            terms.append(f'({c})*{mon}{po_str}')
        return 'PolyOpEx(' + ' + '.join(terms) + ')'


# ---------------------------------------------------------------------------
# Moyal bracket (commutator for Weyl symbols)
# ---------------------------------------------------------------------------

def _moyal_bracket_sum(a1, b1, a2, b2, N):
    """Compute the inner bracket sum for Moyal order N.  Pure Python ints."""
    total = 0
    for k in range(N + 1):
        ff1x = _falling_factorial(a1, N - k)
        ff1p = _falling_factorial(b1, k)
        ff2x = _falling_factorial(a2, k)
        ff2p = _falling_factorial(b2, N - k)
        if ff1x == 0 or ff1p == 0 or ff2x == 0 or ff2p == 0:
            continue
        total += (-1)**k * math.comb(N, k) * ff1x * ff1p * ff2x * ff2p
    return total


def moyal_commutator(A, B):
    """Commutator [A, B] via the Moyal bracket of Weyl symbols.

    For Weyl symbols f and g, the commutator Weyl symbol is:

        [f, g]_M = Σ_{s≥0}  2·(iℏ/2)^{2s+1} / (2s+1)!
                   × Σ_k (-1)^k C(2s+1,k) (∂_x^{N-k} ∂_p^k f)(∂_x^k ∂_p^{N-k} g)

    where N = 2s+1.  The series terminates for polynomial f, g.

    Uses a fast pure-Python path when both inputs have rational coefficients
    (no symbolic variables other than hbar), avoiding SymPy overhead entirely.
    """
    md, mpo, ps = A._merge_truncation(B)

    # Fast path: both operands have rational + hbar representation
    if A._fast is not None and B._fast is not None:
        return _moyal_commutator_fast(A, B, md, mpo, ps)

    # Slow path: general sympy coefficients
    return _moyal_commutator_sympy(A, B, md, mpo, ps)


def _moyal_commutator_fast(A, B, md, mpo, ps=None):
    """Fast Moyal commutator using pure Python Fraction arithmetic."""
    result_fast = defaultdict(Fraction)
    result_po = {}

    # Group A's fast coeffs by (a, b) for iteration
    a_by_ab = defaultdict(list)
    for (a1, b1, h1), frac1 in A._fast.items():
        a_by_ab[(a1, b1)].append((h1, frac1))

    b_by_ab = defaultdict(list)
    for (a2, b2, h2), frac2 in B._fast.items():
        b_by_ab[(a2, b2)].append((h2, frac2))

    for (a1, b1), a_terms in a_by_ab.items():
        po1 = A.get_pert_order((a1, b1))
        for (a2, b2), b_terms in b_by_ab.items():
            po2 = B.get_pert_order((a2, b2))
            po_out = po1 + po2

            if mpo is not None and po_out > mpo:
                continue

            max_N = min(a1 + b1, a2 + b2)

            for N in range(1, max_N + 1, 2):
                rx = a1 + a2 - N
                rp = b1 + b2 - N
                if rx < 0 or rp < 0:
                    continue
                if md is not None and rx + rp > md:
                    continue

                bracket_sum = _moyal_bracket_sum(a1, b1, a2, b2, N)
                if bracket_sum == 0:
                    continue

                # Prefactor as Fraction: 2 * bracket_sum / (2^N * N!)
                rat = Fraction(bracket_sum, (1 << (N - 1)) * math.factorial(N))

                # Multiply each pair of (h1, frac1) × (h2, frac2)
                for h1, frac1 in a_terms:
                    for h2, frac2 in b_terms:
                        h_out = h1 + h2 + N   # total (i*hbar) power
                        result_fast[(rx, rp, h_out)] += rat * frac1 * frac2

                # Track pert_order
                key = (rx, rp)
                if key in result_po:
                    result_po[key] = min(result_po[key], po_out)
                else:
                    result_po[key] = po_out

    return PolyOpEx._from_fast(dict(result_fast), result_po, md, mpo,
                               pert_symbol=ps)


def _moyal_commutator_sympy(A, B, md, mpo, ps=None):
    """Moyal commutator using SymPy symbolic arithmetic (general case)."""
    contributions = defaultdict(list)
    result_po = {}

    for (a1, b1), c1 in A.coeffs.items():
        po1 = A.get_pert_order((a1, b1))
        for (a2, b2), c2 in B.coeffs.items():
            po2 = B.get_pert_order((a2, b2))
            po_out = po1 + po2

            if mpo is not None and po_out > mpo:
                continue

            max_N = min(a1 + b1, a2 + b2)

            for N in range(1, max_N + 1, 2):
                rx = a1 + a2 - N
                rp = b1 + b2 - N
                if rx < 0 or rp < 0:
                    continue
                if md is not None and rx + rp > md:
                    continue

                bracket_sum = _moyal_bracket_sum(a1, b1, a2, b2, N)
                if bracket_sum == 0:
                    continue

                rat = Fraction(bracket_sum, (1 << (N - 1)) * math.factorial(N))
                contributions[(rx, rp, N)].append((rat, c1, c2))

                key = (rx, rp)
                if key in result_po:
                    result_po[key] = min(result_po[key], po_out)
                else:
                    result_po[key] = po_out

    result = defaultdict(lambda: sy.Integer(0))
    ihbar_cache = {}
    for (rx, rp, N) in contributions:
        if N not in ihbar_cache:
            ihbar_cache[N] = (I * hbar)**N

    for (rx, rp, N), contribs in contributions.items():
        ihbar_N = ihbar_cache[N]
        coeff_sum = sy.Integer(0)
        for rat, c1, c2 in contribs:
            coeff_sum += sy.Rational(rat.numerator, rat.denominator) * c1 * c2
        if coeff_sum != 0:
            result[(rx, rp)] += ihbar_N * coeff_sum

    return PolyOpEx(dict(result), md, pert_order=result_po, max_pert_order=mpo,
                    pert_symbol=ps)


# ---------------------------------------------------------------------------
# BCH word tables  (verified against Arnal, Casas, Chiralt 2020)
# ---------------------------------------------------------------------------

BCH_WORDS = {
    2: [('xy', Rational(1, 2))],
    3: [('xxy', Rational(1, 12)),
        ('yxy', Rational(-1, 12))],
    4: [('xyxy', Rational(-1, 24))],
    5: [('xxxxy', Rational(-1, 720)),
        ('xyxxy', Rational(-1, 120)),
        ('xyyxy', Rational(-1, 360)),
        ('yxxxy', Rational(1, 360)),
        ('yyxxy', Rational(1, 120)),
        ('yyyxy', Rational(1, 720))],
    6: [('xxyyxy', Rational(-1, 720)),
        ('xyyxxy', Rational(1, 240)),
        ('xyyyxy', Rational(1, 1440)),
        ('yxxxxy', Rational(1, 1440))],
    7: [('xxxxxxy', Rational(1, 30240)),
        ('xxyxxxy', Rational(1, 5040)),
        ('xxyyxxy', Rational(-1, 10080)),
        ('xyxxxxy', Rational(1, 10080)),
        ('xyxyxxy', Rational(1, 1008)),
        ('xyxyyxy', Rational(1, 5040)),
        ('xyyxxxy', Rational(-1, 7560)),
        ('xyyyxxy', Rational(1, 3360)),
        ('xyyyyxy', Rational(1, 10080)),
        ('yxxxxxy', Rational(-1, 10080)),
        ('yxyxxxy', Rational(-1, 1260)),
        ('yxyyxxy', Rational(-1, 1680)),
        ('yyxxxxy', Rational(1, 3360)),
        ('yyxyxxy', Rational(-1, 3360)),
        ('yyxyyxy', Rational(-1, 2520)),
        ('yyyxxxy', Rational(1, 7560)),
        ('yyyyxxy', Rational(1, 10080)),
        ('yyyyyxy', Rational(-1, 30240))],
    8: [('xxxyyyxy', Rational(-5, 24192)),
        ('xxyxyyxy', Rational(1, 2520)),
        ('xxyyyyxy', Rational(1, 20160)),
        ('xyxxyyxy', Rational(1, 15120)),
        ('xyxyyxxy', Rational(-1, 2016)),
        ('xyxyyyxy', Rational(-1, 20160)),
        ('xyyxyxxy', Rational(1, 20160)),
        ('xyyxyyxy', Rational(-1, 10080)),
        ('xyyyyyxy', Rational(-1, 60480)),
        ('yxxxxxxy', Rational(-1, 60480)),
        ('yxxxyxxy', Rational(1, 20160)),
        ('yxxyxxxy', Rational(-1, 5040)),
        ('yyxxxxxy', Rational(1, 20160))],
}


def _eval_word(word, X, Y, comm, cache):
    """Recursively evaluate a right-nested commutator word with memoization."""
    if word in cache:
        return cache[word]

    if len(word) == 2:
        assert word == 'xy', f"Unexpected base word: {word}"
        result = comm(X, Y)
    else:
        inner = _eval_word(word[1:], X, Y, comm, cache)
        arg = X if word[0] == 'x' else Y
        result = comm(arg, inner)

    cache[word] = result
    return result


# ---------------------------------------------------------------------------
# Generic BCH expansion
# ---------------------------------------------------------------------------

def bchn(X, Y, n_order=8, comm=None):
    """Baker-Campbell-Hausdorff expansion: exp(X)·exp(Y) = exp(Z).

    Computes Z up to the given order using right-nested commutators
    with coefficients from Arnal, Casas, Chiralt (2020).

    Works with any algebra element type that supports +, *, and a
    commutator function.

    Parameters
    ----------
    X, Y : PolyOpEx (or any type supporting +, *, simplify)
    n_order : int
        BCH expansion order (1-8).
    comm : callable, optional
        Commutator function comm(A, B).  Defaults to moyal_commutator.

    Returns
    -------
    Z : same type as X, Y
    """
    if comm is None:
        comm = moyal_commutator
    if n_order < 1 or n_order > 8:
        raise ValueError(f"n_order must be 1-8, got {n_order}")

    result = X + Y
    cache = {}

    for order in range(2, n_order + 1):
        for word, coeff in BCH_WORDS[order]:
            term = _eval_word(word, X, Y, comm, cache)
            result = result + term * coeff

    return result


# ---------------------------------------------------------------------------
# Interaction picture
# ---------------------------------------------------------------------------

def interaction_picture(H0, V, t):
    """Compute V_I(t) = e^{iH₀t/ℏ} V e^{-iH₀t/ℏ} for quadratic H₀.

    For quadratic H₀, the Heisenberg evolution of x and p is a linear
    symplectic transformation given by the classical equations of motion.
    The interaction-picture operator is obtained by substituting the
    evolved x(t), p(t) into V's Weyl symbol.

    Perturbative order is preserved from V.

    Parameters
    ----------
    H0 : PolyOpEx
        Quadratic Hamiltonian (degree ≤ 2).
    V : PolyOpEx
        Perturbation (arbitrary polynomial degree).
    t : sympy expression
        Evolution time.

    Returns
    -------
    V_I : PolyOpEx
        Interaction-picture perturbation with time-dependent coefficients.
    """
    if H0.degree > 2:
        raise ValueError("H0 must be quadratic (degree ≤ 2) for exact interaction picture")

    M_exp, v_shift = classical_evolution_matrix(H0, t)

    x_t = (M_exp[0, 0], M_exp[0, 1], v_shift[0])
    p_t = (M_exp[1, 0], M_exp[1, 1], v_shift[1])

    x_sym, p_sym = sy.symbols('_x_tmp _p_tmp', commutative=True)

    result = defaultdict(lambda: sy.Integer(0))
    result_po = {}

    for (a, b), c in V.coeffs.items():
        po = V.get_pert_order((a, b))

        x_expr = x_t[0] * x_sym + x_t[1] * p_sym + x_t[2]
        p_expr = p_t[0] * x_sym + p_t[1] * p_sym + p_t[2]

        poly_expr = sy.expand(c * x_expr**a * p_expr**b)

        poly = sy.Poly(poly_expr, x_sym, p_sym)
        for monom, coeff in poly.as_dict().items():
            result[monom] = result[monom] + coeff
            # Preserve pert_order from input monomial
            if monom in result_po:
                result_po[monom] = min(result_po[monom], po)
            else:
                result_po[monom] = po

    return PolyOpEx(dict(result), V.max_degree,
                    pert_order=result_po, max_pert_order=V.max_pert_order,
                    pert_symbol=V.pert_symbol)


# ---------------------------------------------------------------------------
# Classical phase-space evolution for quadratic Hamiltonians
# ---------------------------------------------------------------------------

def classical_evolution_matrix(H0, t):
    """Compute the classical phase-space evolution matrix for a quadratic H₀.

    For H₀ with Weyl symbol a·p² + b·p + 2c·xp + d·x + e + f·x²,
    Hamilton's equations are:

        dx/dt = 2a·p + b + 2c·x
        dp/dt = -(d + 2c·p + 2f·x)

    i.e.  d/dt [x; p] = M [x; p] + v

    where M = [[2c, 2a], [-2f, -2c]], v = [b, -d].

    Returns
    -------
    (M_exp, v_shift) where:
        M_exp : 2x2 sympy Matrix, exp(M·t)
        v_shift : 2x1 sympy Matrix, displacement
        such that [x(t); p(t)] = M_exp · [x(0); p(0)] + v_shift
    """
    if H0.degree > 2:
        raise ValueError("H0 must be quadratic")

    a = H0.coeffs.get((0, 2), 0)   # p²
    b = H0.coeffs.get((0, 1), 0)   # p
    c = H0.coeffs.get((1, 1), 0)   # xp (Weyl)
    d = H0.coeffs.get((1, 0), 0)   # x
    f = H0.coeffs.get((2, 0), 0)   # x²

    M = sy.Matrix([[c, 2*a], [-2*f, -c]])
    v = sy.Matrix([b, -d])

    M_exp = sy.simplify(sy.exp(M * t))

    if M.det() != 0:
        v_shift = sy.simplify((M_exp - sy.eye(2)) * M.inv() * v)
    else:
        v_shift = sy.zeros(2, 1)
        Mk = sy.eye(2)
        for k in range(20):
            v_shift += Mk * v * t**(k+1) / factorial(k+1)
            Mk = Mk * M
            if sy.simplify(Mk) == sy.zeros(2):
                break
        v_shift = sy.simplify(v_shift)

    return M_exp, v_shift
