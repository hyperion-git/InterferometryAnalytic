import sympy as sy

try:
    from IPython.display import display as _ipython_display
except ImportError:
    _ipython_display = None

from poly_operator import PolyOpEx, moyal_commutator, bchn, hbar


# Operator expression class  -----------------------------------------------

class OpEx():
    def __init__(self, param):
        '''Class is basic element to implement operator algebra.
        An instance of OpEx consists of a list [a,b,c,d,e,f] representing the expression
        ap^2 + b*p +c(px+xp) + d*x + e + f*x^2.

        Input: list [a,b,c,d,e,f] representing the expression, parameters are sympy or numpy variables
        '''
        if not len(param)==6:
            raise ValueError('length of parameter list must be 6')
        self.a=param[0]
        self.b=param[1]
        self.c=param[2]
        self.d=param[3]
        self.e=param[4]
        self.f=param[5]

    def __repr__(self):
        if _ipython_display is not None:
            p=sy.symbols('\hat{p}')
            px_xp=sy.symbols('(\hat{p}\hat{x}+\hat{x}\hat{p})')
            x=sy.symbols('\hat{x}')
            qid=sy.symbols('\hat{1}')
            _ipython_display(self.a*p**2+self.b*p+self.c*px_xp+self.d*x+self.f*x**2+self.e*qid)
            return ''
        return (f'OpEx([{self.a}, {self.b}, {self.c}, '
                f'{self.d}, {self.e}, {self.f}])')

    def __add__(self, other):
        #allows to add two operator expressions
         return OpEx([self.a+other.a,  self.b+other.b, self.c+other.c, self.d+other.d,self.e+other.e,self.f+other.f] )

    def __mul__(self, other):
        # allows to multiply an operator expression from the right with a numpy or sympy variable
        return OpEx([self.a*other,self.b*other,self.c*other,self.d*other,self.e*other,self.f*other] )

    def simplify(self):
        # simplifies the expression for each parameter using sympy.simplify()
        return OpEx([sy.simplify(self.a),sy.simplify(self.b),sy.simplify(self.c), sy.simplify(self.d), sy.simplify(self.e), sy.simplify(self.f)] )

    def to_poly(self, **kwargs):
        """Convert to PolyOpEx for use with the generalized algebra."""
        return PolyOpEx.from_opex(self, **kwargs)


# Interferometer classes----------------------

class Hamiltonian(OpEx):
    '''Class to represent a Hamiltonian.
    Input: list [a,b,c,d,e,f] representing the Hamiltonian H=ap^2 + b*p +c(px+xp) + d*x + e + f*x^2.,
    parameters are sympy or numpy variables'''
    def __init__(self, param):
        super().__init__(param)


class Pulse():
    def __init__(self, wave_vector):
        '''Class to represent a momentum kick
        Input: wave_vector: sympy expression, representing the operator exp(i*wave_vector*x)'''
        hbar_sym=sy.symbols('hbar')
        self.k=wave_vector
        self.H=OpEx([0,0,0,-hbar_sym*wave_vector,0,0])
        self.time=1 # gives exp(i*k*z) in time-evolution operator
    def __repr__(self):
        from sympy import I
        if _ipython_display is not None:
            x=sy.symbols('\hat{x}')
            print('implements:')
            _ipython_display(sy.exp(I*self.k*x))
            return ''
        return f'Pulse(k={self.k})'

class U():
    '''Class representing a time evolution operator.
    Input: Hamiltonian: A Hamiltonian object
           time: a sympy variable for the evolution time.
           This object represents the operator exp(-I/hbar*H*time)'''
    def __init__(self,Hamiltonian, time):
        '''Hamiltonian: Operator Expression with real! parameters'''
        self.H=Hamiltonian
        self.time=time
    def __repr__(self):
        from sympy import I
        if _ipython_display is not None:
            hbar_sym=sy.symbols('hbar')
            H=sy.symbols('\hat{H}')
            print('implements:')
            _ipython_display(sy.exp(-I*H*self.time/hbar_sym))
            return ''
        return f'U(H, t={self.time})'


class Interferometer():
    '''Class for an interferometer.

    Input: UpperSequence, LowerSequence: list of U and Pulse objects.
    BCHOrder: order of BCH expansion (1-8, default 8).
    use_poly: if True, use the generalized PolyOpEx algebra (supports
              non-quadratic Hamiltonians and perturbative tracking).
    max_degree: if set with use_poly, truncate polynomial degree.
    max_pert_order: if set with use_poly, truncate perturbative order.
    '''
    def __init__(self, UpperSequence, LowerSequence, BCHOrder=8,
                 use_poly=False, max_degree=None, max_pert_order=None):
        self.Upper=UpperSequence
        self.Lower=LowerSequence
        self.UpperSequenceLength=len(UpperSequence)
        self.LowerSequenceLength=len(LowerSequence)
        self.BCHExpansionOrder=BCHOrder
        self.use_poly = use_poly
        self.max_degree = max_degree
        self.max_pert_order = max_pert_order

    def overlap(self):
        '''Calculates the phase of an interferometer object.

        Returns (res_dic, result) where:
          - res_dic: dict with keys 'p2', 'p', 'px_xp', 'x', 'const', 'x2'
          - result: OpEx (if use_poly=False) or PolyOpEx (if use_poly=True)
        '''
        if self.use_poly:
            return self._overlap_poly()
        return self._overlap_opex()

    def overlap_graded(self):
        '''Compute the overlap and return its graded decomposition.

        Returns (grades, result) where:
          - grades: dict {n: PolyOpEx} with grade-n component of the overlap
            (0=phase, 1=displacement, 2=distortion, 3+=aberrations)
          - result: the full PolyOpEx overlap operator

        Always uses the PolyOpEx pipeline internally. The OpEx result
        is promoted to PolyOpEx if needed.
        '''
        if self.use_poly:
            _, result = self._overlap_poly()
        else:
            _, opex_result = self._overlap_opex()
            result = PolyOpEx.from_opex(opex_result)
        return result.by_grade(), result

    def _overlap_opex(self):
        """Original quadratic overlap using OpEx + BCHN."""
        hbar_sym = sy.symbols('hbar')
        from sympy import I
        Oup = OpEx([0,0,0,0,0,0])
        Olow = OpEx([0,0,0,0,0,0])

        for u in self.Upper:
            Oup = BCHN(Oup, u.H*(-I/hbar_sym*u.time), self.BCHExpansionOrder).simplify()

        for u in self.Lower:
            Olow = BCHN(u.H*(I/hbar_sym*u.time), Olow, self.BCHExpansionOrder).simplify()

        OpEx_res = BCHN(Olow, Oup, self.BCHExpansionOrder).simplify()
        res_dic = { 'p2':    OpEx_res.a,
                    'p':     OpEx_res.b,
                    'px_xp': OpEx_res.c,
                    'x':     OpEx_res.d,
                    'const': OpEx_res.e,
                    'x2':    OpEx_res.f}
        return (res_dic, OpEx_res)

    def _overlap_poly(self):
        """Generalized overlap using PolyOpEx + bchn (Moyal bracket)."""
        from sympy import I

        md = self.max_degree
        mpo = self.max_pert_order

        def _to_poly(u):
            """Convert a sequence element's H to PolyOpEx."""
            if isinstance(u.H, PolyOpEx):
                return u.H
            return PolyOpEx.from_opex(u.H, max_degree=md, max_pert_order=mpo)

        Oup = PolyOpEx.zero(max_degree=md, max_pert_order=mpo)
        for u in self.Upper:
            H_poly = _to_poly(u)
            arg = H_poly * (-I / hbar * u.time)
            Oup = bchn(Oup, arg, self.BCHExpansionOrder).simplify()

        Olow = PolyOpEx.zero(max_degree=md, max_pert_order=mpo)
        for u in self.Lower:
            H_poly = _to_poly(u)
            arg = H_poly * (I / hbar * u.time)
            Olow = bchn(arg, Olow, self.BCHExpansionOrder).simplify()

        result = bchn(Olow, Oup, self.BCHExpansionOrder).simplify()

        # Build result dict (extract quadratic terms if present)
        res_dic = { 'p2':    result.coeffs.get((0, 2), 0),
                    'p':     result.coeffs.get((0, 1), 0),
                    'px_xp': result.coeffs.get((1, 1), 0) / 2 if (1, 1) in result.coeffs else 0,
                    'x':     result.coeffs.get((1, 0), 0),
                    'const': result.coeffs.get((0, 0), 0),
                    'x2':    result.coeffs.get((2, 0), 0)}
        return (res_dic, result)



# Helper Functions  --------------------------------------------

def C(OpEx1, OpEx2):
    """input: OpEx1,OpEx2: Of two quadratic operator expressions
       output: Operator expression for the commutator C(OpEx1,OpEx2) """
    from sympy import I
    hbar_sym=sy.symbols('hbar')
    a1=OpEx1.a
    b1=OpEx1.b
    c1=OpEx1.c
    d1=OpEx1.d
    e1=OpEx1.e
    f1=OpEx1.f

    a2=OpEx2.a
    b2=OpEx2.b
    c2=OpEx2.c
    d2=OpEx2.d
    e2=OpEx2.e
    f2=OpEx2.f

    aI= 4*I*hbar_sym*(a2*c1-a1*c2)                  # p**p
    bI= 2*I*hbar_sym*(b2*c1-b1*c2+a2*d1-a1*d2)      # p
    cI= 2*I*hbar_sym*(a2*f1-a1*f2)                  # px_xp
    dI= 2*I*hbar_sym*(c2*d1-c1*d2+b2*f1-f2*b1)      # x
    eI= I*hbar_sym*(b2*d1-b1*d2)                    # 1
    fI= 4*I*hbar_sym*(f1*c2-c1*f2)                  # x**x
    return OpEx([aI,bI,cI,dI,eI,fI])



def BCHN(X, Y, nOrder=8):
    '''Computes the BCH up to 8th order for two operator expression objects X and Y
        exp(X)*exp(Y)=exp(Z(X,Y)).
        Input: X,Y: Operator expression; nOrder optional desired BCH order
        Output: Operator expression for Z(X,Y)
        Implementation: - Recursive right-nested commutators with maximal reuse
                        - Minimally right-nested operator basis
                        - Details in https://link.springer.com/article/10.1007/s00009-020-01681-6'''

    if nOrder >= 1:
        # do nothing here
        E1=X+Y
        phi1=E1
        phiTemp=phi1
    if nOrder >= 2:
        E2xy=C(X,Y)
        phi2=E2xy*sy.Rational('1/2')
        phiTemp=phiTemp+phi2
    if nOrder >= 3:
        E3xxy=C(X,E2xy)
        E3yxy=C(Y,E2xy)

        phi3=(E3xxy+E3yxy*(-1))*sy.Rational('1/12')
        phiTemp=phiTemp+phi3

    if nOrder >= 4:
        E4xyxy=C(X,E3yxy)

        phi4=E4xyxy*sy.Rational('-1/24')
        phiTemp=phiTemp+phi4

    if nOrder >= 5:
        E5xxxxy=C(X,C(X,E3xxy))
        E5xyxxy=C(X,C(Y,E3xxy))
        E5xyyxy=C(X,C(Y,E3yxy))
        E5yxxxy=C(Y,C(X,E3xxy))
        E5yyxxy=C(Y,C(Y,E3xxy))
        E5yyyxy=C(Y,C(Y,E3yxy))

        phi5a=(E5xyyxy*(-1)+E5yxxxy)*sy.Rational('1/360')
        phi5b=(E5xyxxy*(-1)+E5yyxxy)*sy.Rational('1/120')
        phi5c=(E5xxxxy*(-1)+E5yyyxy)*sy.Rational('1/720')

        phi5=phi5a+phi5b+phi5c
        phiTemp=phiTemp+phi5

    if nOrder >= 6:
        E6xxyyxy=C(X,E5xyyxy)
        E6xyyxxy=C(X,E5yyxxy)
        E6xyyyxy=C(X,E5yyyxy)
        E6yxxxxy=C(Y,E5xxxxy)

        phi6a=E6xxyyxy*sy.Rational('-1/720')
        phi6b=E6xyyxxy*sy.Rational('1/240')
        phi6c=(E6xyyyxy+E6yxxxxy)*sy.Rational('1/1440')

        phi6=phi6a+phi6b+phi6c
        phiTemp=phiTemp+phi6

    if nOrder >=7:
        # Additional order 6 elements needed for order 7
        E6xxxxxy=C(X,E5xxxxy)
        E6xyxxxy=C(X,E5yxxxy)

        E7xyxyxxy=C(X,C(Y,E5xyxxy))   # 1/1008
        E7yxyxxxy=C(Y,E6xyxxxy)       #-1/1260
        E7yxyyxxy=C(Y,E6xyyxxy)       #-1/1680
        E7yyxyyxy=C(Y,C(Y,E5xyyxy))   #-1/2520

        E7xyyyxxy=C(X,C(Y,E5yyxxy))   # 1/3360
        E7yyxxxxy=C(Y,E6yxxxxy)       # 1/3360
        E7yyxyxxy=C(Y,C(Y,E5xyxxy))   #-1/3360

        E7xxyxxxy=C(X,E6xyxxxy)       # 1/5040
        E7xyxyyxy=C(X,C(Y,E5xyyxy))   # 1/5040

        E7xyyxxxy=C(X,C(Y,E5yxxxy))   #-1/7560
        E7yyyxxxy=C(Y,C(Y,E5yxxxy))   # 1/7560

        E7xyyyyxy=C(X,C(Y,E5yyyxy)) #  1/10080
        E7yyyyxxy=C(Y,C(Y,E5yyxxy)) #  1/10080
        E7xyxxxxy=C(X,E6yxxxxy)     #  1/10080
        E7yxxxxxy=C(Y,E6xxxxxy)     # -1/10080
        E7xxyyxxy=C(X,E6xyyxxy)     # -1/10080

        E7xxxxxxy=C(X,E6xxxxxy)     #  1/30240
        E7yyyyyxy=C(Y,C(Y,E5yyyxy)) # -1/30240

        phi7a=E7xyxyxxy*sy.Rational('1/1008')+E7yxyxxxy*sy.Rational('-1/1260')+E7yxyyxxy*sy.Rational('-1/1680')+E7yyxyyxy*sy.Rational('-1/2520')
        phi7b=(E7xyyyxxy+E7yyxxxxy+E7yyxyxxy*(-1))*sy.Rational('1/3360')
        phi7c=(E7xxyxxxy+E7xyxyyxy)*sy.Rational('1/5040')
        phi7d=(E7xyyxxxy*(-1)+E7yyyxxxy)*sy.Rational('1/7560')
        phi7e=(E7xyyyyxy+E7yyyyxxy+E7xyxxxxy+(E7yxxxxxy+E7xxyyxxy)*(-1))*sy.Rational('1/10080')
        phi7f=(E7xxxxxxy+E7yyyyyxy*(-1))*sy.Rational('1/30240')

        phi7=phi7a+phi7b+phi7c+phi7d+phi7e+phi7f
        phiTemp=phiTemp+phi7

        if nOrder >=8:
            E6xxyxxy=C(X,E5xyxxy)
            E7xxxyxxy=C(X,E6xxyxxy)

            # Brute force implementation again
            E8xxxyyyxy=C(X,C(X,E6xyyyxy))           # -5/24192
            E8xxyxyyxy=C(X,E7xyxyyxy)               #  1/2520
            E8xxyyyyxy=C(X,E7xyyyyxy)               #  1/20160
            E8xyxxyyxy=C(X,C(Y,E6xxyyxy))           #  1/15120
            E8xyxyyxxy=C(X,E7yxyyxxy)               # -1/2016
            E8xyxyyyxy=C(X,C(Y,E6xyyyxy))           # -1/20160
            E8xyyxyxxy=C(X,E7yyxyxxy)               #  1/20160
            E8xyyxyyxy=C(X,E7yyxyyxy)               # -1/10080
            E8xyyyyyxy=C(X,E7yyyyyxy)               # -1/60480
            E8yxxxxxxy=C(Y,E7xxxxxxy)               # -1/60480
            E8yxxxyxxy=C(Y,C(X,C(X,E5xyxxy)))       #  1/20160
            E8yxxyxxxy=C(Y,E7xxyxxxy)               # -1/5040
            E8yyxxxxxy=C(Y,E7yxxxxxy)               #  1/20160

            phi8a=E8xxxyyyxy*sy.Rational('-5/24192')
            phi8b=E8xxyxyyxy*sy.Rational('1/2520')
            phi8c=(E8xxyyyyxy+E8xyxyyyxy*(-1)+E8xyyxyxxy+E8yxxxyxxy+E8yyxxxxxy)*sy.Rational('1/20160')
            phi8d=E8xyxxyyxy*sy.Rational('1/15120')
            phi8e=E8xyxyyxxy*sy.Rational('-1/2016')
            phi8f=(E8xyyyyyxy+E8yxxxxxxy)*sy.Rational('-1/60480')
            phi8g=E8yxxyxxxy*sy.Rational('-1/5040')
            phi8h=E8xyyxyyxy*sy.Rational('-1/10080')

            phi8=phi8a+phi8b+phi8c+phi8d+phi8e+phi8f+phi8g+phi8h
            phiTemp=phiTemp+phi8

    phiFinal=phiTemp

    return phiFinal
