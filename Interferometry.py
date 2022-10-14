import sympy as sy
from IPython.display import display

# Operator expression class  -----------------------------------------------

class OpEx():
    def __init__(self, param):
        '''Class is basic element to implement operator algebra.
        An instance of OpEx consists of a list [a,b,c,d,e,f] representing the expressiopn 
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
        p=sy.symbols('\hat{p}')
        px_xp=sy.symbols('(\hat{p}\hat{x}+\hat{x}\hat{p})')
        x=sy.symbols('\hat{x}')
        display(self.a*p**2+self.b*p+self.c*px_xp+self.d*x+self.f*x**2+self.e)
        return ''
        
    def __add__(self, other):
        #allows to add two operator expressions
         return OpEx([self.a+other.a,  self.b+other.b, self.c+other.c, self.d+other.d,self.e+other.e,self.f+other.f] )
        
    def __mul__(self, other):
        # allows to multiply an operator expression from thr right with a numpy or sympy variable
        return OpEx([self.a*other,self.b*other,self.c*other,self.d*other,self.e*other,self.f*other] )

    def simplify(self):
        # simplifies the expression for each parameter using sympy.simplify()
        return OpEx([sy.simplify(self.a),sy.simplify(self.b),sy.simplify(self.c), sy.simplify(self.d), sy.simplify(self.e), sy.simplify(self.f)] )
        


# Interferometer classes----------------------

class Hamiltonian(OpEx):
    '''Class to represent a Hamiltonian.
    Input: list [a,b,c,d,e] representing the Hamiltonian H=ap^2 + b*p +c(px+xp) + d*x + e + f*x^2.,
    parameters are sympy or numpy variables'''
    def __init__(self, param):
        super().__init__(param)

        
class Pulse():
    def __init__(self, wave_vector):
        '''Class to represent a momentum kick
        Input: wave_vector: sympy expression, representing the operator exp(i*wave_vector*x)'''
        hbar=sy.symbols('hbar')
        self.k=wave_vector
        self.H=OpEx([0,0,0,-hbar*wave_vector,0,0])
        self.time=1 # gives exp(i*k*z) in time-evolution operator
    def __repr__(self):
        from sympy import I
        x=sy.symbols('\hat{x}')
        print('implements:')
        display(sy.exp(I*self.k*x))
        return ''

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
        hbar=sy.symbols('hbar')
        H=sy.symbols('\hat{H}')
        print('implements:')
        display(sy.exp(-I*H*self.time/hbar))
        return ''
        
        
class Interferometer():
    '''Class for an interferometer
    Input: UpperSequence, LowerSequence: list of Hamiltonian and Pulse objects'''
    def __init__(self, UpperSequence, LowerSequence):
        '''UpperSequence, LowerSequence: List of Pulse and U objects denoting the sequence
           ordering: In the correct way the time evolution operators are writen down for the full time evolution'''
        self.Up=UpperSequence
        self.Low=LowerSequence
        
    def phase(self):
        '''Calculates the phase of an interferometer object'''
        hbar=sy.symbols('hbar')
        from sympy import I
        Oup=OpEx([0,0,0,0,0,0])
        Olow=OpEx([0,0,0,0,0,0])
        
        for u in self.Up:
            Oup=BCHN(Oup,u.H*(-I/hbar*u.time)).simplify()
            
        for u in self.Low:
            Olow=BCHN(u.H*(I/hbar*u.time), Olow).simplify()
            
            
            OpEx_res=BCHN(Olow, Oup).simplify()
            a_res=OpEx_res.a
            b_res=OpEx_res.b
            c_res=OpEx_res.c
            d_res=OpEx_res.d
            e_res=OpEx_res.e
            f_res=OpEx_res.f
            res_dic = { 'p2':    a_res,
                        'p':     b_res,
                        'px_xp': c_res,
                        'x':     d_res,
                        'const': e_res, 
                        'x2':    f_res}                 
        return (res_dic, OpEx_res)


        
# Helper Functions  --------------------------------------------

def C(OpEx1, OpEx2):
    """input: OpEx1,OpEx2: Two operator expressions
       output: Operator expression for the commutator C(OpEx1,OpEx2) """
    from sympy import I
    hbar=sy.symbols('hbar')
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
    
    # Get someone to check this for typos
    aI= 4*I*hbar*(a2*c1-a1*c2)                  # p**p
    bI= 2*I*hbar*(b2*c1-b1*c2+a2*d1-a1*d2)      # p
    cI= 2*I*hbar*(a2*f1-a1*f2)                  # px_xp
    dI= 2*I*hbar*(c2*d1-c1*d2+b2*f1-f2*b1)      # x
    eI= I*hbar*(b2*d1-b1*d2)                    # 1
    fI= 4*I*hbar*(f1*c2-c1*f2)                  # x**x
    return OpEx([aI,bI,cI,dI,eI,fI])



def BCH4(X,Y):
    '''Computes the BCH to 4th order for two operator expression objects X and Y
        exp(X)*exp(Y)=exp(Z(X,Y)).
        Input: X,Y: Operator expression
        Output: Operator expression for Z(X,Y)'''
    E1=X+Y
    E2xy=C(X,Y)
    E3x=C(X,E2)
    E3y=C(Y,E2)
    E4=C(X,E3y)
    return E1+E2xy*sy.Rational('1/2')+(E3x+E3y*(-1))*sy.Rational('1/12')+E4*sy.Rational('-1/24')


def BCHN(X,Y,nOrder=7):
    '''Computes the BCH up to 6th order for two operator expression objects X and Y
        exp(X)*exp(Y)=exp(Z(X,Y)).
        Input: X,Y: Operator expression
        Output: Operator expression for Z(X,Y)'''
    
    # Implement this to order 10
    # https://link.springer.com/article/10.1007/s00009-020-01681-6    
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
        phiTemp=phiTemp+phi5a+phi5b+phi5c
        
    if nOrder >= 6:
        E6xxyyxy=C(X,E5xyyxy)
        E6xyyxxy=C(X,E5yyxxy)
        E6xyyyxy=C(X,E5yyyxy)
        E6yxxxxy=C(Y,E5xxxxy)
        
        phiTemp=phiTemp+E6xxyyxy*sy.Rational('-1/720')+E6xyyxxy*sy.Rational('1/240')+(E6xyyyxy+E6yxxxxy)*sy.Rational('1/1440')
   
    if nOrder >=7:
        # Order 7 needs to be checked... we need a better way to implement higher orders...
        # However, we still want to reuse commutators!
        E6xxxxxy=C(X,E5xxxxy)
        E6xyxxxy=C(X,E5yxxxy)
        
        # E7xxxxxxy=C(X,E6xxxxxy) # 1/30240 A
        # E7xxyxxxy=C(X,E6xyxxxy) # 1/5040
        # E7xxyyxxy=C(X,E6xyyxxy)     #-1/10080
       
        # E7xyxxxxy=C(X,E6yxxxxy)       #1/10080
        # E7xyxyxxy=C(X,C(Y,E5xyxxy))   #1/1008
        # E7xyxyyxy=C(X,C(Y,E5xyyxy))   #1/5040
        
        # E7xyyxxxy=C(X,C(Y,E5yxxxy)) #-1/7560
        # E7xyyyxxy=C(X,C(Y,E5yyxxy)) #1/3360
        # E7xyyyyxy=C(X,C(Y,E5yyyxy)) #1/10080
        
        # E7yxxxxxy=C(Y,E6xxxxxy) #-1/10080 A
        # E7yxyxxxy=C(Y,E6xyxxxy) #-1/1260
        # E7yxyyxxy=C(Y,E6xyyxxy) #-1/1680
        
        # E7yyxxxxy=C(Y,E6yxxxxy)     #1/3360
        # E7yyxyxxy=C(Y,C(Y,E5xyxxy)) #-1/3360
        # E7yyxyyxy=C(Y,C(Y,E5xyyxy)) #-1/2520
        
        # E7yyyxxxy=C(Y,C(Y,E5yxxxy)) #1/7560
        # E7yyyyxxy=C(Y,C(Y,E5yyxxy)) #1/10080
        # E7yyyyyxy=C(Y,C(Y,E5yyyxy)) #-1/30240

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
    
    return phiTemp





