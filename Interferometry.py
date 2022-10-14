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
            Oup=BCH4(Oup,u.H*(-I/hbar*u.time)).simplify()
            
        for u in self.Low:
            Olow=BCH4(u.H*(I/hbar*u.time), Olow).simplify()
            
            
            OpEx_res=BCH4(Olow, Oup).simplify()
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
    E2=C(X,Y)
    E3=C(X,E2)
    E4=C(Y,E2)
    E5=C(Y,E3)
    return E1+E2*sy.Rational('1/2')+(E3+E4*(-1))*sy.Rational('1/12')+E5*sy.Rational('-1/24')


def BCH6(X,Y):
    '''Computes the BCH to 6th order for two operator expression objects X and Y
        exp(X)*exp(Y)=exp(Z(X,Y)).
        Input: X,Y: Operator expression
        Output: Operator expression for Z(X,Y)'''
    E1=X+Y
    E2=C(X,Y)
    E3=C(X,E2)
    E4=C(Y,E2)
    E5=C(Y,E3)
    E6=C(Y,E3)






