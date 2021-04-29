from sympy import *
import numpy as np
import pandas as pd
from sympy.parsing.latex import parse_latex
class integracion_numerica():
    """
        Aproximación de integrales simples y dobles por los métodos de:
        
        Trapezoidal
            .trapezoidal()
            .trapezoidal_compuesto(particiones)
            .trapecio_compuesto_doble(intervalo2, particiones)
            
        Simpson 1/3
            .simpson1_3()
            .simpson1_3_compuesto(particiones)
            .simpson1_3_compuesto_doble(intervalo2, particiones)
            
        Simpson 3/8
            .simpson3_8()
            .simpson3_8_compuesto(particiones)
            .simpson3_4_compuesto_doble(intervalo2, particiones)
            
        OBS: Cada vez que se ejecute un método nuevo, se tiene que reinstanciar el objeto.
            
        Parámetros
        -----------------------
        limites: list 
            Lista con dos elementos representando los limites de la integral. Si es una
            integral doble, entonces son los limites de la integral de adentro.
        funcion_texto: str
            Representa la función escrita con los operadores de Python.
        
        Atributos
        -----------------------
        a: float
            Limite a de la integral. 
        b: float
            Limite b de la integral. 
        c: float
             Cuando hay dos integrales, representa el límite a de la integral
             de afuera.
        d: float
            Cuando hay dos integrales, representa el límite b de la integral
             de afuera.
             
        solucion: sympy.core.numbers.Float
            Aproximacion a la integral.
        exp: sympy.parse_expr
            Objeto que representa la función simbólica.
        metodo: str
            Nombre del método utilizado para aproximar la solución. Se actualiza
            cada vez que se cambia el método.
        
        aproximado: sympy.core.numbers.Float
            Error aproximado
        estimado: sympy.core.numbers.Float
            Error estimado (es el mismo que cota)
        cota: sympy.core.numbers.Float
            Cota del error (es el mismo que el estimado)
        total: sympy.core.numbers.Float
            Error total
        relativo: sympy.core.numbers.Float
            Error relativo
        verdadero: sympy.core.numbers.Float
            Error verdero
        
        """
    def __init__(self, limites, funcion_texto):
        
        self.a = limites[0]
        self.b = limites[1]
        self.c = None
        self.d = None
       
        self.solucion = None
        self.exp = parse_expr(str(parse_latex(funcion_texto)))
        #self.h = (self.b-self.a)/2
        self.metodo = None
        
        self.aproximado = None
        self.estimado = None
        self.cota = None 
        self.total = None
        self.relativo = None
        self.verdadero = None
        
        
        
    ####----- SIMPLES: ------####  
    def trapezoidal(self):
        x = symbols('x')
        self.solucion =  ((self.b-self.a)/2)*(self.exp.subs(x, self.a) + self.exp.subs(x, self.b))
        self.metodo = "Trapezoidal"

        segunda = integrate(diff(self.exp, x, x), (x, self.a, self.b))/(self.b-self.a)
        self.aproximado = (-((self.b - self.a)**3)/12)*segunda 
        self.estimado =  ((self.b-self.a)**3/12)*self.maximo(3, segunda)
        self.total = (-((self.b - self.a)**3)/12)*diff(self.exp, x, x).subs(x, (self.b-self.a)//2)
        return N(self.solucion)
    
    def simpson1_3(self):
        x = symbols('x')
        h = (self.b-self.a)/2
        self.solucion = h/3 * (self.exp.subs(x,self.a) + 4*self.exp.subs(x,(self.a+self.b)/2) + self.exp.subs(x,self.b) )
        self.metodo = "Simpson 1/3"
        
        tprima = diff(func, x, x, x)
        cuatriprima = diff(func, x, x, x, x)
        
        self.aproximado =  integrate(cuatriprima, (x, self.a,self.b)) 
        self.cota =  ((self.b-self.a)**3/12)*self.maximo(4,tprima)
        
        return N(self.solucion)
    
    def simpson3_8(self):
        h = (self.b-self.a)/3
        x = symbols('x')
        x0 = self.a
        x1 = x0+h
        x2 = x1+h
        x3 = self.b
        
        self.solucion = (x3-x1)*((self.exp.subs(x,x0)+3*self.exp.subs(x,x1)+3*self.exp.subs(x,x2)+self.exp.subs(x,x3)))/8
        self.metodo = "Simpson 3/8"
        
        triprima = diff(fun, x, x, x)
        cuatriprima = diff(fun, x, x, x, x)
        
        self.aproximado = -1*(h**5/90)*((integrate(cuatriprima, (x, x0, x2)))/(x3-x0))
                    
        self.cota = abs((((x3-x0)**5)/6480)*self.maximo(4, triprima))
        
        return self.solucion
    
    ####----- COMPUESTOS: ------####
    
    def trapezoidal_compuesto(self, particiones):
        x = symbols('x')
        h = (self.b-self.a)/particiones
        puntos_soporte = np.linspace(self.a + h,self.b - h, particiones-1 )
        
        funcion_evaluada = [self.exp.subs(x,t) for t in puntos_soporte]

        aprox_inter = sum(funcion_evaluada)

        self.solucion = h*((1/2)*(self.exp.subs(x,self.a) + self.exp.subs(x,self.b)) + aprox_inter) 
        self.metodo = "Trapezoidal compuesto"
        
        f_biprima = diff(self.exp,x,x)
        
        self.total = abs((-((self.b-self.a)*h**2)/12)*f_biprima.subs(x,(self.b-self.a)/2))
        self.aproximado = (-(h**2)/12)*integrate(f_biprima, (x, self.a, self.b))
        self.cota = (((self.b-self.a)*h**2)/12)*self.maximo(3, f_biprima)
        
        
        return N(self.solucion)
    
    def simpson1_3_compuesto(self, particiones):
        x = symbols('x')
        h = (self.b-self.a)/(2*particiones)
        soportes = np.linspace(self.a + h, self.b-h, 2*particiones - 1)
        
        cuatri = diff(self.exp, x, x, x, x)
        Rt = - ((h**5)/90)*cuatri.subs(x,((self.b-self.a)/2))
       
        S1 = sum([self.exp.subs(x,soportes[i]) for i in range(0,2*particiones,2)])
        S2 = sum([self.exp.subs(x,soportes[i]) for i in range(1,2*particiones-1,2)])
        
        self.solucion = (h/3)*(self.exp.subs(x, self.a) + 4*S1 + 2*S2 + self.exp.subs(x, self.b)) + (Rt if particiones > 1 else 0)
        self.metodo = "Simpson 1/3 compuesto"
        
        self.total = -(((self.b-self.a)**5)/(180*particiones**4))*cuatri.subs(x,(self.b-self.a)/2) 
        self.aproximado = -((h**4)/180)*(integrate(cuatri,(x, self.a, self.b)))
        self.cota = ((self.b-self.a)*h**4)/180*self.maximo(5,cuatri)
        
        return N(self.solucion)
    
    def simpson3_8_compuesto(self, particiones):
        x = symbols('x')
        h = (self.b-self.a)/(3*particiones)
        soportes = np.linspace(self.a+h, self.b-h, 3*particiones - 1)
        
        
        S1 = sum([self.exp.subs(x,soportes[i]) for i in range(0,3*particiones,3)])# 1,4,7,10
        S2 = sum([self.exp.subs(x,soportes[i]) for i in range(1,3*particiones,3)])# 2,5,8,11
        S3 = sum([self.exp.subs(x,soportes[i]) for i in range(2,3*particiones-1,3)])# 3,6,9
        #S4 = sum([self.exp.subs(x,soportes[i]) for i in range(3,3*particiones,3)]) 4,
        
        #print(f'{self.exp.subs(x,self.a)} {self.exp.subs(x,self.b)} 3*{[self.exp.subs(x,soportes[i]) for i in range(0,3*particiones,3)]} 3*{[self.exp.subs(x,soportes[i]) for i in range(1,3*particiones,3)]}  2*{[self.exp.subs(x,soportes[i]) for i in range(2,3*particiones-1,3)]}' )
        
        self.solucion = (3*h/8)*(self.exp.subs(x,self.a) + 3*S1 + 3*S2 + 2*S3  + self.exp.subs(x,self.b))
        self.metodo = "Simpson 3/8 compuesto"
        
        cuatriprima = diff(self.exp, x, x, x, x)
        
        self.total = (-((self.b-self.a)/80)*h**4)*cuatriprima.subs(x,(self.b-self.a)/2)
        self.aproximado = (-(h**4/80))*(integrate(cuatriprima, (x, self.a, self.b)))
        self.cota = ((self.b-self.a)*h**4)/80*self.maximo(5,cuatriprima)
        
        return  self.solucion
    
    ####----- DOBLES: ------####
    def trapecio_compuesto_doble(self, intervalo2, particiones):
        """
        Recibe los limites de la integral de afuera.
        """
        x = symbols('x')
        y = symbols('y')
        try: 
            a = float(self.a)
            b = float(self.b)
        except:
            bandera = 0
        else:
            bandera = 1
        
        if bandera:
            h = (self.b-self.a)/particiones

            puntos_soporte = np.linspace(self.a + h,self.b - h, particiones-1)

            funcion_evaluada = [self.exp.subs(x,t) for t in puntos_soporte]

            aprox_inter = sum(funcion_evaluada)
            primera_integral =  h*((1/2)*(self.exp.subs(x,self.a) + self.exp.subs(x,self.b)) + aprox_inter) 

            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/particiones
            self.c = c
            self.d = d
            puntos_soporte = np.linspace(c + h, d - h, particiones-1)
            funcion_evaluada = [primera_integral.subs(y,t) for t in puntos_soporte]

            aprox_inter = sum(funcion_evaluada)
            self.solucion =  h*((1/2)*(primera_integral.subs(y,c) + primera_integral.subs(y,d)) + aprox_inter) 

            self.metodo = "Trapezoidal compuesto doble numérico"

            f_biprima = diff(self.exp,x,x)

            self.total = abs((-((self.b-self.a)*h**2)/12)*f_biprima.subs(x,(self.b-self.a)/2))
            self.aproximado = (-(h**2)/12)*integrate(f_biprima, (x, self.a, self.b))
            self.cota = (((self.b-self.a)*h**2)/12)*self.maximo(3, f_biprima)

            return N(self.solucion)
        else:
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/particiones
            
            a = parse_expr(self.a)
            b = parse_expr(self.b)
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            puntos_soporte = np.linspace(c + h,d - h, particiones-1)
            
            funcion_evaluada = [integracion_numerica([float(a.subs(x, t)), float(b.subs(x, t))], str(self.exp.subs(x,t).subs(y, x))).trapezoidal_compuesto(particiones) for t in puntos_soporte]
            aprox_inter = sum(funcion_evaluada)
            
            self.solucion =  h*((1/2)*(integracion_numerica([float(a.subs(x, c)), float(b.subs(x,c))], str(self.exp.subs(x,c).subs(y, x))).trapezoidal_compuesto(particiones)  + integracion_numerica([float(a.subs(x, d)), float(b.subs(x, d))], str(self.exp.subs(x,d).subs(y, x))).trapezoidal_compuesto(particiones))  + aprox_inter) 
            self.metodo = "Trapezoidal compuesto doble"
            return N(self.solucion)
        
    def simpson1_3_compuesto_doble(self, intervalo2, particiones):
        """
        Recibe los limites de la integral de afuera.
        """
        x = symbols('x')
        y = symbols('y')

        try: 
            a = float(self.a)
            b = float(self.b)
        except:
            bandera = 0
        else:
            bandera = 1
        
        if bandera:
            h = (self.b-self.a)/(2*particiones)

            soportes = np.linspace(self.a + h, self.b-h, 2*particiones - 1)
            
            S1 = sum([self.exp.subs(x,soportes[i]) for i in range(0,2*particiones,2)])
            S2 = sum([self.exp.subs(x,soportes[i]) for i in range(1,2*particiones-1,2)])

            primera_integral =  (h/3)*(self.exp.subs(x, self.a) + 4*S1 + 2*S2 + self.exp.subs(x, self.b)) #+Rt
            
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(2*particiones)
            self.c = c
            self.d = d
            soportes = np.linspace(c + h, d - h, 2*particiones - 1)

            S1 = sum([primera_integral.subs(y,soportes[i]) for i in range(0,2*particiones,2)])
            S2 = sum([primera_integral.subs(y,soportes[i]) for i in range(1,2*particiones-1,2)])
            
            self.solucion =  (h/3)*(primera_integral.subs(y, c) + 4*S1 + 2*S2 + primera_integral.subs(y, d)) #+Rt

            self.metodo = "Simpson 1/3 compuesto doble numérico"

            return N(self.solucion)
            
        else:
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(2*particiones)
            
            a = parse_expr(self.a)
            b = parse_expr(self.b)
            
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            
            soportes = np.linspace(c + h, d - h, 2*particiones - 1)
            
            S1 = sum([integracion_numerica([float(a.subs(x, soportes[t])), float(b.subs(x, soportes[t]))], str(self.exp.subs(x,soportes[t]).subs(y, x))).simpson1_3_compuesto(particiones) for t in range(0,2*particiones,2)])
            S2 = sum([integracion_numerica([float(a.subs(x, soportes[t])), float(b.subs(x, soportes[t]))], str(self.exp.subs(x,soportes[t]).subs(y, x))).simpson1_3_compuesto(particiones) for t in range(1,2*particiones-1,2)])
            
            
            self.solucion =  (h/3)*(integracion_numerica([float(a.subs(x, c)), float(b.subs(x,c))], str(self.exp.subs(x,c).subs(y, x))).simpson1_3_compuesto(particiones) + 4*S1 + 2*S2 + integracion_numerica([float(a.subs(x, d)), float(b.subs(x, d))], str(self.exp.subs(x,d).subs(y, x))).simpson1_3_compuesto(particiones) ) #+Rt
            
            self.metodo = "Simpson 1/3 compuesto doble"
            return N(self.solucion)            
        
    def simpson3_8_compuesto_doble(self, intervalo2, particiones):
        """
        Recibe los limites de la integral de afuera.
        """
        x = symbols('x')
        y = symbols('y')
        try: 
            a = float(self.a)
            b = float(self.b)
        except:
            bandera = 0
        else:
            bandera = 1
        
        if bandera:
           
            h = (self.b-self.a)/(3*particiones)

            soportes = np.linspace(self.a+h, self.b-h, 3*particiones - 1)
            
            S1 = sum([self.exp.subs(x,soportes[i]) for i in range(0,3*particiones,3)])# 1,4,7,10
            S2 = sum([self.exp.subs(x,soportes[i]) for i in range(1,3*particiones,3)])# 2,5,8,11
            S3 = sum([self.exp.subs(x,soportes[i]) for i in range(2,3*particiones-1,3)])# 3,6,9
            
            primera_integral =  (3*h/8)*(self.exp.subs(x,self.a) + 3*S1 + 3*S2 + 2*S3  + self.exp.subs(x,self.b))
            
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(3*particiones)
            self.c = c
            self.d = d
            soportes = np.linspace(c + h, d - h, 3*particiones - 1)

            S1 = sum([primera_integral.subs(y,soportes[i]) for i in range(0,3*particiones,3)])
            S2 = sum([primera_integral.subs(y,soportes[i]) for i in range(1,3*particiones-1,3)])
            S3 = sum([primera_integral.subs(y,soportes[i]) for i in range(2,3*particiones-1,3)])
            
            self.solucion =  (3*h/8)*(primera_integral.subs(y,c) + 3*S1 + 3*S2 + 2*S3  + primera_integral.subs(y,d))

            self.metodo = "Simpson 3/8 compuesto doble numérico"

            return N(self.solucion)
            
        else:
            c = intervalo2[0]
            d = intervalo2[1]
            h = (d-c)/(3*particiones)
            
            a = parse_expr(self.a)
            b = parse_expr(self.b)
            
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            
            soportes = np.linspace(c + h, d - h, 3*particiones - 1)
            
            S1 = sum([integracion_numerica([float(a.subs(x, soportes[t])), float(b.subs(x, soportes[t]))], str(self.exp.subs(x,soportes[t]).subs(y, x))).simpson3_8_compuesto(particiones) for t in range(0,3*particiones,3)])
            S2 = sum([integracion_numerica([float(a.subs(x, soportes[t])), float(b.subs(x, soportes[t]))], str(self.exp.subs(x,soportes[t]).subs(y, x))).simpson3_8_compuesto(particiones) for t in range(1,3*particiones,3)])
            S3 = sum([integracion_numerica([float(a.subs(x, soportes[t])), float(b.subs(x, soportes[t]))], str(self.exp.subs(x,soportes[t]).subs(y, x))).simpson3_8_compuesto(particiones) for t in range(3,3*particiones-1,3)])
            
            self.solucion =  (3*h/8)*(integracion_numerica([float(a.subs(x, c)), float(b.subs(x,c))], str(self.exp.subs(x,c).subs(y, x))).simpson3_8_compuesto(particiones) + 3*S1 + 3*S2 + 2*S3 + integracion_numerica([float(a.subs(x, d)), float(b.subs(x, d))], str(self.exp.subs(x,d).subs(y, x))).simpson3_8_compuesto(particiones) ) 
            
            self.metodo = "Simpson 3/8 compuesto doble"
            return N(self.solucion)
        
    ####----- ERRORES: ------####
    def maximo(self, grado, f = None):
        if not f:
            f = self.exp
        x = symbols = 'x'
        g = self.exp
        for _ in range(grado):
            g = diff(g, x)

        try:             
            puntos_criticos = solve(g,x)
        except:
            puntos_criticos = []
        puntos_criticos = [p for p in puntos_criticos if p >= self.a and p<=self.b] + [self.a,self.b]
        
        maximo = float(abs(f.subs(x,puntos_criticos[0])))
        for punto in puntos_criticos:
            valor = float(abs(f.subs(x, N(punto))))
            if valor > maximo:
                maximo = valor
                
        return maximo
        
    def errores(self, error = None):
        x = symbols('x')
        y = symbols('y')
        if error:
            if error == 'total':
                return self.total
            elif error == 'verdadero':
                return integrate(self.exp, (x, self.a, self.b)) - self.solucion
            elif error == 'relativo ':
                return (1 - self.solucion/integrate(self.exp, (x, self.a, self.b)))*100
            elif error == 'aproximado':
                return self.aproximado            
            elif error == 'estimado':
                return self.estimado
            elif error == 'cota':
                return self.cota
            elif error == 'total':
                return self.total
            else:
                raise ValueError('No existe ese error')
        else: 
            if 'doble' in self.metodo:
                if 'numérico' not in self.metodo:
                    integral_y = integrate(self.exp, (y, self.a, self.b))
                    valor_verdadero = N(integrate (integral_y, (x, self.c, self.d)))
                else:
                    integral_x = integrate(self.exp, (x, self.a, self.b))
                    valor_verdadero = N(integrate (integral_x, (y, self.c, self.d)))
            else:
                valor_verdadero = integrate(self.exp, (x, self.a, self.b))
                
            self.verdadero = valor_verdadero - self.solucion
            self.relativo = (1 - self.solucion/valor_verdadero)*100

            error = {'Total': self.total, 'Verdadero': self.verdadero, 'Relativo': self.relativo, 'Aproximado':self.aproximado, 
                    'Estimado':self.estimado, 'Cota':self.cota}
            return pd.DataFrame(error.values(),index = error.keys(), columns = ['Valor']).reset_index().rename({'index':'Error'}, axis = 1).set_index('Error')