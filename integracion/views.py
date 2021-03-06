from django.shortcuts import render, redirect
from sympy import parse_expr, latex, integrate, symbols
from symboesfm.metodos import integracion_numerica
from sympy.parsing.latex import parse_latex

def view(request):
    
    equation =  request.session['eq'] 
    #equation = latex(parse_expr( request.session['eq']  )).replace('\\','*').replace('*',chr(92) )
    
    if request.session['tipo'] == "simple":
        integral = integracion_numerica([float(request.session['a']), float(request.session['b'])], request.session['eq'] )
        metodos = { '1': integral.trapezoidal_compuesto, 
                    '2': integral.simpson1_3_compuesto, 
                    '3': integral.simpson3_8_compuesto}

        nombres = { '1': 'Trapezoidal', 
                    '2': 'Simpson 1/3', 
                    '3': 'Simpson 3/8'}

        aproximacion = metodos[request.session['metodo']](int(request.session['particiones']))
        
        errores = integral.errores().reset_index().to_dict('records')
        pasos = integral.pasos
        return render(request, 'integracion/view.html', {'equation':latex(parse_expr(request.session['eq'])), 'aproximacion': aproximacion, 'metodo':nombres[request.session['metodo']], 'tipo':'simple', 'errores':errores, 'pasos':pasos})

    elif request.session['tipo'] == "doble":

        try:
            a = float(request.session['a'])
            b = float(request.session['b'])
           
        except:
            a =  request.session['a']
            b =  request.session['b']
            tipo = "doble2"
            aa = latex(parse_expr(  a  )).replace('\\','*').replace('*',chr(92) )
            bb = latex(parse_expr(  b  )).replace('\\','*').replace('*',chr(92) )
        else:
            tipo = "doble1"
            aa = None
            bb = None
            
        integral = integracion_numerica([a, b], request.session['eq'])
        metodos = { '1': integral.trapecio_compuesto_doble, 
                    '2': integral.simpson1_3_compuesto_doble, 
                    '3': integral.simpson3_8_compuesto_doble}

        nombres = { '1': 'Trapezoidal Doble', 
                    '2': 'Simpson 1/3 Doble', 
                    '3': 'Simpson 3/8 Doble'}

        aproximacion = metodos[request.session['metodo']]([float(request.session['c']), float(request.session['d'])], int(request.session['particiones']))
        errores = integral.errores().reset_index().to_dict('records')
        pasos = integral.pasos
        return render(request, 'integracion/view.html', {'equation':latex(parse_expr(request.session['eq'])), 'aproximacion': aproximacion, 'metodo':nombres[request.session['metodo']], 'tipo':tipo, 'aa':aa, 'bb':bb, 'errores':errores, 'pasos':pasos})

    elif request.session['tipo'] == "extrapolacion":

        nombres = { '1': 'Romberg con Trapezoidal', 
                    '2': 'Romberg con Simpson 1/3', 
                    '3': 'Romberg con Simpson 3/8'}

        integral = integracion_numerica([float(request.session['a']), float(request.session['b'])], request.session['eq'] )

        aproximacion = integral.romberg(n = int(request.session['particiones']), metodo = request.session['metodo'])
        errores = integral.errores().reset_index().to_dict('records')
        return render(request, 'integracion/view.html', {'equation':latex(parse_expr(request.session['eq'])), 'aproximacion': aproximacion, 'metodo':nombres[request.session['metodo']], 'tipo':'romberg', 'errores':errores})
    
    elif request.session['tipo'] == "indefinida":
        x = symbols('x')
        
        equation = parse_expr(request.session['eq'])
        aproximacion = integrate(equation, x)
        return render(request, 'integracion/view.html', {'equation':latex(parse_expr(request.session['eq'])), 'aproximacion': latex(aproximacion), 'metodo': None, 'tipo': 'Indeinida', 'errores': None})
def indefinida(request):
    return render(request, 'integracion/indefinida.html')
def simple(request):

    return render(request, 'integracion/simple.html')

def doble(request):

    return render(request, 'integracion/doble.html')
    
def extrapolacion(request):
    return render(request, 'integracion/extrapolacion.html')

def submit(request):
    if 'indefinida' in request.META.get('HTTP_REFERER'):
        request.session['eq'] =  str(parse_latex(request.POST['eq'].replace('\\left','').replace('\\right','').replace('e', 'E')))
        request.session['tipo'] = "indefinida"
        return redirect('view')
    if request.method == 'POST':
        
        request.session['eq'] =  str(parse_latex(request.POST['eq'].replace('\\left','').replace('\\right','').replace('e', 'E')))
        request.session['eql'] = request.POST['eq']
        request.session['a'] =  request.POST['a']
        request.session['b'] =  request.POST['b']
        request.session['metodo'] =  request.POST['metodo']
        request.session['particiones'] =  request.POST['particiones']
        request.session['tipo'] = request.POST['tipo']
        if request.session['tipo'] == "doble":
            request.session['c'] =  request.POST['c']
            request.session['d'] =  request.POST['d']
        return redirect('view')
