from django.shortcuts import render, redirect
from sympy import parse_expr, latex
from symboesfm.metodos import integracion_numerica
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

        return render(request, 'integracion/view.html', {'equation':equation, 'aproximacion': aproximacion, 'metodo':nombres[request.session['metodo']], 'tipo':'simple'})

    else:

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
        return render(request, 'view.html', {'equation':equation, 'aproximacion': aproximacion, 'metodo':nombres[request.session['metodo']], 'tipo':tipo, 'aa':aa, 'bb':bb})


def simple(request):

    return render(request, 'integracion/simple.html')

def doble(request):

    return render(request, 'integracion/doble.html')

def submit(request):
    if request.method == 'POST':
        request.session['eq'] =  request.POST['eq']
        request.session['a'] =  request.POST['a']
        request.session['b'] =  request.POST['b']
        request.session['metodo'] =  request.POST['metodo']
        request.session['particiones'] =  request.POST['particiones']
        request.session['tipo'] = request.POST['tipo']
        if request.session['tipo'] == "doble":
            request.session['c'] =  request.POST['c']
            request.session['d'] =  request.POST['d']
        return redirect('view')