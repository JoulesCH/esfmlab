from django.shortcuts import render, redirect


def creditos(request):
    creditss = [
        {'nombre': 'Rodolfo Lagunas J.', 'semestre': 4, 'carrera': 'Ingeniería Matemática - IPN'},
        {'nombre':'Julio Hernández G. ', 'semestre':4, 'carrera': 'Ingeniería Matemática - IPN, Matemáticas Aplicadas y Computación - UNAM'}
    ]
    return render(request, 'creditos.html', {'creditos':creditss})

def home(request):
    return render(request, 'home.html')

def construccion(request):
    return render(request, 'construccion.html')