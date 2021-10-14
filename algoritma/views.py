from django.shortcuts import render
from . import views

def index(request):
    context = {
        'greetings':'Welcome to forward modelling tool'
        }
   
   
    return render(request,'algoritma/index.html', context)
