from django.shortcuts import render
from . import views

def index(request):
    context = {
        'greetings':'Welcome to forward modelling tool'
        }
   
   
    return render(request,'latarbelakang/index.html', context)
