from django.shortcuts import render
from . import views

def index(request):
    context = {
        'nav':[['Montecarlo/','Monte Carlo Simulation'],['inverse/','Inverse Modeling']]
        }
    return render(request,'index.html', context)