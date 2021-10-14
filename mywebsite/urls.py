
from django.contrib import admin
from django.urls import path, include
from . import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('',views.index),
    path('Montecarlo/',include('forward.urls')),
    path('inverse/',include('inverse.urls')),
    path('Resis/',include('Resistivity_Section.urls')),
    path('algoritma/',include('algoritma.urls')),
    path('latarbelakang/',include('latarbelakang.urls'))
]
