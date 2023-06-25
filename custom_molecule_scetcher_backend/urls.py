from django.urls import path
from api import views

urlpatterns = [
    path('coordinates/', views.get_coordinates, name='get_coordinates'),
]

