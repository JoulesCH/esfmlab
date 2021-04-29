
from django.contrib import admin
from django.urls import path
from symboesfm import views as home_view

urlpatterns = [
    path('admin/', admin.site.urls),
    path('view/', home_view.view, name = 'view'),
    path('simple/', home_view.simple, name =  "simple"),
    path('doble/', home_view.doble, name =  "doble"),
    path('submit/', home_view.submit, name = 'submit'),
    path('creditos/', home_view.creditos, name = 'creditos'),
    path('', home_view.home, name = 'home')
]+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
