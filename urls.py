
from django.conf.urls import patterns, url
from django.views.decorators.csrf import csrf_exempt   


from disimweb import views

urlpatterns = patterns('',
            url(r'^$', views.index, name='index'),
            url(r'^explore/(\d{6})/$', views.explore, name='explore'),

            url(r'^getNeighbourhood/(\d{6})/(\d+)/(\d+)/$',views.getNeighbourhood_ajax,name='neighbourhood'),
            url(r'^getNeighbourhood/(\d{6})/$',views.getNeighbourhood_ajax,name='neighbourhood'),
            url(r'^getNeighbourhood/$',views.getNeighbourhood_ajax,name='neighbourhood'),

            url(r'^fillNetwork/(\w+)$',csrf_exempt(views.fillNetwork_ajax),name='fillnetwork'),
            url(r'^fillNetwork/$',csrf_exempt(views.fillNetwork_ajax),name='fillnetwork'),

            url(r'^details/(\d{6})/$', views.getDetails_ajax, name='details'),
            url(r'^score/(\d{6})/(\d{6})/$', views.score, name='score'),
            url(r'^error/$', views.error, name='error'),
            url(r'^disclaimer/$', views.disclaimer, name='disclaimer'),
            url(r'^get_entities/$', views.get_entities, name='get_entities'),
        )

