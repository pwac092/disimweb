from django.shortcuts import render

from disimweb.models import similarityScores, mimtoprot, ppi, mesh, omim_details
from disimweb.models import *
from django.db.models import Q

from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext, loader, Context
from django.shortcuts import get_object_or_404, render, render_to_response
from django.core.servers.basehttp import FileWrapper
from django.core.urlresolvers import reverse
from django.conf import settings

import json, math
import mimetypes
import numpy as np
from collections import defaultdict

def index(request):
    return render(request, 'index.html')

#########################
# Download              #
#########################

def respond_as_attachment(request, whatToDownload):


    response = prepareFile(whatToDownload)

    response['Content-Type'] = 'text/plain'
    response['Content-Length'] = str(os.stat(file_path).st_size)
    response['Content-Encoding'] = '7bit' #this id the default encoding.

    # To inspect details for the below code, see http://greenbytes.de/tech/tc2231/
    if u'WebKit' in request.META['HTTP_USER_AGENT']:
        # Safari 3.0 and Chrome 2.0 accepts UTF-8 encoded string directly.
        filename_header = 'filename=%s' % original_filename.encode('utf-8')
    elif u'MSIE' in request.META['HTTP_USER_AGENT']:
        # IE does not support internationalized filename at all.
        # It can only recognize internationalized URL, so we do the trick via routing rules.
        filename_header = ''
    else:
        # For others like Firefox, we follow RFC2231 (encoding extension in HTTP headers).
        filename_header = 'filename*=UTF-8\'\'%s' % urllib.quote(original_filename.encode('utf-8'))
    response['Content-Disposition'] = 'attachment; ' + filename_header
    return response

def prepareFile(what):

    if what == 1: #"simMatrix":
        path = '/static/files/combined_similarity_matrix'
    elif what == 2: #"simTriplet":
        path = '/static/files/combined_similarity_triplet'
    elif what == 3: #"omim2mesh":
        path = '/static/files/mim2mesh'
    elif what == 4: #"omim2pubmed":
        path = '/static/files/omim2pubmed'

    fp = open(path, 'rb')
    response = HttpResponse(fp.read())
    fp.close()
    return response

#########################
# Catch all             #
#########################
def disclaimer(request):
    return render(request, 'ack.html')

def error(request):
    return render(request, 'error.html',{'message':var})

#########################
# Result score          #
#########################
def ordinal(value):
    try:
        value = int(value)
    except ValueError:
        return value

    if value % 100//10 != 1:
        if value % 10 == 1:
            ordval = u"%d%s" % (value, "st")
        elif value % 10 == 2:
            ordval = u"%d%s" % (value, "nd")
        elif value % 10 == 3:
            ordval = u"%d%s" % (value, "rd")
        else:
            ordval = u"%d%s" % (value, "th")
    else:
        ordval = u"%d%s" % (value, "th")

    return ordval

def score(request, omim_A, omim_B):
    json_data = defaultdict()
    #get similarity 
    sim = similarityScores.objects.get(omim1__exact = min(omim_A,omim_B), omim2__exact = max(omim_A,omim_B))
    #get MeSH terms
    #http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term=Sensitivity+and+Specificity&field=entry
    A_mesh_unique = set([i.mesh_term for i in mesh.objects.filter(omim__exact=omim_A)])
    B_mesh_unique  = set([i.mesh_term for i in mesh.objects.filter(omim__exact=omim_B)])
    A_mesh = [(i,"http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term="+i.replace(' ', '+')+"&field=entry") for i in A_mesh_unique]
    B_mesh = [(i,"http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term="+i.replace(' ', '+')+"&field=entry") for i in B_mesh_unique]
    #get proteins
    A_proteins = [i.uniprot_id for i in mimtoprot.objects.filter(omim__exact=omim_A)]
    B_proteins = [i.uniprot_id for i in mimtoprot.objects.filter(omim__exact=omim_B)]
    #get details
    detailsA = get_object_or_404(omim_details, omim__exact=omim_A)
    #with prefix
    #name_disease_A = detailsA.prefix + " " +detailsA.title 
    name_disease_A = detailsA.title 
    detailsB = get_object_or_404(omim_details, omim__exact=omim_B)
    #with previx
    #name_disease_B = detailsB.prefix + " " +detailsB.title
    name_disease_B = detailsB.title
    return render(request, 'score.html', {'A_mesh':A_mesh, 'B_mesh':B_mesh, 'A_proteins':A_proteins,  'B_proteins':B_proteins, 'name_disease_A':name_disease_A, 'name_disease_B':name_disease_B, 'disease_A': detailsA.omim, 'disease_B':detailsB.omim , 'similarity':sim.similarity, 'percentile': ordinal(sim.percentile)})

#########################
# Explore neighbourhood #
#########################



def explore(request, disease):
    #details = omim_details.objects.get(omim__exact=disease)
    details = get_object_or_404(omim_details, omim__exact=disease)
    (mim, name, prefix) = (details.omim, details.title, details.prefix)
    url = reverse('neighbourhood')
    return render(request, 'explore.html', {'chosen_disease_mim':mim,'chosen_disease_name':name, 'URL':url, 'omim':disease })

def passThroughSigmoid(score):
     return (1.0 / (1.0 + np.exp(-2.5 * score + 5)))

def getNeighbourhood_ajax(request, disease=None):
    #json prototype
    #nodes: [
            #{ data: { id: 'a', foo: 3, bar: 5, baz: 7 } },
            #{ data: { id: 'b', foo: 7, bar: 1, baz: 3 } },
            #], 
    #
    #edges: [
            #{ data: { id: 'ae', weight: 1, source: 'a', target: 'e' } },
            #{ data: { id: 'ab', weight: 3, source: 'a', target: 'b' } },
            #]
    #};

    nodes = list()
    edges = list()
    direct_k_most = [(i.omim2, i.similarity) for i in similarityScores.objects.filter(omim1__exact=disease).order_by('-similarity')[:int(30)]]
    #append the target node.
    details = get_object_or_404(omim_details, omim__exact=disease)
    nodes.append({'data': { 'id' :str(disease) , 'level': 290, 'colour': '#FFFFF', 'title': details.title}})
    for (j,i) in enumerate(direct_k_most):
        if str(disease) != str(i[0]):
            details = get_object_or_404(omim_details, omim__exact=str(i[0]))
            nodes.append({'data': { 'id' : i[0], 'level' : 160, 'colour': '#B3767E', 'title': details.title}})
            edges.append({ 'data':{ 'id': str(disease)+"_"+str(i[0]), 'similarity': float(i[1]), 'colour' : passThroughSigmoid((float(i[1]))), 'source': str(disease), 'target': str(i[0])}})

    ##get the 10 most similar neighbours of each of the 50 original neighbours
    for (neighbour_disease,sim) in direct_k_most:
        #fetch the neighbours
        second_k_most = [(i.omim2, i.similarity) for i in similarityScores.objects.filter(omim1__exact=neighbour_disease).order_by('-similarity')[:int(5)]]
        #add the current node
        #nodes.append({'data': { 'id' : str(neighbour_disease) , 'level':'2'}})
        for (index,pair) in enumerate(second_k_most):
            if str(neighbour_disease) != str(pair[0]):
                details = get_object_or_404(omim_details, omim__exact=str(pair[0]))
                nodes.append({'data': {'id' : str(pair[0]), 'level' : 50 - ((index%2) * 40), 'colour': '#B3767E', 'title': details.title}})
                edges.append({ 'data':{'id': str(neighbour_disease)+"_"+str(pair[0]), 'similarity' : float(pair[1]), 'colour': passThroughSigmoid(float(pair[1])), 'source': str(neighbour_disease), 'target': str(pair[0])}})
    ##set the json data to use
    final_set = defaultdict()
    final_set['nodes'] = nodes
    final_set['edges'] = edges
    return HttpResponse(json.dumps(final_set), content_type = "application/json")


#########################
# Disease details       #
#########################

def getDetails_ajax(request, omim):
    if request.is_ajax():
        json_data = defaultdict()
        #get details of disease A
        details = get_object_or_404(omim_details, omim__exact=omim)
        #get MeSH terms
        disease_mesh = [i.mesh_term for i in mesh.objects.filter(omim__exact=omim)]
        #get proteins
        disease_proteins = [i.uniprot_id for i in mimtoprot.objects.filter(omim__exact=omim)]
        #create json data
        json_data = {'name':{'mim_no':details.omim, 'title':details.title, 'prefix':details.prefix}, 'proteins': disease_proteins, 'mesh': disease_mesh}
        return HttpResponse(json.dumps(json_data), content_type = "application/json")
    else:
        return render(request, 'error.html',{'message':''})



#########################
# Autocomplete          #
#########################

def get_entities(request):
    if request.is_ajax():
        data = ''
        # get data from the request.
        query = request.GET.get('term', '')
        # filter names
        results = omim_names.objects.filter(Q(label__icontains=query) | Q(value__icontains=query)).order_by('value')[:20]
        results = results.values_list('value', 'label')
        data = json.dumps([{"label": str(i[1]) + "("+str(i[0]) + ")", "value": str(i[0])} for i in results])
    else:
        data = 'fail'
    mimetype = 'application/json'
    return HttpResponse(data, mimetype)

