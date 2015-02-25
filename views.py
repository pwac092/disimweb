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
from itertools import combinations


def index(request):
    return render(request, 'index.html')

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
    A_mesh_unique = (list(set([i.mesh_term for i in mesh.objects.filter(omim__exact=omim_A)])))

    A_mesh = dict()
    for i in A_mesh_unique:
        try:
            A_mesh[i[0]]
        except:
            A_mesh[i[0]] = list()

        A_mesh[i[0]].append((i, "http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term="+i.replace(' ', '+')+"&field=entry"))


    B_mesh_unique  = (list(set([i.mesh_term for i in mesh.objects.filter(omim__exact=omim_B)])))

    B_mesh = dict()
    for i in B_mesh_unique:
        try:
            B_mesh[i[0]]
        except:
            B_mesh[i[0]] = list()

        B_mesh[i[0]].append((i, "http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term="+i.replace(' ', '+')+"&field=entry"))


    #A_mesh = [(i[0], i,"http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term="+i.replace(' ', '+')+"&field=entry") for i in A_mesh_unique]
    #B_mesh = [(i[0], i,"http://www.nlm.nih.gov/cgi/mesh/2014/MB_cgi?mode=&term="+i.replace(' ', '+')+"&field=entry") for i in B_mesh_unique]

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
    url_fill = reverse('fillnetwork')
    return render(request, 'explore.html', {'chosen_disease_mim':mim,'chosen_disease_name':name, 'URL':url, 'URL_fill':url_fill,'omim':disease })

def getNeighbourhood_ajax(request, disease=None, max_direct=10, max_indirect=5):
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

    #quick check to verify the limits are not exceeded.
    if int(max_direct) > 50:
        max_direct = 50
    if int(max_indirect) > 20:
        max_indirect = 20

    nodes = list()
    edges = list()

    #store the nodes that were found, to check for double loops.
    direct_k_most = [(o2, float(s),lca) for o2,s,lca in neighbourhood.objects.filter(omim1=disease).exclude(omim2=disease).order_by('-similarity').values_list('omim2', 'similarity','lca')][:int(max_direct)]

    # we get the titles, for each disease
    titles = dict(omim_details.objects.values_list('omim', 'title'))
    
    #append the pivot disease.
    nodes.append({'data': { 'id' :str(disease) , 'level': 230, 'colour': '#FFFFF', 'title': titles[str(disease)]}})

    #set weights to store the mininmum and maximum weights. Easier than traversing the entire tree.
    min_sim = 1000
    max_sim = -1

    #store the found edges
    found_edges = set()

    for omim,score, lca in direct_k_most:
        nodes.append({'data': { 'id' : omim, 'level' : 160, 'title': titles[omim]}})
        source = str(min(str(disease),str(omim)))
        target = str(max(str(disease),str(omim)))
        edges.append({'data':{ 'id': str(disease)+"_"+str(omim), 'similarity': float(score), 'source': source, 'target': target, 'LCA':lca}})
        found_edges.add((source,target))
        #store teh similarities
        min_sim = min(min_sim,float(score)) 
        max_sim = max(max_sim,float(score)) 

    ##get the 50 most similar neighbours of each of the 50 original neighbours
    for (direct_neighbour,sim, lca) in direct_k_most:
        #fetch the neighbours
        second_k_most = [(o2, float(s),lca) for o2,s,lca in neighbourhood.objects.filter(omim1=direct_neighbour).exclude(omim2=direct_neighbour).order_by('-similarity').values_list('omim2', 'similarity','lca')][:int(max_indirect)]
        #add the current node
        index = 0
        for (omim, score, lca) in second_k_most:
            source = str(min(str(direct_neighbour),str(omim)))
            target = str(max(str(direct_neighbour),str(omim)))
            #check for double loops with the first level and double loops with the pivot disease.
            if (source,target) not in found_edges:
                nodes.append({'data': {'id' : str(omim), 'level' : 50 - ((index%2) * 40), 'title': titles[str(omim)]}})
                edges.append({'data':{'id': str(source) + "_"+target, 'similarity' : float(score), 'source': source, 'target': target, 'LCA':lca}})
                found_edges.add((source,target))
                min_sim = min(min_sim,float(score)) 
                max_sim = max(max_sim,float(score)) 
            index += 1

    ##we normalise the values. This is to do with the javascript, were the linear mapping does not allow for variable limits, 
    #so we need everything between 0 and 1.
    for edge in edges:
        if max_sim == min_sim:
            edge['data']['colour'] = 1
        else:
            edge['data']['colour'] = (float(edge['data']['similarity']) - min_sim) / (max_sim - min_sim)
    ##set the json data to use
    final_set = defaultdict()
    #Just "round" the values to the closest decimal. It is just to make it look nicer.
    final_set['max_sim'] = float(max_sim)
    final_set['min_sim'] = float(min_sim)

    final_set['nodes'] = nodes
    final_set['edges'] = edges
    return HttpResponse(json.dumps(final_set), content_type = "application/json")

def fillNetwork_ajax(request):
    #get the network
    network = json.loads(request.read())
    nodes = [i['data']['id'] for i in network['nodes']]
    edges = [(min(i['data']['source'],i['data']['target']),max(i['data']['source'],i['data']['target'])) for i in network['edges']]
    #set weights to store the mininmum and maximum weights. Easier than traversing the entire tree.
    min_sim = 1000
    max_sim = -1
    #look for all the edges that are not in the list of edges.
    for putative_edge in combinations(nodes, 2):
        if (min(putative_edge),max(putative_edge)) not in edges and min(putative_edge) != max(putative_edge):
            #here we fethc this from the database.
            sim = similarityScores.objects.get(omim1__exact = min(putative_edge), omim2__exact = max(putative_edge))
            network['edges'].append({'data' : {'id' : str(min(putative_edge))+"_"+str(max(putative_edge)) , 'similarity' : float(sim.similarity), 'source' : min(putative_edge), 'target' : max(putative_edge)}})
            #max similarity should not change, so we only fetch the min simialrity
            network['min_sim'] = min(float(network['min_sim']), float(sim.similarity))
            #add the found edges to avoid double edges.
            edges.append((min(putative_edge),max(putative_edge)))

    ##we normalise the values. This is to do with the javascript, were the linear mapping does not allow for variable limits, 
    #so we need everything between 0 and 1.
    for edge in network['edges']:
        if network['max_sim'] == network['min_sim']:
            edge['data']['colour'] = 1
        else:
            edge['data']['colour'] = (float(edge['data']['similarity']) - network['min_sim']) / (network['max_sim'] - network['min_sim'])

    return HttpResponse(json.dumps(network), content_type = "application/json")


#########################
# Disease details       #
#########################

def getDetails_ajax(request, omim):
    if request.is_ajax():
        json_data = defaultdict()
        #get details of disease A
        details = get_object_or_404(omim_details, omim__exact=omim)
        #get MeSH terms
        disease_mesh = set([i.mesh_term for i in mesh.objects.filter(omim__exact=omim)])
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

