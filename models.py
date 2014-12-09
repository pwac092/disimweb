from django.db import models

# Create your models here.

#c.execute('''CREATE TABLE similarity (omim1 text, omim2 text, similarity text)''')
class similarityScores (models.Model):
    omim1 = models.CharField(max_length=10)
    omim2 = models.CharField(max_length=10)
    similarity = models.CharField(max_length=10)

#c.execute('''CREATE TABLE mimtoprot (omim text, uniprot_id text)''')
class mimtoprot(models.Model):
    omim = models.CharField(max_length=10)
    uniprot_id = models.CharField(max_length=10)

#c.execute('''CREATE TABLE ppi (interactor_a text, interactor_b text)''')
class ppi(models.Model):
    interactor_a = models.CharField(max_length=10)
    interactor_b = models.CharField(max_length=10)

#    c.execute('''CREATE TABLE mesh (omim text, mesh text)''')
class mesh(models.Model):
    omim = models.CharField(max_length=10)
    mesh_term = models.CharField(max_length=10)


#    c.execute('''CREATE TABLE omim_details (omim text, title text, prefix text)''')
class omim_details(models.Model):
    omim = models.CharField(max_length=10)
    title = models.CharField(max_length=500)
    prefix = models.CharField(max_length=2)


#CREATE TABLE "disimweb_omim_names" ("id" varchar(6) NOT NULL PRIMARY KEY, "value" varchar(100), "label" varchar(100));
class omim_names(models.Model):
    id = models.CharField(max_length=6)
    value = models.CharField(max_length=100)
    label = models.CharField(max_length=100)


