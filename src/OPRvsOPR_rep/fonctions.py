#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author : Shogofa MORTAZA
IBPC Internship
"""

#import modules
import os
import sys
from math import *
import glob

#Fonctions

def convert_blast_to_cytoscape(blast_file, cytoscape, seuil):
	"""
	Fonction permettant de preparer un fichier pour Cytoscape (ou autre logiciel de visualisation).
	input : 
		- blast-file : nom fichier texte contenant les colonnes P1 P2 Evalue à partir du résultat Blast
		- cytoscape : nom fichier texte de sortie pour Cytoscape
	output : 
		- fichier texte pret a etre traite par Cytoscape (sans repetitions ...)
	"""

	#enregistrement de toutes les lignes du fichier blast
	filin = open(blast_file,'r')
	lines = filin.readlines()
	filin.close()

	filout = open(cytoscape,'w')
	filout.write("OPR1\tOPR2\tEvalue\n")

	PPI_dico = {}

	#pour chaque ligne, mettre les infos dans un dico avec cle = PPI et valeur = Evalue (!= 0)
	for l in lines :
		#print l
		P1, P2, Evalue = l.split()
		if P1 != P2 :
			#print P1[1:-1], P2[1:-1], Evalue, log_Evalue, -float(log_Evalue)
			if P1[1:-1]+"vs"+P2[1:-1] not in PPI_dico :
				PPI_dico[P1+"vs"+P2] = []
			PPI_dico[P1+"vs"+P2].append(Evalue)
		#if P1+"vs"+P2 not in PPI_dico and Evalue != "0.0":
		#	PPI_dico[P1+"vs"+P2] = []
		#if Evalue != "0.0" :
		#	PPI_dico[P1+"vs"+P2].append(-log(float(Evalue)))

	#pour chaque PPI, ne selectionner que celles la Evalue min
	#au fur et a mesure, enregistrement des résultats dans un fichier
	for key in PPI_dico.keys():
		P1,P2 = key.split("vs")
		#print PPI_dico[key]
		#print P1 + "\t" + P2 + "\t" + str(max(PPI_dico[key]))
		if float(min(PPI_dico[key])) <= float(seuil) :
			#print float(min(PPI_dico[key]))
			#print str(P1) + "\t" + str(P2) + "\t" + str(min(PPI_dico[key])) + "\t"
			#filout.write(str(P1) + "\t" + str(P2) + "\t" + str(min(PPI_dico[key])) + "\t" + "\n")
			filout.write(P1 + "\t" + P2 + "\t" + str(max(PPI_dico[key])) + "\n")

	filout.close()


def convert_cyto_to_2csv_Gephi(cyto_file, nodes_csv, edges_csv):
	"""
	Fonction permettant de créer deux fichiers csv pour Gephi.
	input :
		- cyto_file : nom fichier texte contenant les colonnes P1 P2 Evalue à partir du résultat Blast
		- nodes_csv : nom fichier csv de sortie pour les nodes de Gephi
		- edges_csv : nom fichier csv de sortie pour les edges de Gephi
	output :
		- deux fichiers format csv prets a etre traite par Gephi
	"""
	
	#enregistrement de toutes les lignes du fichier blast
	filin = open(cyto_file, 'r')
	lines = filin.readlines()
	filin.close()

	filout1 = open(edges_csv, 'w')
	filout1.write("Source\tTarget\tWeight\n")
	
	#compteur pour donner indices aux OPRs
	OPR_id = 1
	#dico avec cle = nom de l'OPR et valeur = indice donnée à l'OPR
	dico_OPR_id = {}
	
	for l in lines[1:] : 
		P1, P2, Evalue = l.split()
		if P1 not in dico_OPR_id :
			dico_OPR_id[P1] = OPR_id
			OPR_id = OPR_id + 1
		if P2 not in dico_OPR_id :
			dico_OPR_id[P2] = OPR_id
			OPR_id = OPR_id + 1
		#print dico_OPR_id[P1], dico_OPR_id[P2], Evalue
		filout1.write(str(dico_OPR_id[P1]) + "\t" + str(dico_OPR_id[P2]) + "\t" + str(Evalue) + "\n")

	filout2 = open(nodes_csv, 'w')
	filout2.write("Id\tLabel\n")


	for OPR in dico_OPR_id :
		#print dico_OPR_id[OPR], OPR
		filout2.write(str(dico_OPR_id[OPR]) + "\t" + str(OPR) + "\n")

	filout1.close()
	filout2.close()

def cluster_count(cluster_file) :
	"""
	Fonction qui permet de donner la taille du cluster, le nb proteines ncl, regular et collante;
	input :
		- cluster_file : Fichier exporté à partir de cytoscape et qui contient tous les elements d'un seul cluster
	output :
		- print les infos qu'on veut
	"""
	#Liste des gènes OPR par catégorie
	ncl_list = ["Chlre_OPR106", "Chlre_OPR70", "Chlre_OPR77", "Chlre_OPR84", "Chlre_OPR91", "Chlre_OPR98", "Chlre_OPR110",\
	 "Chlre_OPR71", "Chlre_OPR78", "Chlre_OPR85", "Chlre_OPR92", "Chlre_OPR99", "Chlre_OPR111", "Chlre_OPR72", "Chlre_OPR79", \
	 "Chlre_OPR86", "Chlre_OPR93", "Chlre_OPR112", "Chlre_OPR73", "Chlre_OPR80", "Chlre_OPR87", "Chlre_OPR94", "Chlre_OPR20", \
	 "Chlre_OPR74", "Chlre_OPR81", "Chlre_OPR88", "Chlre_OPR95", "Chlre_OPR21", "Chlre_OPR75", "Chlre_OPR82", "Chlre_OPR89", \
	 "Chlre_OPR96", "Chlre_OPR69", "Chlre_OPR76", "Chlre_OPR83", "Chlre_OPR90", "Chlre_OPR97"]
	collantes_list = ["Chlre_OPR105", "Chlre_OPR107", "Chlre_OPR28", "Chlre_OPR29", "Chlre_OPR32", "Chlre_OPR62"]
	regular_list = ["Chlre_OPR100", "Chlre_OPR118", "Chlre_OPR1", "Chlre_OPR36", "Chlre_OPR49", "Chlre_OPR60", \
	"Chlre_OPR101", "Chlre_OPR119", "Chlre_OPR22", "Chlre_OPR37", "Chlre_OPR4", "Chlre_OPR61",\
	"Chlre_OPR102", "Chlre_OPR11", "Chlre_OPR23", "Chlre_OPR38", "Chlre_OPR50", "Chlre_OPR63", \
	"Chlre_OPR103", "Chlre_OPR120", "Chlre_OPR24", "Chlre_OPR39", "Chlre_OPR51", "Chlre_OPR64",\
	"Chlre_OPR104", "Chlre_OPR121", "Chlre_OPR25", "Chlre_OPR40", "Chlre_OPR52", "Chlre_OPR65", \
	"Chlre_OPR108", "Chlre_OPR12", "Chlre_OPR26", "Chlre_OPR41", "Chlre_OPR53", "Chlre_OPR66", \
	"Chlre_OPR109", "Chlre_OPR13", "Chlre_OPR27", "Chlre_OPR42", "Chlre_OPR54", "Chlre_OPR67", \
	"Chlre_OPR10", "Chlre_OPR14", "Chlre_OPR2", "Chlre_OPR43", "Chlre_OPR55", "Chlre_OPR68", \
	"Chlre_OPR113", "Chlre_OPR15", "Chlre_OPR30", "Chlre_OPR44", "Chlre_OPR56", "Chlre_OPR6", \
	"Chlre_OPR114", "Chlre_OPR16", "Chlre_OPR31", "Chlre_OPR45A", "Chlre_OPR57", "Chlre_OPR7",\
	"Chlre_OPR115", "Chlre_OPR17", "Chlre_OPR33", "Chlre_OPR46", "Chlre_OPR58", "Chlre_OPR8", \
	"Chlre_OPR116", "Chlre_OPR18", "Chlre_OPR34", "Chlre_OPR47", "Chlre_OPR59", "Chlre_OPR9", \
	"Chlre_OPR117", "Chlre_OPR19", "Chlre_OPR35", "Chlre_OPR48", "Chlre_OPR5", "Chlre_OPR67=Cre14.g633150.t1.3"]
	# dico contenant les clusters (chaque cluster ayant 3 valeurs (nn, nc et nr) comptant le nombre de genes ds chaque gpe)
	dico_cluster = {}
	#ouverture du fichier contenant les clusters 
	filin = open(cluster_file, 'r')
	l = filin.readline()
	for l in filin :
		num_cluster = l.split(',')[1][1:-1]
		gene = l.split(',')[2][1:-1]
		#si le cluster n'est pas ds le dico, l'ajouter
		if num_cluster not in dico_cluster :
			dico_cluster[num_cluster] = [0, 0, 0]
		#si il est dans le dico, commencer le comptage de chaque gpe de gene
		if num_cluster in dico_cluster :
			if gene in ncl_list :
				dico_cluster[num_cluster][0] = dico_cluster[num_cluster][0] + 1
			elif gene in collantes_list :
				dico_cluster[num_cluster][1] = dico_cluster[num_cluster][1] + 1
			elif gene in regular_list :
				dico_cluster[num_cluster][2] = dico_cluster[num_cluster][2] + 1
			else : 
				print gene , " = gene non classé"
	filin.close()
	#print len(dico_cluster)
	#print dico_cluster
	print cluster_file
	print "cluster [ncl, collantes, regular"
	for i in dico_cluster :
		print i, dico_cluster[i]

def mcl_command(filin, Istart = 2, Iend = 7, Istep = 1) :
	"""
	Fonction permettant de lancer des lignes de commandes mcl 
	input :
		- filin : fichier donné en entrée avec 3 colonnes séparées par des \t (OPR1, OPR2 et Evalue)
		- Istart : valeur initiale de I (par defaut, commence à 2.0)
		- Iend : valeur finale de I (par defaut, termine à 6.0)
		- Istep : pas de I
	output :
		- creation de fichiers par defaut contenant les resultats de mcl
	"""
	for i in range(int(Istart), int(Iend)) :
		os.system('/home/mortaza/local/bin/mcl ' + filin +' --abc --abc-neg-log10 -I '+ str(i))

def mcl_analysis(filin, rep_name) :
	"""
	Fonction permettant d'analyser les resultats de mcl_command
	input : 
		- filin : fichier format abc donné en entrée de mcl_command
		- rep_name : nom du rep où sera stocké tous les fichiers resultats de mcl
	output :
		- 
	"""
	os.system("mkdir " + rep_name)
	os.system("mv out."+filin+".I* " + rep_name)
	files = glob.glob(rep_name+'/*')
	filout = open(rep_name+"/"+rep_name+"_result", "w")
	for i in files : 
		#print i
		filout.write(i+"\n")
		#print "cluster [ncl, collantes, regular]"
		filout.write("cluster [ncl, collantes, regular]\n")
		dico_cluster = mcl_cluster_count(i)
		for cluster in dico_cluster :
			#print cluster, dico_cluster[cluster]
			filout.write(str(cluster)+"\t"+str(dico_cluster[cluster])+"\n")
	filout.close()

def mcl_cluster_count(mcl_outfile) :
	"""
	Fonction qui permet de donner la taille du cluster, le nb proteines ncl, regular et collante d'un fichier mcl
	input : 
		- mcl_outfile : fichier output de mcl_command
	output : 
		- retourne la taille du cluster et sa composition
	"""
	#Liste des gènes OPR par catégorie
	ncl_list = ["Chlre_OPR106", "Chlre_OPR70", "Chlre_OPR77", "Chlre_OPR84", "Chlre_OPR91", "Chlre_OPR98", "Chlre_OPR110",\
	 "Chlre_OPR71", "Chlre_OPR78", "Chlre_OPR85", "Chlre_OPR92", "Chlre_OPR99", "Chlre_OPR111", "Chlre_OPR72", "Chlre_OPR79", \
	 "Chlre_OPR86", "Chlre_OPR93", "Chlre_OPR112", "Chlre_OPR73", "Chlre_OPR80", "Chlre_OPR87", "Chlre_OPR94", "Chlre_OPR20", \
	 "Chlre_OPR74", "Chlre_OPR81", "Chlre_OPR88", "Chlre_OPR95", "Chlre_OPR21", "Chlre_OPR75", "Chlre_OPR82", "Chlre_OPR89", \
	 "Chlre_OPR96", "Chlre_OPR69", "Chlre_OPR76", "Chlre_OPR83", "Chlre_OPR90", "Chlre_OPR97"]
	collantes_list = ["Chlre_OPR105", "Chlre_OPR107", "Chlre_OPR28", "Chlre_OPR29", "Chlre_OPR32", "Chlre_OPR62"]
	regular_list = ["Chlre_OPR100", "Chlre_OPR118", "Chlre_OPR1", "Chlre_OPR36", "Chlre_OPR49", "Chlre_OPR60", \
	"Chlre_OPR101", "Chlre_OPR119", "Chlre_OPR22", "Chlre_OPR37", "Chlre_OPR4", "Chlre_OPR61",\
	"Chlre_OPR102", "Chlre_OPR11", "Chlre_OPR23", "Chlre_OPR38", "Chlre_OPR50", "Chlre_OPR63", \
	"Chlre_OPR103", "Chlre_OPR120", "Chlre_OPR24", "Chlre_OPR39", "Chlre_OPR51", "Chlre_OPR64",\
	"Chlre_OPR104", "Chlre_OPR121", "Chlre_OPR25", "Chlre_OPR40", "Chlre_OPR52", "Chlre_OPR65", \
	"Chlre_OPR108", "Chlre_OPR12", "Chlre_OPR26", "Chlre_OPR41", "Chlre_OPR53", "Chlre_OPR66", \
	"Chlre_OPR109", "Chlre_OPR13", "Chlre_OPR27", "Chlre_OPR42", "Chlre_OPR54", "Chlre_OPR67", \
	"Chlre_OPR10", "Chlre_OPR14", "Chlre_OPR2", "Chlre_OPR43", "Chlre_OPR55", "Chlre_OPR68", \
	"Chlre_OPR113", "Chlre_OPR15", "Chlre_OPR30", "Chlre_OPR44", "Chlre_OPR56", "Chlre_OPR6", \
	"Chlre_OPR114", "Chlre_OPR16", "Chlre_OPR31", "Chlre_OPR45A", "Chlre_OPR57", "Chlre_OPR7",\
	"Chlre_OPR115", "Chlre_OPR17", "Chlre_OPR33", "Chlre_OPR46", "Chlre_OPR58", "Chlre_OPR8", \
	"Chlre_OPR116", "Chlre_OPR18", "Chlre_OPR34", "Chlre_OPR47", "Chlre_OPR59", "Chlre_OPR9", \
	"Chlre_OPR117", "Chlre_OPR19", "Chlre_OPR35", "Chlre_OPR48", "Chlre_OPR5", "Chlre_OPR67=Cre14.g633150.t1.3"]
	# dico contenant les clusters (chaque cluster ayant 3 valeurs (nn, nc et nr) comptant le nombre de genes ds chaque gpe)
	dico_cluster = {}
	#numero de cluster
	num_cluster = 0
	#os.system("wc -l "+mcl_outfile)
	#ouverture du fichier contenant les clusters 
	f = open(mcl_outfile, 'r')
	for line in f :
		num_cluster = num_cluster + 1 
		#print num_cluster
		dico_cluster[num_cluster] = [0, 0, 0]
		nb = line.count('\t')
		#print line[:-1]
		for sep in range(nb + 1) :
			gene = line.split()[sep]
			#print gene
			if gene in ncl_list :
				dico_cluster[num_cluster][0] = dico_cluster[num_cluster][0] + 1
				#print gene
			elif gene in collantes_list :
				dico_cluster[num_cluster][1] = dico_cluster[num_cluster][1] + 1
			elif gene in regular_list :
				dico_cluster[num_cluster][2] = dico_cluster[num_cluster][2] + 1
			else : 
				print gene , " = gene non classé"
	f.close()
	return dico_cluster

def opr_type_in_file():
	"""
	Fonction permettant de creer un fichier de prot et leur type
	"""
	#dico des gènes OPR par catégorie
	dico_type = {'ncl_list' : ["Chlre_OPR106", "Chlre_OPR70", "Chlre_OPR77", "Chlre_OPR84", "Chlre_OPR91", "Chlre_OPR98", "Chlre_OPR110",\
	 "Chlre_OPR71", "Chlre_OPR78", "Chlre_OPR85", "Chlre_OPR92", "Chlre_OPR99", "Chlre_OPR111", "Chlre_OPR72", "Chlre_OPR79", \
	 "Chlre_OPR86", "Chlre_OPR93", "Chlre_OPR112", "Chlre_OPR73", "Chlre_OPR80", "Chlre_OPR87", "Chlre_OPR94", "Chlre_OPR20", \
	 "Chlre_OPR74", "Chlre_OPR81", "Chlre_OPR88", "Chlre_OPR95", "Chlre_OPR21", "Chlre_OPR75", "Chlre_OPR82", "Chlre_OPR89", \
	 "Chlre_OPR96", "Chlre_OPR69", "Chlre_OPR76", "Chlre_OPR83", "Chlre_OPR90", "Chlre_OPR97"],
	'collantes_list' : ["Chlre_OPR105", "Chlre_OPR107", "Chlre_OPR28", "Chlre_OPR29", "Chlre_OPR32", "Chlre_OPR62"],
	'regular_list' : ["Chlre_OPR100", "Chlre_OPR118", "Chlre_OPR1", "Chlre_OPR36", "Chlre_OPR49", "Chlre_OPR60", \
	"Chlre_OPR101", "Chlre_OPR119", "Chlre_OPR22", "Chlre_OPR37", "Chlre_OPR4", "Chlre_OPR61",\
	"Chlre_OPR102", "Chlre_OPR11", "Chlre_OPR23", "Chlre_OPR38", "Chlre_OPR50", "Chlre_OPR63", \
	"Chlre_OPR103", "Chlre_OPR120", "Chlre_OPR24", "Chlre_OPR39", "Chlre_OPR51", "Chlre_OPR64",\
	"Chlre_OPR104", "Chlre_OPR121", "Chlre_OPR25", "Chlre_OPR40", "Chlre_OPR52", "Chlre_OPR65", \
	"Chlre_OPR108", "Chlre_OPR12", "Chlre_OPR26", "Chlre_OPR41", "Chlre_OPR53", "Chlre_OPR66", \
	"Chlre_OPR109", "Chlre_OPR13", "Chlre_OPR27", "Chlre_OPR42", "Chlre_OPR54", "Chlre_OPR67", \
	"Chlre_OPR10", "Chlre_OPR14", "Chlre_OPR2", "Chlre_OPR43", "Chlre_OPR55", "Chlre_OPR68", \
	"Chlre_OPR113", "Chlre_OPR15", "Chlre_OPR30", "Chlre_OPR44", "Chlre_OPR56", "Chlre_OPR6", \
	"Chlre_OPR114", "Chlre_OPR16", "Chlre_OPR31", "Chlre_OPR45A", "Chlre_OPR57", "Chlre_OPR7",\
	"Chlre_OPR115", "Chlre_OPR17", "Chlre_OPR33", "Chlre_OPR46", "Chlre_OPR58", "Chlre_OPR8", \
	"Chlre_OPR116", "Chlre_OPR18", "Chlre_OPR34", "Chlre_OPR47", "Chlre_OPR59", "Chlre_OPR9", \
	"Chlre_OPR117", "Chlre_OPR19", "Chlre_OPR35", "Chlre_OPR48", "Chlre_OPR5", "Chlre_OPR67=Cre14.g633150.t1.3"]}
	filout = open("OPR_type.txt", "w")
	filout.write("name\ttype\n")
	#pour chaque type
	for t in dico_type :
		#pour chaque opr 
		for opr in dico_type[t] :
			filout.write(opr + "\t" + t.split('_')[0] + "\n")
	filout.close()

def OPR_dom_list(file_dom, output = "domain_list.txt") :
	"""
	Fonction permettant de donner la liste des domaines possibles dans un fichier
	input :
		- file_dom : format abc
		- output : fichier de sortie contenant la liste des domaines
	output :
		- retourne la liste des dom existants dans un fichier
	"""
	#liste des OPR
	OPR_list = []
	#lecture des fichiers
	f = open(file_dom, 'r')
	l = f.readline()
	for l in f :
		if l.split()[0] not in OPR_list :
			OPR_list.append(l.split()[0])
		if l.split()[1] not in OPR_list :
			OPR_list.append(l.split()[1])
	f.close()
	#print len(OPR_list)
	#liste des dom
	dom_list = []
	for OPR in OPR_list :
		n = OPR.count("_")
		#print OPR
		dom = ""
		for sep in range(2,n+1) :
			dom = dom + '_' + OPR.split('_')[sep]
		dom = dom[1:]
		if dom not in dom_list :
			#print dom
			dom_list.append(dom)
	filout = open(output, "w")
	for i in dom_list :
		filout.write(i+"\n")
	filout.close()

def what_type_OPR(OPR_name) :
	"""
	Fonction donnant le type / la catégorie de l'OPR en question
	input : 
		- OPR_name : nom de l'OPR
	output :
		- retourne le type (n - ncl-, c -collante- ou r-regular-)
	"""
	#dico des gènes OPR par catégorie
	dico_type = {'ncl_list' : ["Chlre_OPR106", "Chlre_OPR70", "Chlre_OPR77", "Chlre_OPR84", "Chlre_OPR91", "Chlre_OPR98", "Chlre_OPR110",\
	 "Chlre_OPR71", "Chlre_OPR78", "Chlre_OPR85", "Chlre_OPR92", "Chlre_OPR99", "Chlre_OPR111", "Chlre_OPR72", "Chlre_OPR79", \
	 "Chlre_OPR86", "Chlre_OPR93", "Chlre_OPR112", "Chlre_OPR73", "Chlre_OPR80", "Chlre_OPR87", "Chlre_OPR94", "Chlre_OPR20", \
	 "Chlre_OPR74", "Chlre_OPR81", "Chlre_OPR88", "Chlre_OPR95", "Chlre_OPR21", "Chlre_OPR75", "Chlre_OPR82", "Chlre_OPR89", \
	 "Chlre_OPR96", "Chlre_OPR69", "Chlre_OPR76", "Chlre_OPR83", "Chlre_OPR90", "Chlre_OPR97"],
	'collantes_list' : ["Chlre_OPR105", "Chlre_OPR107", "Chlre_OPR28", "Chlre_OPR29", "Chlre_OPR32", "Chlre_OPR62"],
	'regular_list' : ["Chlre_OPR100", "Chlre_OPR118", "Chlre_OPR1", "Chlre_OPR36", "Chlre_OPR49", "Chlre_OPR60", \
	"Chlre_OPR101", "Chlre_OPR119", "Chlre_OPR22", "Chlre_OPR37", "Chlre_OPR4", "Chlre_OPR61",\
	"Chlre_OPR102", "Chlre_OPR11", "Chlre_OPR23", "Chlre_OPR38", "Chlre_OPR50", "Chlre_OPR63", \
	"Chlre_OPR103", "Chlre_OPR120", "Chlre_OPR24", "Chlre_OPR39", "Chlre_OPR51", "Chlre_OPR64",\
	"Chlre_OPR104", "Chlre_OPR121", "Chlre_OPR25", "Chlre_OPR40", "Chlre_OPR52", "Chlre_OPR65", \
	"Chlre_OPR108", "Chlre_OPR12", "Chlre_OPR26", "Chlre_OPR41", "Chlre_OPR53", "Chlre_OPR66", \
	"Chlre_OPR109", "Chlre_OPR13", "Chlre_OPR27", "Chlre_OPR42", "Chlre_OPR54", "Chlre_OPR67", \
	"Chlre_OPR10", "Chlre_OPR14", "Chlre_OPR2", "Chlre_OPR43", "Chlre_OPR55", "Chlre_OPR68", \
	"Chlre_OPR113", "Chlre_OPR15", "Chlre_OPR30", "Chlre_OPR44", "Chlre_OPR56", "Chlre_OPR6", \
	"Chlre_OPR114", "Chlre_OPR16", "Chlre_OPR31", "Chlre_OPR45A", "Chlre_OPR57", "Chlre_OPR7",\
	"Chlre_OPR115", "Chlre_OPR17", "Chlre_OPR33", "Chlre_OPR46", "Chlre_OPR58", "Chlre_OPR8", \
	"Chlre_OPR116", "Chlre_OPR18", "Chlre_OPR34", "Chlre_OPR47", "Chlre_OPR59", "Chlre_OPR9", \
	"Chlre_OPR117", "Chlre_OPR19", "Chlre_OPR35", "Chlre_OPR48", "Chlre_OPR5", "Chlre_OPR67=Cre14.g633150.t1.3"]}
	for i in dico_type :
		#print OPR_name
		#print dico_type[i]
		if OPR_name in dico_type[i] : 
			#print "ok"
			if i == 'ncl_list' :
				return 'n'
			elif i == 'collantes_list' :
				return 'c'
			else : 
				return 'r'
	return 'no idea'



def analyse_MCL_dom_file(output_MCL) :
	"""
	Fonction permettant d'analyser les resultats de MCL obtenus avec les domaines OPR pr un seul fichier
	input :
		- output_MCL : fichier de sortie de commande MCL
	output :
		- un fichier texte regroupant le resultat de l'analyse du output_MCL
	"""
	#fichier avec nombre de n ou r ou c dans dom
	filout2 = open(output_MCL+"_type_prot_by_dom.txt", "w")
	#header
	filout2.write("***\nNumero_cluster --> Nb_dom_diff / Nombre_dom_total\n***\n")
	filout2.write("Domain : NCL [nb], Regular [nb], Collantes [nb]\n\n")
	#fichier avec descriptif complet
	#filout1 = open(output_MCL+"_complete_descriptif.txt", "w")
	#header
	#filout1.write("***\nNumero_cluster --> Nb_dom_diff / Nombre_dom_total\n***\n")
	#filout1.write("NCL [%], Regular [%], Collantes [%]\n\n")
	#fichier de sortie
	#filout = open(output_MCL+"_descriptif.txt", "w")
	#numero de cluster
	num_cluster = 0
	#lecture du fichier
	f = open(output_MCL, 'r')
	for l in f :
		if l != "OPR1\n" and l != "OPR2\n" :
			#pour un cluster
			num_cluster += 1
			#print num_cluster
			nb_dom = l.count('\t')
			#print nb_dom
			#print "***\n" + str(num_cluster),
			#dico de occ de chaque dom
			dico_dom = {}
			#pour chaque proteine du cluster
			for i in range(nb_dom + 1) :
				prot_dom = l.split()[i]
				n = prot_dom.count("_")
				#pour chaque opr
				prot = ""
				#pour chaque dom
				dom = ""
				for sep in range(0,n+1) :
					if sep == 0 or sep == 1 :
						prot = prot + '_' + prot_dom.split('_')[sep]
					else :
						dom = dom + '_' + prot_dom.split('_')[sep]

				prot = prot[1:]
				#print prot
				dom = dom[1:]
				#print dom
				#si le dom n'est pas ds le dico
				if dom not in dico_dom :
					#initialiser liste des prot
					list_prot = [prot]
					#initialiser occ du dom 
					dico_dom[dom] = [1, list_prot]
				#sinon si le dom est ds le dico
				else :
					#incrementer de 1 la valeur de son occurence
					dico_dom[dom][0] += 1
					#concernant la prot, si existe deja l'incrementer
					if prot not in dico_dom[dom][1] :
						dico_dom[dom][1].append(prot)
					#sinon l'initialiser
					else : 
						list = 1
			#print "(" + str(len(dico_dom)) + ") / " + str(sum(dico_dom.values())) + "\n***"
			#filout.write("***\n" + str(num_cluster) + "(" + str(len(dico_dom)) + ") / " + str(sum(dico_dom.values())) + "\n***\n")
			#le nombre de dom au total
			tot_dom = 0
			for s in dico_dom.values() :
				tot_dom = tot_dom + s[0]
			#print tot_dom
			#print len(dico_dom) #le nb de dom diff
			#print "***\n" + str(num_cluster) + " --> " + str(len(dico_dom)) + " / " + str(tot_dom) + "\n***\n"
			filout2.write("***\n" + str(num_cluster) + " --> " + str(len(dico_dom)) + " / " + str(tot_dom) + "\n***\n")
			#filout1.write("***\n" + str(num_cluster) + " --> " + str(len(dico_dom)) + " / " + str(tot_dom) + "\n***\n")
			for i in dico_dom :
				#nb de chaque type occ
				nn = 0
				nc = 0
				nr = 0
				#print str(i) + " : ",
				filout2.write(str(i) + " : ")
				#print str(i) + " : " + str(dico_dom[i][0]) + "\n"
				for j in dico_dom[i][1] :
					#print "\t-" + str(j) + " " + what_type_OPR(j) + "\n" 
					if what_type_OPR(j) == "n" :
						nn += 1
					elif what_type_OPR(j) == "r" :
						nr += 1
					else :
						nc += 1
				#print str(i) + " : " + str(dico_dom[i])
				#filout.write(str(i) + " : " + str(dico_dom[i]) + "\n")	
				#print "N [" + str(nn) + "], R [" + str(nr) + "], C [" + str(nc) + "]\n"
				filout2.write("N [" + str(nn) + "], R [" + str(nr) + "], C [" + str(nc) + "]\n")
			#filout1.write("N [" + str(100*nn/tot_dom) + "], R [" + str(100*nr/tot_dom) + "], C [" + str(100*nc/tot_dom) + "]\n")
	f.close()
	#filout.close()
	#filout1.close()
	filout2.close()


def analyse_MCL_dom(prefix_file) :
	"""
	Fonction permettant d'analyser l'ensemble des fichiers output MCL pour les domaines
	input :
		- prefix_file : nom id pour l'ens des fichiers output MCL domaines
		/!\ Se placer dans le repertoire contenant l'ens de ces fichiers
	output : 
		- pour chaque output du MCL, un fichier descriptif de l'analyse est créé
	"""
	files = glob.glob("./" + prefix_file + "*")
	for file in files :
		#analyse_MCL_dom_file(file)
		#classif1_in_classif2(file)
		occ_dom_by_domcluster(file)


def classif1_in_classif2(classif2, classif1 = "/data/mortaza/Lucifer/oprvsopr/genes/out1/out.opr_vs_opr_0.001.txt.I50") :
	"""
	Fonction qui permet d'analyser les domaines à partir du clustering des gènes
	input :
		- classif1 : fichier mcl classif des genes OPR (I = 5.0 par defaut)
	output :
	"""
	#Infos classif1
	#dico des infos de classif1
	dico_classif1 = {}
	#initialisation du num de cluster
	num_cluster = 0
	filin1 = open(classif1, "r")
	for l in filin1 :
		#changement cluster
		num_cluster += 1
		if l != "OPR1\n" and l != "OPR2\n"  :
			#print num_cluster
			dico_classif1[num_cluster] = []
			#print l
			nb_prot_1 = l.count('\t')
			#print nb_prot_1
			#separation des noms des prot
			for sep in range(nb_prot_1 + 1) :
				prot = l.split()[sep]
				dico_classif1[num_cluster].append(prot)
	filin1.close()
	#verif1
	#for i in dico_classif1 :
		#print i, dico_classif1[i]
	#Infos classif2
	#dico des infos de classif2
	dico_classif2 = {}
	#initialisation du num de cluster
	num_cluster = 0
	filin2 = open(classif2, "r")
	for l in filin2 :
		#changement cluster
		num_cluster += 1
		if l != "OPR1\n" and l != "OPR2\n"  :
			#print num_cluster
			dico_classif2[num_cluster] = []
			#print l
			nb_prot_2 = l.count('\t')
			#print nb_prot_2
			#separation des noms des prot
			for sep in range(nb_prot_2 + 1) :
				prot = l.split()[sep]
				prot = prot.split('_')[0] + '_' + prot.split('_')[1]
				if prot not in dico_classif2[num_cluster] :
					#print prot
					dico_classif2[num_cluster].append(prot)
	filin2.close()
	#verif2
	#for i in dico_classif2 :
		#print i, dico_classif2[i]
	#Classif1 vs Classif2
	#dico contenant le resultat sous forme : dicoc1 and dicoc2
	#{N°clusterC2 : {N°clusterC1 : [liste_prot]}}
	dicoc2 = {}
	for c2_cluster in dico_classif2 :
		#print c2_cluster
		dicoc1 = {}
		dicoc2[c2_cluster] = dicoc1 
		for prot in dico_classif2[c2_cluster] :
			#print prot
			for c1_cluster in dico_classif1 :
				if prot in dico_classif1[c1_cluster] :
					#print prot, c1_cluster
					if c1_cluster not in dicoc1 :
						dicoc1[c1_cluster] = []
						dicoc1[c1_cluster].append(prot)
					else : 
						#eviter de compter une meme prot plusieurs fois
						if prot not in dicoc1[c1_cluster] :
							dicoc1[c1_cluster].append(prot)
					break
		#print dicoc1
	#verif et ecriture dans un fichier
	output_file = classif2 + "_gene" + classif1[-4:]
	filout = open(output_file, 'w')
	filout.write("N°_Cluster_Dom ==> N°_Cluster_Gene [%_prot]\n")
	for i in dicoc2 :
		if dicoc2[i] != {} :
			#print i
			#print str(i) + " ==> ",
			filout.write(str(i) + " ==> ")
			#calculer le total des prot dans le cluster
			total = 0.0
			for j in dicoc2[i] :
				total = total + len(dicoc2[i][j]) 
			for j in dicoc2[i] :
				#print str(j) + " [" + str(len(dicoc2[i][j])) + "] ",
				filout.write(str(j) + " [" + str(100 * len(dicoc2[i][j]) / total) + "] ")
			#print "\n"
			filout.write("\n")
	filout.close()

def c1c2_analyse(c1c2_file, nb_cluster_c1 = 30, seuil = 100) :
	"""
	Fonction permettant d'associer à chaque cluster issu de gene_I50 un cluster domaine_Ifile
	input :
		- c1c2_file : fichier issu de la fonction classif1_in_classif2()
		- nb_cluster_c1 : nombre total de cluster avec c1. par defaut à 30 (correspondant à I = 5.0)
		- seuil : seuil au-dela de laquelle on peut considerer le cluster. par defaut à 100.
	output :
	"""
	#creation du dico contenant pour chaque cluster_C1 les clusters_C2
	c1_clusters = {}
	for i in range(1, int(nb_cluster_c1) + 1) :
		c1_clusters[i] = []
	#lecture du fichier
	filin = open(c1c2_file, 'r')
	l = filin.readline()
	for l in filin : 
		#print l
		c2_num = int(l.split("==>")[0])
		c1_res = l.split("==>")[1]
		#print c1_res
		for i in range(c1_res.count("]")) :
			c1_num = int(c1_res.split("]")[i].split("[")[0])
			taux = float(c1_res.split("]")[i].split("[")[1])
			#print c2_num, c1_num, taux
			if taux >= int(seuil) :
				c1_clusters[c1_num].append(c2_num)
	filin.close()
	#verif
	filout = open(c1c2_file + "_" + str(seuil) + "_domclusters_in_geneclusters.txt", "w")
	filout.write("N_gene_cluster\t[N_domaine_cluster]\n")
	for i in c1_clusters :
		if c1_clusters[i] != [] :
			#print i, c1_clusters[i]
			filout.write(str(i) + '\t' + str(c1_clusters[i]) + '\n')
	filout.close()

def occ_dom_by_domcluster(output_MCL, domain_list_file = "/data/mortaza/Lucifer/oprvsopr/domains/domain_list.txt") :
	"""
	Fonction permettant de donner pour chaque domaine les N° de cluster_dom dans lequels il apparaît et de donner également son occ (en nb ou pourcentage dans cluster)
	input :
		- output_MCL : fichier de sortie resultat MCL
		- domain_list_file : fichier contenant le nom de tous les domaines OPR
	output :
		- fichier txt contenant les resultats (pourcentage (filout2) ou occ (filout1))
	"""
	#dictionnaire des domaines contenant comme cle nom dom et comme valeur dico des clusters (ce dico aura pour cle le N° du cluster et pour valeur l'occ du dom dans le cluster considéré)
	DOMAINE = {}
	dom_file = open(domain_list_file, 'r')
	for l in dom_file :
		DOMAINE[l[:-1]] = {}
	dom_file.close()
	#traitement du fichier MCL
	mcl_file = open(output_MCL, 'r')
	#N° du cluster
	num_cluster = 0
	for l in mcl_file :
		#print l
		num_cluster += 1
		# separation par des tab ou espaces donne le nombre de mots (dom opr) sur la ligne = ds cluster
		for i in range(l.count('\t') + 1) :
			#l.split()[i] donne le nom complet du dom (collé au nom de la proteine OPR) enlever les artefacts (provenant du header)
			if l.split()[i] != "OPR1" and l.split()[i] != "OPR2" :	
				complete_name = l.split()[i]
				#print complete_name
				#nombre de separations par underscore
				#print complete_name.count('_')
				dom_name = ""
				#avoir seulement le nom du domaine
				for sep in range(2, complete_name.count('_') + 1) :
					dom_name = dom_name + complete_name.split('_')[sep] + "_"
				dom_name = dom_name[:-1]
				#print dom_name, num_cluster
				if num_cluster not in DOMAINE[dom_name] :
					DOMAINE[dom_name][num_cluster] = 1
				else :
					DOMAINE[dom_name][num_cluster] += 1
	#Verif des resultats du dico : DOMAINE = {dom_name : {num_cluster : occ_ds_cluster}}
	#filout1 = open(output_MCL + "_occ_dom_by_dom_cluster.txt", "w")
	#filout1.write("Domain_name\t>>>\tN°_Dom_Cluster [Occurence]\n")
	filout2 = open(output_MCL + "_perc_dom_by_dom_cluster.txt", "w")
	filout2.write("Domaine_name\t>>>\tN°_Dom_Cluster [Percentage]\n")
	for i in DOMAINE :
		#print str(i) + " >>> ",
		#filout1.write(str(i) + " >>> ")
		filout2.write(str(i) + " >>> ")
		#nombre total de domaines dans le cluster considéré
		total = 0.0
		for j in DOMAINE[i] : 
			total = total + DOMAINE[i][j]
		for j in DOMAINE[i] :
			#print str(j) + " [" + str(100 * DOMAINE[i][j] / total) + "] ",
			#filout1.write(str(j) + " [" + str(DOMAINE[i][j]) + "] ")
			filout2.write(str(j) + " [" + str(100 * DOMAINE[i][j] / total) + "] ")
		#print
		#filout1.write("\n")
		filout2.write("\n")
	mcl_file.close()
	#filout1.close()
	filout2.close()

def nb_cluster_by_domain(occ_dom_by_domcluster_rep) :
	"""
	Fonction permettant de donner le nombre de clusters dans lequel se trouve le domaine.
	Permet de voir si la composition des clusters est identique ou pas entre les diff I.
	input :
		- occ_dom_by_domcluster_rep : repertoire contenant seulement les fichiers .txt obtenu avec la fonction occ_dom_by_domcluster
	output :
		- 1 fichier txt avec 1 colonne pour les domaines et une colonne pour le nombre de clusters par I.
	"""
	#initialisation d'un dico contenant comme cle les noms des domaines et comme valeurs une liste 
	DOMAINE = {}
	files = glob.glob(occ_dom_by_domcluster_rep + '/*')
	print "Files used \n"
	for file in files :
		print file  #(pour voir l'ordre des I)
		f = open(file, 'r')
		for l in f :
			dom_name = l.split(">>>")[0]
			#print dom_name
			if dom_name not in DOMAINE :
				DOMAINE[dom_name] = []
			n_cluster = l.split(">>>")[1].count("[")
			DOMAINE[dom_name].append(n_cluster)
		f.close()
	print "\nDomain_name\tNb_of_clusters_containing_Domain_for_diff_I\n"
	for i in DOMAINE :
		print i, DOMAINE[i]


#programme principal
if __name__ == "__main__":

	#convert_blast_to_cytoscape('/home/mortaza/Documents/OPRvsOPR/opr_vs_opr_logEvalue.txt', '/home/mortaza/Documents/OPRvsOPR/OPR_genes_visual.txt')
	#convert_cyto_to_2csv_Gephi('/home/mortaza/Documents/OPRvsOPR/cyto.txt', '/home/mortaza/Documents/OPRvsOPR/nodes.csv', '/home/mortaza/Documents/OPRvsOPR/edges.csv')
	#convert_blast_to_cytoscape('/home/mortaza/Documents/OPRvsOPR/Cytoscape/opr_vs_opr_cytoscape.txt', '/home/mortaza/Documents/OPRvsOPR/opr_vs_opr.txt')
	#convert_blast_to_cytoscape('/data/mortaza/Lucifer/oprvsopr/genes/opr_vs_opr_cytoscape.txt', 'autre')
	#cluster_count('/data/mortaza/Lucifer/oprvsopr/genes/cluster1.csv')
	#cluster_count('/data/mortaza/Lucifer/oprvsopr/genes/mcl2.0-logEvalue.csv')
	#cluster_count('/data/mortaza/Lucifer/oprvsopr/genes/mcl2.5-logEvalue.csv')
	#cluster_count('/data/mortaza/Lucifer/oprvsopr/genes/mcl3.0-logEvalue.csv')
	#cluster_count('/data/mortaza/Lucifer/oprvsopr/genes/MCL_8_SCPS.csv')
	#*****************************
	#workflow pr OPR genes
	#convert_blast_to_cytoscape(sys.argv[1], sys.argv[2], sys.argv[3])
	#mcl_command(sys.argv[1], Istart = sys.argv[2], Iend = sys.argv[3], Istep = sys.argv[4])
	#mcl_cluster_count(sys.argv[1])
	#mcl_analysis(sys.argv[1], sys.argv[2]) #/home/mortaza/Documents/src/OPRvsOPR_rep/fonctions.py opr_vs_opr_0.001.txt out1
	#*****************************
	#workflow pr OPR domains
	#convert_blast_to_cytoscape(sys.argv[1], sys.argv[2], sys.argv[3])
	#mcl_command(sys.argv[1], Istart = sys.argv[2], Iend = sys.argv[3], Istep = sys.argv[4]) 
	#mcl_analysis(sys.argv[1], sys.argv[2]) 
	#*****************************
	#creation d'un fichier contenant les prot OPR et leur type
	#opr_type_in_file()
	#*****************************
	#analyse domaines OPR
	#1°) creer un fichier ayant tous les domaines existants
	#OPR_dom(sys.argv[1], sys.argv[2])
	#2°) 
	#analyse_MCL_dom_file(sys.argv[1]) 
	#classif1_in_classif2(sys.argv[1])
	#c1c2_analyse(sys.argv[1])
	#occ_dom_by_domcluster(sys.argv[1], domain_list_file = "/data/mortaza/Lucifer/oprvsopr/domains/domain_list.txt")
	#analyse_MCL_dom(sys.argv[1]) 
	nb_cluster_by_domain(sys.argv[1]) #/home/mortaza/Documents/src/OPRvsOPR_rep/fonctions.py type_dom/by_dom_clusters/by_percentage/ >> Clusters_Compostition_id_or_not.txt