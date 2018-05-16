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
from matplotlib import pyplot as plt
import numpy as np
sys.path.append('/data/mortaza/src/OPRvsgenome/')
from oprvsgenome import *
from db_building_part2 import *

#Fonctions

def dico_blast_results_pt(list_file, csv_files) :
	"""
	Fonction permettant de retourner les resultats de blast sous forme d'un 
	dico (decommenter les print pr afficher les resultats à l'écran)
	input :
		- list_file : fichier contenant la liste des omes (p ou t)
		- csv_files : ens des resultats blasts pour ce ome
	output :
		- retourne un dico de resultats de type
		ome, contig, sstart, send, strand, Evalue, query
	"""
	#initialisation dico
	dico = {}
	filin = open(list_file, "r")
	for line in filin :
		#pr chaque esp
		species = line[:-1]
		dico[species] = {}
	filin.close()
	#print dico
	#print csv_files
	for csv_file in csv_files :
		#print csv_file
		#espece
		sp = csv_file.split("/")[-1]
		#l'ens des resultats blasts du rep
		files = glob.glob(csv_file + "/*csv")
		#print files
		#pour chaque fichier csv de ce rep
		for file in files : 
			#print file
			filin = open(file, "r")
			for line in filin :
				#type de separateur 
				sep = line.split("\t")
				query = sep[0]
				contig = sep[4]
				sstart = int(sep[6])
				send = int(sep[7])
				qstart = int(sep[2])
				qend = int(sep[3])
				sframe = int(sep[9])
				sseq = sep[5].replace("-","").replace("*","")
				if sframe < 0 :
					strand = "-"
					if send < sstart :
						tmp = send
						send = sstart
						sstart = tmp
				else :
					strand = "+"
				Evalue = float(sep[16])
				if Evalue < 1e-06 :
					if contig not in dico[sp] :
						dico[sp][contig] = []
						dico[sp][contig].append((sstart, send, strand, Evalue, query, qstart, qend, sseq))
						#dico[sp][contig].append((sstart, send, strand, Evalue, query, qstart, qend))
						#dico[sp][contig].append((sstart, send, sframe, Evalue, query, sseq))
					else :
						dico[sp][contig].append((sstart, send, strand, Evalue, query, qstart, qend, sseq))
						#dico[sp][contig].append((sstart, send, strand, Evalue, query, qstart, qend))
						#dico[sp][contig].append((sstart, send, sframe, Evalue, query, sseq))
			filin.close()
	#Affichage des resultats
	#print "ome\tcontig\tsstart\tsend\tstrand\tEvalue\tquery\tsseq"
	for i in dico :
		#print i
		for j in dico[i] :
			#print j
			for k in dico[i][j] :
				#print i, j, k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]
				a = 1 # pr ne pas laisser la boucle vide (sans commandes)
	#utilisation du dico ds d'autres fct pr la suite de l'analyse
	return dico

def resume_regions_bis(list_opr_region) :
	"""
	Fonction permettant de retourner des coord de regions non chevauchantes
	input :
		- list_opr_region : liste de regions sous forme de tuples(chevauchantes
		 ou pas)
	output :
		- retourne une liste contenant des regions non chevauchantes
	"""
	rmin = 'X'
	rmax = 'X'
	rstrand = 'X'
	rEvalue = 20000
	rquery = []
	#les gdes regions dites opr
	real_regions = {}
	#initialisation nb de region
	nb_region = 1
	real_regions[nb_region] = (rmin, rmax, rstrand, rEvalue, rquery)
	#m = start, M = end, S = strand
	for (m, M, strand, Evalue, query) in list_opr_region :
		m = int(m); M = int(M)
		#min et max
		flag = 0
		#pour les regions dejà recensées
		for n in real_regions :
			rmin = real_regions[n][0]
			rmax = real_regions[n][1]
			rstrand = real_regions[n][2]
			rEvalue = real_regions[n][3]
			rquery = real_regions[n][4]
			#initialisation du rmin et rmax et rstrand avec les vraies valeurs
			if rmin == 'X' and rmax == 'X' and rstrand == 'X' :
				rmin = m
				rmax = M
				rstrand = strand
				flag = 1
			#sinon tester les differents cas
			elif m < rmin and rmin <= M <= rmax and strand == rstrand :
				rmin = m
				flag = 1
			elif M > rmax and rmin <= m <= rmax and strand == rstrand :
				rmax = M
				flag = 1
			elif m < rmin and M > rmax and strand == rstrand :
				rmin = m
				rmax = M
				flag = 1
			elif rmin <= m <= rmax and rmin <= M <= rmax and strand == rstrand:
				flag = 1
			if (flag == 0 and M < rmin) or (flag == 0 and m > rmax) or \
			(flag == 0 and strand != rstrand) :
				rmin = m
				rmax = M
				if strand != rstrand :
					rstrand = strand
			#print nb_region, rmin, rmax
			#si region chevauchante
			if flag == 1 :
				#choix Evalue la plus faible
				if float(Evalue) < float(rEvalue) :
					rEvalue = Evalue
				rquery.append(query)
				real_regions[n] = (rmin, rmax, rstrand, rEvalue, rquery)
			#print nb_region, m, M, rmin, rmax
		#sinon nouvelle region non chevauchante
		if flag == 0 :
				nb_region += 1
				#print nb_region, m, M, rmin, rmax
				rquery = [query]
				rEvalue = Evalue
				real_regions[nb_region] = (rmin, rmax, rstrand, rEvalue, rquery)
	#liste de tuples de coord de regions
	list_R = []
	for i in real_regions :
		#verif des differentes regions si pb
		#print i, real_regions[i]
		list_R.append(real_regions[i]) 
	return list_R

def dico_g_query(rep_resume1_files) :
	"""
	Fonction permettant de retourner les infos de ce fichier sous forme d'un 
	dico pr etre utilisable ds d'autres fct.
	input :
		- g_regions_file : fichier de resume des regions oprlike obtenues avec 
		le blast des proteines opr contre les génomes
	output :
		- retourne un dico avec regions de similarité à prendre en compte
	"""
	#recup dico
	dico_region = dico_info_resume1_file(rep_resume1_files)
	#initialisation dico pr infos nec
	dico_oprlike = {}
	for genome in dico_region :
		for contig in dico_region[genome] :
			#print dico_region[genome][contig]
			contig_length = int(dico_region[genome][contig][0])
			sstart = int(dico_region[genome][contig][1].split("-")[0])
			send = int(dico_region[genome][contig][1].split("-")[1])
			region_length = float(dico_region[genome][contig][2])
			#print contig, sstart, send, region_length
			dico_oprlike[contig] = (sstart, send)
	return dico_oprlike

def g_contig_selection(m, M, qstart, qend) :
	"""
	Fonction permettant de sélectionner le contig en fct d'une region de 
	similarité.
	input :
		- m : start de la region de similarité. entier.
		- M : end de la region de similarité. entier.
		- qstart : start du hit. entier.
		- qend : end du hit. entier.
	output :
		- retourne 0 (pas ok) ou 1 (ok)
	"""
	if qstart >= m and qend <= M :			return 1
	elif qstart < m and qend >= m + 100 :	return 1
	elif qend > M and qstart <= M - 100 :	return 1
	elif qstart < m and qend > M :			return 1
	else : 									return 0

def dico_g_query_result_analysis(list_file, csv_files, rep_resume1_files, type_ome = "X") :
	"""
	Fonction permettant de sélectionner les résultats génomiques 'corrects'. 
	si enregistrement des données, ne pas oublier de decommenter.
	input :
		- list_file : fichier contenant la liste des omes (p ou g ou t)
		- csv_files : ens des resultats blasts pour ce ome
		- rep_resume1_files : fichier de resume des regions oprlike obtenues avec 
		le blast des proteines opr contre les génomes
		- type_ome : p. par defaut "X" pr g ou t.
	output :
		- retourne un dico contenant les infos (comme dans la fonction 
	"""
	#recup dico des regions genomiques utilisées en tant que query
	dico_regions = dico_g_query(rep_resume1_files)
	#initialisation dico
	dico = {}
	#lecture du fichier des esp pour relever leur nom
	filin = open(list_file, "r")
	for line in filin :
		#initialisation dico de l'esp
		species = line[:-1]
		dico[species] = {}
	filin.close()
	#prise en compte des fichiers resultats de blast dans les rep csv_files
	for csv_file in csv_files :
		#espece
		sp = csv_file.split("/")[-1]
		#l'ens des resultats blasts du rep de l'esp sp
		files = glob.glob(csv_file + "/*csv")
		#pour chaque fichier csv de ce rep
		for file in files : 
			#lecture du fichier
			filin = open(file, "r")
			#de chaque ligne
			for line in filin :
				#type de separateur de ligne
				sep = line.split("\t")
				#infos de la ligne
				query = sep[0]
				contig = sep[4]
				qstart = int(sep[2])
				qend = int(sep[3])
				sstart = int(sep[6])
				send = int(sep[7])
				sframe = int(sep[9])
				sseq = sep[5]
				#de facon a avoir qstart < qend
				if qend < qstart :
					tmp = qend
					qend = qstart
					qstart = tmp
				#de facon a avoir sstart < send et avoir le strand
				if sframe < 0 :
					strand = "-"
					if send < sstart :
						tmp = send
						send = sstart
						sstart = tmp
				else :
					strand = "+"
				Evalue = float(sep[16])
				#selection de Evalue pour les proteines 
				flag = 0
				if type_ome == "p" and Evalue < 1e-03 :
					flag = 1
				elif (type_ome == "t" or type_ome == "g") and Evalue <= 1e-09 :
					flag = 1
				if flag == 1 :
					#selection des hits sur la region de similarité génomique
					m = dico_regions[query][0] ; M = dico_regions[query][1]
					if g_contig_selection(m, M, qstart, qend) == 1 :
						#enregistrement des hits sélectionnés ds le dico
						if contig not in dico[sp] :
							dico[sp][contig] = []
							#dico[sp][contig].append((sstart, send, strand, Evalue, query, sseq))
							dico[sp][contig].append((sstart, send, sframe, Evalue, query, sseq))
						else :
							#dico[sp][contig].append((sstart, send, strand, Evalue, query, sseq))
							dico[sp][contig].append((sstart, send, sframe, Evalue, query, sseq))
			filin.close()
	#Affichage des resultats
	#print "ome\tcontig\tsstart\tsend\tstrand\tEvalue\tquery\tsseq"
	for i in dico :
		#print i
		for j in dico[i] :
			#print j
			for k in dico[i][j] :
				print i, j, k[0], k[1], k[2], k[3], k[4], k[5]
				a = 1 # pr ne pas laisser la boucle vide (sans commandes)
	#utilisation du dico ds d'autres fct pr la suite de l'analyse
	return dico

def cleaned_opr_like_regions(list_file, contig_len_rep, csv_files_p_g, csv_files_p_pt, rep_resume1_files) :
	"""
	Fonction permettant de retourner les regions OPR_like non chevauchantes
	input :
		- list_file : fichier contenant la liste des omes (p ou g ou t)
		- contig_len_rep : chemin du rep contenant la taille des contigs
		- csv_files_p_g : ens des resultats blasts from g query
		- csv_files_p_pt : ens des resultats blasts from p or t query
		- rep_resume1_files : fichier de resume des regions oprlike obtenues avec 
		le blast des proteines opr contre les génomes
	output :
		- retourne les regions non chevauchantes
	"""
	dico_contigs = dico_len_contigs(contig_len_rep)
	#verif dico
	#for i in dico_contigs :
	#	for j in dico_contigs[i] :
	#		print i, j, dico_contigs[i][j]
	#print dico_contigs["Scenedesmus_obliquus_6623EST"]
	#header
	print "species\tcontig\tstart\tend\tstrand\tEvalue\tcouverture_contig\tquery"
	#dico resumant les resultats des fichiers csv / de sortie de blast
	#resultats p et t
	all_results_dico1 = dico_blast_results_pt(list_file, csv_files_p_pt)
	#resultats genomiques
	ome = contig_len_rep.split("/")[-1][0]
	all_results_dico2 = dico_g_query_result_analysis(list_file, csv_files_p_g,\
		rep_resume1_files, ome)
	l = [all_results_dico1, all_results_dico2]
	#l = [all_results_dico1]
	for all_results_dico in l :
		#pr chaque esp
		for species in all_results_dico :
			#print species
			#nom de l'esp pr pouvoir acceder aux données de dico_contigs
			sp_contig = species.split(".")[0]
			#print sp_contig
			#pr chaque contig 
			for contig in all_results_dico[species] :
				#print contig, sp_contig
				#modif au niv ttt des prot
				#if "|" in contig :
				#	ctg = contig.split("|")[-2]
				#else :
				#	ctg = contig
				ctg = contig
				regions_list = sorted(all_results_dico[species][contig], \
					key=lambda data: data[0])
				opr_list = resume_regions_bis(regions_list)
				#print opr_list
				for i in opr_list :
					#couverture du  contig
					cov = round(100 * (i[1] - i[0]) / float(dico_contigs[sp_contig][ctg]), 2)
					#cov = "X"
					print species + "\t" + contig + "\t" + \
					str(i[0]) + "\t" + str(i[1]) + "\t" + i[2] + \
					"\t" + str(i[3]) + "\t" + str(cov) + "\t", 
					print i[4]# + "\t" + i[4]				

def dico_oprlike(oprlike_file) :
	"""
	Fonction permettant de renvoyer les infos du fichier ss forme de dico
	input :
		- oprlike_file : fichier resultat obtenu avec la fct 
		cleaned_opr_like_regions() de type
	output : 
		- retourne un dico
	"""
	#initialisation dico 
	dico = {}
	#lecture du fichier pr mettre les données ds un dico
	filin = open(oprlike_file, "r")
	#ne pas lire le header
	line = filin.readline()
	#pr les autres lignes du fichier :
	for line in filin : 
		#print line
		#infos
		sep = line.split("\t")
		ome = sep[0]
		if ome not in dico : 
			dico[ome] = {}
		contig = sep[1]
		if contig not in dico[ome] :
			dico[ome][contig] = []
		sstart = int(sep[2])
		send = int(sep[3])
		strand = sep[4]
		Evalue = float(sep[5])
		cov = sep[6]
		#avoir sous forme de liste l'ens des query
		query = []
		for i in range(sep[7].count(",") + 1) :
			q = sep[7].split(",")[i].replace("'","").replace(" ","")\
			.replace("[","").replace("]\n","")
			query.append(q)
		#enregistrement des données dans le dico
		dico[ome][contig].append((sstart, send, strand, Evalue, cov, query))
	filin.close()
	return dico

def selection_by_opr_len(oprlike_file, contig_len_rep, th = 76, th_cov = 50) :
	"""
	Fonction permettant de selectionner les proteines oprlike repondant aux 
	criteres de repetitions (1 repet OPR = 38 res = 96 nt). que pr analyse pt.
	selection par couverture aussi. Selection par la Evalue (seuil etait à 1e-03,
	maintenant à 1e-06).
	input :
		- oprlike_file : fichier resultat obtenu avec la fct 
		cleaned_opr_like_regions() de type
		species|contig|start|end|strand|Evalue|cov|query
		- contig_len_rep = rep de taille des contigs
		- th : seuil pr selectionner les vraies opr. par defaut, 
		4 rep de 32 res.
	output : 
		- affiche les resultats de la selection
	"""
	#header
	print "species\tcontig\tstart\tend\tstrand\tEvalue\tcouverture_contig\tquery"
	#recup dico de taille des contigs
	dico_contigs = dico_len_contigs(contig_len_rep)
	#conversion en entier
	th = int(th)
	th_cov = int(th_cov)
	#recuperation dico resumant les infos du fichier
	dico = dico_oprlike(oprlike_file)
	#pr chaque esp
	for ome in dico :
		sp = ome.split(".")[0]
		#pr chaque contig
		for contig in dico[ome] :
			if "|" in contig :
				ctg = contig.split("|")[-2]
			else :
				ctg = contig
			taille = 0; rEvalue = 2000; rstrand = "X"; rquery = []
			for region in dico[ome][contig] :
				#releve des infos à partir de la region
				sstart = region[0]; send = region[1]; strand = region[2]
				Evalue = region[3]; cov = region[4]; query = region[5]
				#mise à jour des variables
				taille += send - sstart
				if Evalue < rEvalue : rEvalue = Evalue
				rquery += query
				rquery = list(set(rquery))
				if rstrand == "X" :
					rstart = sstart
					rend = send
					rstrand = strand
				else :
					if sstart < rstart : rstart = sstart
					if send > rend : rend = send
					if strand != rstrand : rstrand = "?" 
			#selection par la taille de la region
			if taille >= th :
				rcov = round(100.0 * (rend - rstart) / dico_contigs[sp][ctg], 2)
				if rcov >= th_cov :
					if rEvalue < float(1e-06) :
						print ome + "\t" + contig + "\t" + str(rstart) + "\t" + \
						str(rend) + "\t" + rstrand + "\t" + str(rEvalue) + "\t" + \
						str(rcov) + "\t", 
						print rquery

def dico_region_blast1(blast1_file) :
	"""
	Fonction permettant de retourner un dico des regions issues du blast1. les 
	verif donneront le nb de oprlike pr chaque ome.
	input : 
		- blast1_file : fichier resumé des régions des omes considérés
	output :
		- retourne un dico {ome : [contig]}
	"""
	#initialisation du dico 
	dico_blast1 = {}
	filin = open(blast1_file, "r")
	for line in filin :
		#initialisation de la liste de regions de l'ome 
		if ">" in line :
			ome = line[2:-1]
			dico_blast1[ome] = []
		#completer la liste des regions de cet ome
		else : 
			region = line[:-1]
			dico_blast1[ome].append(region)
	filin.close()
	#verif dico
	for i in dico_blast1 :
		print i
		print len(dico_blast1[i])
	return dico_blast1

def dico_region_blast2(blast2_file) :
	"""
	Fonction permettant de retourner un dico des regions issues du blast2. les
	verif donneront le nb de oprlike pr chaque ome. modif du sep "\t" ou " ".
	input :
		- blast2_file : fichier issu de la fct selection_by_opr_len()
	output :
		- retourne un dico {ome : [contig]}
	"""
	#initialisation du dico
	dico_blast2 = {}
	filin = open(blast2_file, "r")
	#ne pas prendre en compte le header
	line = filin.readline()
	for line in filin :
		#print line
		ome = line.split(" ")[0]
		region = line.split(" ")[1].replace("\n", "")
		if "|" in region :
			region = region.split("|")[-2]
		if ome not in dico_blast2 : 
			dico_blast2[ome] = []
		dico_blast2[ome].append(region)
	filin.close()
	#verif dico
	#for i in dico_blast2 :
	#	print i
	#	print len(dico_blast2[i])
		#print dico_blast2[i]
	return dico_blast2

def dico_region_blast2_FTRep(FTRep_OPRcandidates_rep, list_file) :
	"""
	Fonction permettant de retourner dans un dico l'ens des regions pr chaque
	ome.
	input :
		- FTRep_OPRcandidates_file : rep des fichiers obtenus à partir de la 
		fct OPR_candidates_name_from_FTRep_result()
		- list_file : fichier contenant la liste des proteomes
	output :
		- retourne un dico
	"""
	#creation dico vide ({ome : []})
	dico = {}
	filin1 = open(list_file, "r")
	for line1 in filin1 : 
		dico[line1[:-1]] = []
	filin1.close()
	#print dico
	#ens des fichiers contenant les OPRcandidates
	files = glob.glob(FTRep_OPRcandidates_rep + "/*OPRcandidates")
	for file in files :
		ome = file.split("/")[-1].split("_OPRcandidates")[0].split("FTRep_")[1]
		for keys in dico :
			if ome in keys :
				#print ome
				filin = open(file, "r")
				for line in filin :
					dico[keys].append(line[:-1])
				filin.close()
	#for i in dico :
	#	print i
	#	print dico[i]
	return dico

def blast1_blast2_comparison(blast1_file, blast2_file, list_file) :
	"""
	Fonction permettant de comparer les regions issues du blast1 et celles 
	issues du blast2. affiche les regions du blast1 non retrouvées parmi les
	regions du blast2.
	input : 
		- blast1_file : fichier resumé des régions des omes considérés
		- blast2_file : fichier issu de la fct selection_by_opr_len()
		- list_file : fichier contenant la liste des omes
	output :
		- affiche les regions du blast1 non retrouvées dans blast2
	"""
	#recup des regions oprlike
	blast1 = dico_region_blast1(blast1_file)
	#blast2 = dico_region_blast2(blast2_file)
	blast2 = dico_region_blast2_FTRep(blast2_file, list_file)
	#verif des dico et permet aussi de donner le nb de regions oprlike
	#issues des deux blasts
	#for i in blast1 :
	#	print i
	#	print len(blast1[i]),
	#	print len(blast2[i])
	#header
	print "ome\tnb_total_blast2\tnb_regions_found\tnb_regions_not_found\tregions_not_found"
	#pr l'esp
	for ome in blast1 :
		#print "\nome : ", ome
		#print "regions not found : ", 
		#nb total de contig non retrouvé
		t = 0
		#liste contenant les contigs non retrouvés
		l = []
		#pr le contig / la region du blast1
		for contig in blast1[ome] :
			if contig not in blast2[ome] :
				t += 1
				l.append(contig)
		#print "\nnb regions not found : ", t
		#print "nb regions found : ", len(blast1[ome]) - t
		print ome + "\t" + str(len(blast2[ome])) + "\t" + str(len(blast1[ome]) - t) + "\t" + str(t) + "\t",
		print l
	#print blast2["GCA_001584585.1_ASM158458v1_G_pecorale_protein.faa"]

def Evalue_blast1_vs_blast2(ome1_oprlike, ome2_oprlike) :
	"""
	Fonction permettant de comparer les Evalues des regions retrouvées dans les
	resultats de blast1 et de blast2.
	input :
		- proteome1 : fichier csv resumant les resultats de blast1. créé à 
		partir de la fonction cleaned_opr_like_regions()
		- proteome2 : fichier csv resumant les resultats de blast2. créé à
		partir de proteome1.csvla fonction cleaned_opr_like_regions()
	output :
		- affiche les resultats de la comparaison
	"""
	#dico qui contiendra les valeurs de Evalues du blast1 et du blast2
	#dico_comparison = {}
	#dico des ome de blast1 et de blast2
	dico1 = dico_oprlike(ome1_oprlike)
	dico2 = dico_oprlike(ome2_oprlike)
	#header
	print "ome\tcontig\tEvalue1\tEvalue2"
	#nb de Evalue1 <= Evalue2
	nb1 = 0
	#nb de Evalue1 > Evalue2
	nb2 = 0
	for ome in dico1 :
		#print ome
		#dico_comparison[ome] = {}
		for contig in dico1[ome] :
			#print contig
			#if contig not in dico_comparison[ome] : 
				#dico_comparison[ome][contig] = {"Evalue1" : [], "Evalue2" : []}
			for region1 in dico1[ome][contig] :
				evalue1 = float(region1[3])
				#dico_comparison[ome][contig]["Evalue1"].append(evalue1)
			for region2 in dico2[ome][contig] :
				evalue2 = float(region2[3])
				#dico_comparison[ome][contig]["Evalue2"].append(evalue2)
			print ome, contig, evalue1, evalue2
			if evalue1 <= evalue2 :
				nb1 += 1
			else : 
				nb2 += 2
	print nb1, nb2

def distr_Evalue(rep_ome) :
	"""
	Fonction permettant de donner la distribution des Evalues. 
	input : 
		- rep_ome : repertoire contenant les repertoires des diff omes avec
		les resultats csv de blast dans chacun de ces rep
	output : 
		- plot de la distribution des Evalues. 
	"""
	#liste contenant toutes les valeurs de Evalues
	Evalue = []
	rep_omes = glob.glob(rep_ome + "/*f*")
	for ome in rep_omes :
		files = glob.glob(ome + "/*csv")
		for file in files : 
			filin = open(file, "r")
			for l in filin :
				evalue = float(l.split("\t")[16])
				#print evalue
				if evalue == 0.0 :
					evalue = 1e-200
				if evalue <= 1e-06 :
					Evalue.append(-log(evalue))
			filin.close()
	#tracer le plot de la distribution des Evalues
	plt.hist(np.array(Evalue), bins = 100)
	plt.title("Evalue Distribution")
	plt.xlabel("-log(Evalue)")
	plt.ylabel("Number of Evalue")
	#pour avoir meme echelle sur l'axe des ordonnées
	plt.ylim(0, 11300)
	plt.show()

def multifasta_per_ome_for_FTRep(selected_opr_file, rep_fasta_files) :
	"""
	Fonction permettant de créer des fichiers multifasta par espece.
	input : 
		- selected_opr_file : fichier csv obtenu avec la fonction 
		selection_by_opr_len()
		- rep_fasta_files : repertoire contenant l'ens des fichiers fasta d'un
		ome
	output : 
		- création de fichiers multifasta.
	"""
	dico = dico_region_blast2(selected_opr_file)
	dico_fasta = dico_info_fasta(rep_fasta_files)
	for ome in dico :
		sp = "FTRep_" + ome
		print sp
		filout = open(sp, "w")
		for contig in dico[ome] :
			#print contig
			if ome in dico_fasta :
				seq = dico_fasta[ome][contig] + "\n"
				#print ">" + contig + "\n",
				#print seq, 
				filout.write(">" + contig + "\n")
				filout.write(seq)
		filout.close()

def OPR_candidates_name_from_FTRep_result(FTRep_result) :
	"""
	Fonction qui donne le nom des OPR candidates à partir d'un resultat de FTRep
	input : 
		- FTRep_result : repertoire contenant les resultats de FTRep (file.res)
	output :
		- affichage des noms de OPR candidates
	"""
	#recuperation de tous les fichiers .res dans une liste
	files = glob.glob(FTRep_result + "/*.res")
	#pour chacun des fichier de resultat
	for file in files :
		ome = file.split("/")[-1][:-4] + "_OPRcandidates"
		print ome
		filin = open(file, "r")
		#ne pas lire le header du fichier
		line = filin.readline()
		#liste contenant tous les opr
		opr_list = []
		#pour les autres lignes, relever les noms des OPR candidates
		for line in filin : 
			opr = line.split("\t")[0]
			test = line.split("\t")[2]
			#si contient des repet OPR, alors c'est un OPR candidates
			if test == "OPR" :
				if opr not in opr_list : 
					opr_list.append(opr)
		filin.close()
		#stockage des resultats ds un fichier
		filout = open(ome, "w")
		for opr in opr_list :
			filout.write(opr + "\n")
		filout.close()




#Programme principal
if __name__ == "__main__":
	#csv_rep = glob.glob("/data/mortaza/Lucifer/evolution/oprvsgenome/results/proteome/*f*")
	#distr_Evalue(sys.argv[1])
	#les deux lignes de commandes à lancer ds le rep distr_Evalues
	#/home/mortaza/Documents/src/OPRvsgenome/result_analysis_part2.py ../../../QueryP/SubjectP
	#/home/mortaza/Documents/src/OPRvsgenome/result_analysis_part2.py /data/mortaza/Lucifer/evolution/oprvsgenome/results/proteome
	#suppression des fichiers vides ds les rep de resultats de blast
	#frame_fichiers_non_vide(sys.argv[1])
	#rangement des fichiers csv de chaque esp dans des rep ầ part
	#results_stockage(sys.argv[1])
	#analyse proteomes
	##regions p et t
	csv_files_p_pt = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryP/SubjectP/*') + \
	glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryT/SubjectP/*')
	##si enregistrement des données
	#dico_blast_results_pt(sys.argv[1], csv_files_p_pt) 
	##regions genomiques
	csv_files_p_g = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryG/SubjectP/*')
	rep_resume1_files = "/data/mortaza/Lucifer/evolution/oprvsgenome/results/genome/regions/resume/part1"
	##si enregistrement des données
	##dico_g_query_result_analysis(sys.argv[1], csv_files_p_g, rep_resume1_files, "p")
	#obtention de regions non chevauchantes
	#cleaned_opr_like_regions(sys.argv[1], sys.argv[2], csv_files_p_g, \
	#	csv_files_p_pt, rep_resume1_files)
	#selection des oprlike en fct du nb de rep (en nt ou res) et de la couverture
	#selection_by_opr_len(sys.argv[1], sys.argv[2])
	#verif que regions blast1 dans blast2
	#blast1_blast2_comparison(sys.argv[1], sys.argv[2])
	#Creation de fichiers multifasta pr FTRep
	#multifasta_per_ome_for_FTRep(sys.argv[1], sys.argv[2])
	#Analyse des resultats de FTRep
	#OPR_candidates_name_from_FTRep_result(sys.argv[1])
	#blast1 dans FTRep ?
	#blast1_blast2_comparison(sys.argv[1], sys.argv[2], sys.argv[3])
	#analyse transcriptomes
	csv_files_t_t = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryT/SubjectT/*')
	csv_files_t_g = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryG/SubjectT/*')
	#lignes de commandes
	#cleaned_opr_like_regions(sys.argv[1], sys.argv[2], csv_files_t_g, csv_files_t_t, rep_resume1_files) 

"""
Notes & Commentaires :
- /data/mortaza/src/OPRvsgenome/result_analysis_part2.py
- liste de notre db perso (g, t, p) ds le rep
/data/mortaza/Lucifer/evolution/OPR_evolution/oprvsgenome/results/list
- chemin pour avoir les rep des tailles des diff contigs
/data/mortaza/Data
"""