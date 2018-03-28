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
sys.path.append('/home/mortaza/Documents/src/OPRvsOPR_rep/')
from cluster_quality import *

#Fonctions

def results_stockage(species_file) :
	"""
	Fonction permettant de ranger les résultats obtenus avec tblastn.
	input :
		- species_file : fichier contenant la liste des fichiers fasta des génomes
	output :
		- rangement des fichiers et premieres analyses du contenu des fichiers (nb de lignes)
	"""
	filin = open(species_file, 'r')
	S = []
	for l in filin :
		S.append(l[:-1])
	filin.close()
	#print len(S)
	#print S
	"""
	for espece in S :
		print espece
		os.chdir(espece)
		os.system("wc -l *csv")
		os.chdir("..")
		#os.system("ls " + espece + " | wc -l")
	"""
	#creer repertoire avec les csv repartis
	files = glob.glob("./*.csv")
	for espece in S :
		#print espece
		os.system("mkdir " + espece)
		for file in files :
			if espece in file :
				#print file
				os.system("mv " + file + " " + espece)
	
def frame_fichiers_non_vide(species_rep) :
	"""
	Fonction permettant d'afficher les fichiers qui contiennent au moins un résultat.
	Suppression des fichiers sans resultats
	input : 
		- species_rep : repertoire contenant les resultats de tblastn pour une espece donnée
	output :
		- fichiers texte par OPR et par frame créés
	"""
	#ens des fichiers tblastn format csv presents dans ce rep
	files = glob.glob(species_rep + "/*.csv")
	#pour chacun de ces fichiers
	for file in files :
		#print file
		#stocker toutes les lignes du fichiers dans une liste
		filin = open(file, 'r')
		lines = filin.readlines()
		filin.close()
		#si le fichier de resultat n'est pas vide (taille de la liste diff de 0)
		if len(lines) == 0 :
			os.system("rm " + file)

def genome_quality(genome_file) :
	"""
	Fonction permettant de donner la distribution de la taille des contigs du genome en fonction de leur nombre
	input : 
		- genome_file : fichier fasta du  genome
	output :
	"""
	#initialisation de la taille
	contig = 0
	length_contigs = [] #taille de cette liste donnera le nombre de contigs
	fasta = open(genome_file, 'r')
	for line in fasta :
		if ">" in line :
			if contig != 0 : # pour ne pas mettre le compteur initialisé
				length_contigs.append(contig)
			contig = 0
		if ">" not in line :
			flag = 1 #presence de seq
			contig = contig + len(line[:-1])
	#Ajout de la taille du dernier contig
	length_contigs.append(contig)
	mediane = np.median(np.array(length_contigs))
	quartile_25 = np.percentile(np.array(length_contigs), 25)
	quartile_75 = np.percentile(np.array(length_contigs), 75)
	#print mediane, quartile_25, quartile_75
	#print len(length_contigs)
	plt.hist(np.array(length_contigs), bins = 100)
	plt.axvline(mediane, color = 'red', label = str(mediane))
	plt.axvline(quartile_25, color = 'green', label = str(quartile_25))
	plt.axvline(quartile_75, color = 'green', label = str(quartile_75))
	plt.title("Distribution of Genome Contigs length")
	plt.xlabel("Genome Contigs Length")
	plt.ylabel("Number of Genome Contigs")
	plt.legend()
	#plt.show()
	plt.savefig(genome_file + '.png')
	plt.close()
	fasta.close()

def dico_len_OPR(OPR_length_file) :
 	"""
 	Fonction permettant de retourner un dico avec cle = OPR et valeur = taille de l'OPR
 	input :
 		- OPR_length_file : fichier txt contenant les OPR et leurs tailles
 	output :
 		- retourne un dico
 	"""
 	OPR_qlen = {}
 	filin = open(OPR_length_file, 'r')
 	for line in filin :
 		OPR_qlen[line.split()[0]] = line.split()[1]
 	filin.close()
 	#verif
 	#for i in OPR_qlen :
 	#	print i, OPR_qlen[i]
 	return OPR_qlen

def OPR_regions(species_rep, OPR_length_file, contig_len_rep) :
	"""
	Fonction permettant de donner les differentes regions dites OPR d'une espece à partir des resultats de tblastn
	input :
		- species_rep : repertoire contenant les resultats de tblastn pour une espece donnée
		- OPR_length_file : fichier txt contenant les OPR et leurs tailles
		- contig_len_rep : repertoire contenant les fichiers (contigs - taille)
	output :
		- fichier resumé de toutes les regions OPR de l'espece (qstart, qend, sstart, send). 
		on considere que frame1=frame2=frame3 et framem1=framem2=framem3
	"""
	###pg gal
	#dico avec cle = OPR et valeurs = liste des regions dans l'espece (qstart, qend, sstart, send)
	OPR_regions = {}
	#dico avec cle = contig et valeurs = liste des regions OPR
	contig_regions = {}
	#ens des fichiers tblastn format csv presents dans ce rep
	files = glob.glob(species_rep + "/*.csv")
	#pour chacun de ces fichiers
	for file in files :
		#print file
		#stocker toutes les lignes du fichiers dans une liste
		filin = open(file, 'r')
		for line in filin :
			OPR = line.split()[0]
			qstart = int(line.split()[2])
			qend = int(line.split()[3])
			contig = line.split()[4]
			sstart = int(line.split()[6])
			send = int(line.split()[7])
			sframe = int(line.split()[9])
			#print OPR, qstart, qend, contig, sstart, send, sframe
			collantes = ["Chlre_OPR105", "Chlre_OPR107", "Chlre_OPR28", "Chlre_OPR29", "Chlre_OPR32", "Chlre_OPR62"]
			#ncl = ["Chlre_OPR110", "Chlre_OPR86", "Chlre_OPR89", "Chlre_OPR97"] #proteines NCL pour genomes et transcriptomes
			ncl = ["NCL", "NCL1", "NCL2"] #profils NCL pour proteomes
			#print collantes
			#print ncl
			#if (OPR not in collantes) and (OPR not in ncl):
			if 1 :
				#print OPR
				#####Contigs
				#initialisation dico pour contig en question
				if contig not in contig_regions :
					contig_regions[contig] = {}
				if OPR not in contig_regions[contig] :
					contig_regions[contig][OPR] = []
				#normalement contig et OPR connu, ajout des infos ds la liste
				contig_regions[contig][OPR].append((qstart, qend, sstart, send, sframe))
				#####OPRs
				#initialisation dico pour OPR en question
				if OPR not in OPR_regions :
					OPR_regions[OPR] = {}
				#initialisation de la liste pour le contig concernant l'OPR en question
				if contig not in OPR_regions[OPR] :
					OPR_regions[OPR][contig] = []
				
				#normalement OPR et contig connu, ajout des infos dans la liste
				OPR_regions[OPR][contig].append((qstart, qend, sstart, send, sframe))
		filin.close()
	#soit creer des fichiers de resultats
	#info_OPRs(OPR_regions, OPR_length_file)
	#print species_rep + "_opr_regions"
	#info_contigs(contig_regions, OPR_length_file, contig_len_rep, species_rep + "_opr_regions")
	#soit utiliser seulement le dico
	return OPR_regions

def info_contigs(contig_regions, OPR_length_file, contig_len_rep, output) :
	"""
	Une partie de la fonction OPR_regions.
	Affichage :
	> Contig
	OPR : sstart send --> % de couverture de l'OPR
	input :
		- contig_regions : dico de la fct generale
		- OPR_length_file : fichier contenant la taille des diff OPRs
		- contig_len_rep : repertoire contenant les fichiers (contigs - taille)
	#verif que pas de pb avec dico
	for contig in contig_regions :
		print contig
		for OPR in contig_regions[contig] :
			print OPR, contig_regions[contig][OPR]
	"""
	len_contig = dico_len_contigs(contig_len_rep)
	OPR_qlen = dico_len_OPR(OPR_length_file)
	#dico contenant les diff regions opr possibles
	opr_regions = {}
	#pour chaque contig
	for contig in contig_regions :
		#print ">", contig
		#ajout contig dans le dico
		opr_regions[contig] = []
		#pour chaque opr
		for opr in contig_regions[contig] :
			#print opr
			#flag pour initialisation des smin et smax 
			flag = 0
			liste = contig_regions[contig][opr]
			#min et max du query (OPR)
			qmin = 100000.0
			qmax = 0.0
			#pour chaque info de la liste :
			for i in liste :
				if i[0] < qmin :
					qmin = i[0]
				if i[1] > qmax :
					qmax = i[1]
				if int(i[4]) > 0 :
					#print "+"
					if flag == 0 :
						smin = int(i[2])
						smax = int(i[3])
						flag = 1
					if int(i[2]) < smin :
						smin = int(i[2])
					if int(i[3]) > smax :
						smax = int(i[3])
				else :
					#print "-"
					if flag == 0 :
						smin = int(i[3])
						smax = int(i[2])
						flag = 1
					if int(i[3]) < smin :
						smin = int(i[3])
					if int(i[2]) > smax :
						smax = int(i[2])
			#print opr, ":", smin, smax, "-->", 100 * (qmax - qmin) / float(OPR_qlen[opr]), "-", qmin, qmax
			#ajouter les infos dans le dico opr_regions
			opr_regions[contig].append((opr, smin, smax))
	filout = open(output, 'w')
	#verif du dico opr_regions {contig : [(opr, smin, smax)]} #affichage des grandes regions opr
	for contig in opr_regions :
		#print contig
		#filout.write(str(contig) + "\n")
		#initialisation de rmin et rmax (r = regions)
		rmin = 'X'
		rmax = 'X'
		#les gdes regions dites opr
		real_regions = {}
		#initialisation nb de region
		nb_region = 1
		real_regions[nb_region] = [rmin, rmax, []]
		for region in opr_regions[contig] :
			#print region
			#min et max
			opr = region[0]
			m = int(region[1])
			M = int(region[2])
			flag = 0
			#pour les regions dejà recensées
			for n in real_regions :
				rmin = real_regions[n][0]
				rmax = real_regions[n][1]
				opr_list = real_regions[n][2]
				#initialisation du rmin et rmax avec les vraies valeurs
				if rmin == 'X' and rmax == 'X' :
					rmin = m
					rmax = M
					flag = 1
				elif m < rmin and rmin <= M <= rmax :
					rmin = m
					flag = 1
				elif M > rmax and rmin <= m <= rmax :
					rmax = M
					flag = 1
				elif m < rmin and M > rmax :
					rmin = m
					rmax = M
					flag = 1
				elif rmin <= m <= rmax and rmin <= M <= rmax :
					flag = 1
				if (flag == 0 and M < rmin) or (flag == 0 and m > rmax) :
					rmin = m
					rmax = M
				#print nb_region, rmin, rmax
				if flag == 1 :
					if opr not in opr_list :
						opr_list.append(opr)
					real_regions[n] = [rmin, rmax, opr_list]
				#print nb_region, m, M, rmin, rmax
			if flag == 0 :
					nb_region += 1
					#print nb_region, m, M, rmin, rmax
					opr_list = [opr]
					real_regions[nb_region] = [rmin, rmax, opr_list]
		for i in real_regions :
			#print contig, i, real_regions[i][0], real_regions[i][1], real_regions[i][2]
			filout.write(str(contig) + "\t" + str(i) + "\t" + str(real_regions[i][0]) + "\t" + str(real_regions[i][1]) + "\t" + str(real_regions[i][2]) \
				+ "\t" + str(int(real_regions[i][1]) - int(real_regions[i][0])) + "\n")
	filout.close()
	
def info_OPRs(OPR_regions, OPR_length_file) :
	"""
	Une partie de la fonction OPR_regions. 
	Affichage :
	> OPR - taille de l'OPR
	contig --> %age de couverture qmin qmax
	"""
	#verif
	#for OPR in OPR_regions :
	#	print OPR
	#	for contig in OPR_regions[OPR] :
	#		for i in range(len(OPR_regions[OPR][contig])) :
	#			print contig, OPR_regions[OPR][contig][i]
	
	OPR_qlen = dico_len_OPR(OPR_length_file)
	for opr in OPR_regions :
		print ">", opr, '-', OPR_qlen[opr]
		for region in OPR_regions[opr] :
			#print region
			#regions rangées selon le query 
			liste = sorted(OPR_regions[opr][region], key=lambda data: data[0]) 
			 #print liste
			#min et max du query (OPR)
			qmin = 100000.0
			qmax = 0.0
			for i in liste :
				if i[0] < qmin :
					qmin = i[0]
				if i[1] > qmax :
					qmax = i[1]
				#print de la liste
				#print i
				print region, "-->", 100 * (qmax - qmin) / float(OPR_qlen[opr]), qmin, qmax, i[2], i[3], i[4]

def analysis_all_files(species_file, results_rep, OPR_length_file = '/data/mortaza/Lucifer/oprvsopr/genes/qcov_evalue_selection/OPR_qlen.txt', contig_len_rep = '/data/mortaza/Data/genomes/Contigs_length') :
 	"""
 	Fonction permettant de traiter les resultats de plusieurs rep (genomes)
 	input : 
 		- species_file : fichier contenant la liste des fichiers fasta des génomes
 		- results_rep : repertoire contentant les repertoires des diff resultats tblastn des genomes
 		- OPR_length_file = fichier de la taille des diff OPRs
 		- contig_len_rep : repertoire contenant les fichiers (contigs - taille)
 	output : 
 		- le resultat (des fichiers frames)
 	"""
 	#lecture de la liste des genomes disponibles (rmq : les diff rep ont le meme nom)
 	filin = open(species_file, 'r')
	S = []
	for l in filin :
		S.append(l[:-1])
	filin.close()
	#print len(S)
	#print S
	#pour chaque esp 
	for espece in S :
		#chemin du rep
		species_rep = results_rep + "/" + espece
		print species_rep
		#creation des fichiers frames
		#frame_fichiers_non_vide(species_rep) 
		#OPR_regions(species_rep, OPR_length_file, contig_len_rep)
		region_analysis(results_rep, species_rep, OPR_length_file, contig_len_rep)
		
def contigs_length(species_fasta_rep) :
	"""
	Fonction permettant de créer un fichier par genome et contenant le nom et la taille des diff contigs
	input :
		- species_fasta_rep : nom du rep contenant les fichiers fasta des genomes des diff especes 
		(remarque : et ne contenant pas de f)
	output :
		- pour chaque genome, fichier texte contenant la taille des contigs le composant
	"""
	fasta_files = glob.glob(species_fasta_rep + "/*f*")
	for file in fasta_files :
		#print file
		filout = open(file.split('/')[-1] + "_length", "w")
		fasta = open(file, "r")
		contig_len = 0
		for line in fasta :
			if ">" in line :
				contig_name = line.split()[0][1:]
				if contig_len != 0 : # pour ne pas mettre le compteur initialisé
					#print contig_len
					filout.write(str(contig_len) + "\n")
				#print contig_name,
				filout.write(str(contig_name) + "\t") 
				contig_len = 0
			if ">" not in line :
				flag = 1 #presence de seq
				contig_len = contig_len + len(line[:-1])
		#print contig_len
		filout.write(str(contig_len) + "\n")
		filout.close()
		fasta.close()

def dico_len_contigs(contig_len_rep) :
	"""
	Fonction permettant de renvoyer un dico contenant pour chaque genome les contigs et leur taille
	input :
		- contig_len_rep : repertoire contenant les fichiers (contigs - taille)
	output :
		- retour d'un dico avec cle = genome et valeur = dico contig avec cle = contig et valeur = taille
		{genome : {contig : taille}}
	"""
	len_contig = {}
	files = glob.glob(contig_len_rep + "/*")
	for file in files :
		file_name = file.split('_length')[0].split('/')[-1]
		#print file_name.split(".")[0]
		#print 
		len_contig[file_name.split(".")[0]] = {}
		filin = open(file, "r")
		#pour chaque ligne du fichier
		for line in filin :
			#ajouter ds dico[fichier][contig] = taille
			#print line.split()[0], line.split()[1]
			len_contig[file_name.split(".")[0]][line.split()[0]] = int(line.split()[1])
		filin.close()
	#verif
	#print len_contig
	return len_contig

def opr_like_regions(opr_type_rep) :
	"""
	Fonction permettant d'afficher le nombre de regions opr_like dans les fichiers resultats
	input :
		- opr_type_rep : repertoire contenant les fichiers *_opr_regions
	output :
		- affiche les resultats voulus
	"""
	files = glob.glob(opr_type_rep + "/*")
	for file in files :
		print file
		filin = open(file, 'r')
		nb_region = 0
		nb_1opr_by_region = 0
		for l in filin :
			nb_region += 1
			#print l
			#print l.split('[')[1][:-2]
			#print l.split('[')[1][:-2].count(",")
			#si nombre de prot OPR == 1 : (correspond a des regions ou y a 1 seule prot OPR)
			if l.split('[')[1][:-2].count(',') == 0 :
				#print l.split('[')[1][:-2]
				nb_1opr_by_region += 1
		filin.close()
		print nb_region, "(", nb_1opr_by_region, ")"

def replace_col1_name() :
	"""
	Fonction permettant de remplacer la premiere colonne du fichier (nom du query) pr les profils NCL
	"""
	proteomes = glob.glob("./*")
	for proteome in proteomes :
		files = glob.glob(proteome + "/*psiblast.csv")
		for file in files :
			col_name = file.split('vs')[0].split('/')[-1]
			os.system("sed -i 's/^1/" + col_name + "/' " + file)

def prot_result_Evalue_selection(proteome_rep, th_evalue) :
	"""
	Fonction permettant de sélectionner les résultats inf à un certain seuil de Evalue.
	input :
		- proteome_rep : repertoire contenant l'ens des rep faa
		- th_evalue : seuil de Evalue 
	output : 
		- un fichier par fichier de resultat blast ou psiblast contenant uniquement les hits selectionnés.
	"""
	proteomes = glob.glob(proteome_rep + "/*.faa")
	for proteome in proteomes :
		#print proteome
		files = glob.glob(proteome + "/*.csv")
		for file in files : 
			new_file = proteome + "/bis_" + file.split('/')[-1]
			#print new_file
			filin = open(file, 'r')
			filout = open(new_file, 'w')
			for l in filin :
				start = l.split('\t')[6]
				end = l.split('\t')[7]
				Evalue = float(l.split('\t')[16])
				#print Evalue
				if Evalue <= float(th_evalue) :
					#print '*', Evalue, start, end
					filout.write(l)
			filin.close()
			filout.close()

def region_analysis(results_rep, species_rep, OPR_length_file, contig_len_rep) :
	"""
	Fonction permettant de donner la composition de chaque region en domaine, gene ou cluster
	input :
		- les input de la fonction OPR_regions() car besoin du dico
		- results_rep : repertoire contentant les repertoires des diff resultats tblastn des genomes
	output :
		- affichage dans le terminal
	"""
	#dico de la taille des oprs
	OPR_qlen = dico_len_OPR(OPR_length_file)
	OPR_qlen["NCL"] = 1500
	OPR_qlen["NCL1"] = 1500
	OPR_qlen["NCL2"] = 1500
	#dico des infos sur les diff regions retrouvées dans le ome
	regions_dico = OPR_regions(species_rep, OPR_length_file, contig_len_rep)
	#for opr in regions_dico :
	#	print opr
	#	for region in regions_dico[opr] :
	#		print region, regions_dico[opr][region]
	#ouverture du fichier contenant toutes les regions avec tous les oprs
	file_path = results_rep + '/regions/' + species_rep.split('/')[-1] + '_opr_regions'
	filin = open(file_path, 'r')
	#ouverture d'un fichier resume par espece
	filename = results_rep + "/regions/resume/resume_" + species_rep.split('/')[-1] + ".csv"
	filout = open(filename, 'w')
	#pour G et T
	filout.write("Contig\tsstart-send\tRegion_length\tType_region\tOPR\t%OPR\t[(qstart-qend-strand)]\n")
	#pour P 
	#filout.write("Contig\tsstart-send\tRegion_length\tType_region\tOPR\t%OPR\t[(qstart-qend)]\n")
	#pour chaque ligne du fichier
	for line in filin :
		#print line
		#infos sur contig
		contig = line.split('\t')[0]
		##print contig
		cstart = line.split('\t')[2]
		cend = line.split('\t')[3]
		#opr_list pour l'analyse des genomes & transcriptomes
		opr_list = line.split('\t')[4][1:-2]
		#opr_list pour l'analyse des proteines
		#opr_list = line.split('\t')[4][1:-1]
		#print opr_list
		#nb d'opr dans la liste
		nb_opr_in_list = opr_list.count(',') + 1
		# affichage des infos du contig
		taille_contig = int(cend) - int(cstart)
		taille_contig_kb = taille_contig / 1000.0
		#type de contig = domaine, gene ou cluster de gene
		if taille_contig_kb < 5 :
			type_contig = "domain"
		elif 5 < taille_contig_kb < 50 :
			type_contig = "gene"
		else : 
			type_contig = "gene or cluster"
		#print "\t>" + contig + " : " + cstart + " - " + cend + " (" + str(taille_contig_kb) + " kb) ==> " + type_contig
		#pour chaque opr present dans ce contig en question
		for i in range(nb_opr_in_list) :
			#nom de l'opr
			opr = opr_list.split(',')[i].replace(" ","")[1:-1]
			#print opr
			#liste des regions retrouvées dans ce contig avec cet opr
			region_list = regions_dico[opr][contig]
			#initialisation d'une liste contenant l'ensemble des regions_opr retrouvées dans les hits et une variable = somme des % des oprs
			total_taille_region = 0
			list_opr_region = []
			#pour chacune de ces regions
			for region in region_list :
				##print region
				#print region[0], region[1]
				#taille_region_contig = 100 * abs(int(region[3])-int(region[2])) / float(taille_contig)
				# attention, taille de la proteine en aa et region donnée en nt 
				#taille_region_opr = 100 * abs(int(region[1])-int(region[0])) / (float((OPR_qlen[opr])))
				#total_taille_region = total_taille_region + taille_region_opr
				list_opr_region.append((region[0],region[1],region[-1]))
				#print opr + " : ",
				#affichage infos contig
				#print str(taille_region_contig) + " %" + " contig [" + str(region[2]) + " - " + str(region[3]) + "] <=> ",
				#affichage infos opr par region
				#print str(taille_region_opr) + " %" + " opr [" + str(region[0]) + " - " + str(region[1]) + "]"
			#affichage resumé des infos de l'ens des regions concernant l'opr en question
			#print opr + " : " + str(total_taille_region) + " % " + str(sorted(list_opr_region, key=lambda data: data[0]))
			#print  opr + " : " + str(sorted(list_opr_region, key=lambda data: data[0]))
			#print opr, OPR_frame
			#print list_opr_region
			real_regions = resume_regions(list_opr_region)
			#print real_regions
			#pour chaque region non chevauchantes
			for R in sorted(real_regions, key=lambda data: data[0]) :
				#calcul de la taille en pourcentage d'opr au niveau de cette region
				taille_region_opr = 100 * abs(int(R[1])-int(R[0])) / float((OPR_qlen[opr])) #R[0] = qstart, R[1] = qend
				#%opr present dans le contig en question 
				total_taille_region = total_taille_region + taille_region_opr
				
			#print contig + " : " + cstart + "-" + cend + " (" + str(taille_contig_kb) + " kb) " + type_contig,
			#print opr + " : " + str(round(total_taille_region,2)) + "% " + str(sorted(real_regions, key=lambda data: data[0]))
			#ecriture dans le fichier des infos contigs pr G and T
			#print contig + " : " + cstart + "-" + cend + "\t" + str(taille_contig_kb) + " kb\t" + type_contig + "\t",
			filout.write(contig + " : " + cstart + "-" + cend + "\t" + str(taille_contig_kb) + " kb\t" + type_contig + "\t")
			#ecriture dans le fichier des infos contigs pr P
			#filout.write(contig + " : " + cstart + "-" + cend + "\t" + str(taille_contig) + " res\t" + type_contig + "\t")
			#ecriture dans le fichier des infos opr /!\ modifier le real_regions pour T & P
			#print opr + " : " + str(round(total_taille_region,2)) + "%\t" + str(sorted(real_regions, key=lambda data: data[0])) + "\n"
			filout.write(opr + " : " + str(round(total_taille_region,2)) + "%\t" + str(sorted(real_regions, key=lambda data: data[0])) + "\n")
	filout.close()
	filin.close()

def resume_regions(list_opr_region) :
	"""
	Fonction permettant de retourner des coord de regions non chevauchantes
	input :
		- list_opr_region : liste de regions sous forme de tuples(chevauchantes ou pas)
	output :
		- retourne une liste contenant des regions non chevauchantes
	"""
	rmin = 'X'
	rmax = 'X'
	rstrand = 'X'
	#les gdes regions dites opr
	real_regions = {}
	#initialisation nb de region
	nb_region = 1
	real_regions[nb_region] = (rmin, rmax, rstrand)
	#m = start, M = end, S = strand
	for (m,M,frame) in list_opr_region :
		#min et max
		flag = 0
		if int(frame) < 0 :
			strand = "-"
		else :
			strand = "+"
		#pour les regions dejà recensées
		for n in real_regions :
			rmin = real_regions[n][0]
			rmax = real_regions[n][1]
			rstrand = real_regions[n][2]
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
			if (flag == 0 and M < rmin) or (flag == 0 and m > rmax) or (flag == 0 and strand != rstrand) :
				rmin = m
				rmax = M
				if strand != rstrand :
					rstrand = strand
			#print nb_region, rmin, rmax
			#si region chevauchante
			if flag == 1 :
				real_regions[n] = (rmin, rmax, rstrand)
			#print nb_region, m, M, rmin, rmax
		#sinon nouvelle region non chevauchante
		if flag == 0 :
				nb_region += 1
				#print nb_region, m, M, rmin, rmax
				real_regions[nb_region] = (rmin, rmax, rstrand)
	#liste de tuples de coord de regions
	list_R = []
	for i in real_regions :
		list_R.append(real_regions[i]) 
	return list_R

def Evalue_choice_proteome(rep_proteome) :
	"""
	Fonction permettant de plotter la distribution des Evalue des fichiers resultats CSV des proteomes
	input :
		- rep_proteome : repertoire contenant les rep des resultats des diff especes
	output :
		- plot de la distribution de l'ensemble des Evalues
	"""
	Evalue = []
	files = glob.glob(rep_proteome + "/*csv")
	for file in files : 
		filin = open(file, "r")
		for l in filin :
			evalue = float(l.split("\t")[16])
			#print evalue
			if evalue == 0.0 :
				evalue = 1e-200
			Evalue.append(-log(evalue))
		filin.close()
	"""
	species_rep = glob.glob(rep_proteome + "/*.faa")
	for species in species_rep :
		print species
			files = glob.glob(species + "/*csv")
			for file in files : 
				filin = open(file, "r")
				for l in filin :
					evalue = float(l.split("\t")[16])
					print evalue
					if evalue == 0.0 :
						evalue = 1e-200
					Evalue.append(-log(evalue))
				filin.close()
	"""
	plt.hist(np.array(Evalue), bins = 100)
	plt.title("Proteome Evalue Distribution")
	plt.xlabel("-log(Evalue)")
	plt.ylabel("Number of Evalue")
	plt.show()

def add_contig_len_column(resume_rep,contig_len_rep) :
	"""
	Fonction permettant d'ajouter au fichier resume une colonne contenant la taille de chaque contig
	"""
	#dico contenant la taille de chaque contig pour chaque esp
	dico_len = dico_len_contigs(contig_len_rep)
	#verif contenu dico
	#for i in dico_len :
	#	print i
	#ens des fichiers resume
	files = glob.glob(resume_rep + "/resume*")
	for file in files :
		#print file 
		#nom de l'espece pour pouvoir acceder directement au contenu du dico
		species = file.split('resume_')[1][:-4]
		filin = open(file, 'r')
		line = filin.readline()
		filout = open("resume1_" + species + ".csv", "w")
		filout.write("Contig\tContig_length\tsstart-send\tRegion_length\tType_region\tOPR\t%OPR\t[(qstart-qend-strand)]\n")
		for line in filin :
			#print line
			contig = line.split(":")[0]
			#pour genome et transcriptome
			last_part = line[len(contig) + 1:]
			#pour genome et proteome
			if "|" in contig :
				contig = contig.split("|")[-2]
			else :
				contig = contig[:-1]
			#pour transcriptome
			#contig = contig[:-1]
			#taille contig & proteine
			contig_len = dico_len[species][contig]
			#print contig + "\t" + str(contig_len) + "\t" + last_part
			#pour genome et transcriptome
			filout.write(contig + "\t" + str(contig_len) + "\t" + last_part)
			#pour proteome
			#part1 = last_part.split("\t")[0] + "\t" + last_part.split("\t")[1] + "\t"
			#type_region = last_part.split("\t")[2]
			#part2 = last_part[len(part1) + len(type_region) :]
			#type_region = int(last_part.split("\t")[0].split("-")[1]) - int(last_part.split("\t")[0].split("-")[0])
			#print type_region, contig_len
			#print type_region / float(contig_len)
			#print contig + "\t" + str(contig_len) + "\t" + part1 + str(round(100 * type_region / float(contig_len), 2)) + "%" + part2
			#filout.write(contig + "\t" + str(contig_len) + "\t" + part1 + str(round(100 * type_region / float(contig_len), 2)) + "%" + part2)
		filout.close()
		filin.close()


#Programme principal
if __name__ == "__main__":
	print 
	#Evalue_choice_proteome(sys.argv[1])
	#frame_fichiers_non_vide(sys.argv[1])
	#analysis_all_files(sys.argv[1], sys.argv[2])
	#liste = glob.glob(sys.argv[1] + '/*')
	#for i in liste :
	#	print i
	#	genome_quality(i)
	#genome_quality(sys.argv[1])
	#dico_len_OPR(sys.argv[1])
	#OPR_regions(sys.argv[1], sys.argv[2])
	#contigs_length(sys.argv[1])
	#dico_len_contigs(sys.argv[1])
	#********************Analyses genomes & transcriptomes
	#rangement des fichiers csv (etre dans le rep en question)
	#results_stockage(sys.argv[1])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py /data/mortaza/Lucifer/oprvsgenome/results/list/genome_list
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py /data/mortaza/Lucifer/oprvsgenome/results/list/transcriptome_list
	#enlever les fichiers vides
	#frame_fichiers_non_vide(sys.argv[1])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py .
	#all_analysis avec frame_fichiers_non_vide à lancer dans results
	#analysis_all_files(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py list/genome_list genome
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py list/transcriptome_list transcriptome
	#all_analysis avec infos sur les contigs : [contig, region, OPR]
	#analysis_all_files(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py list/genome_list genome
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py list/transcriptome_list transcriptome
	#compter le nombre de regions opr-like (nb de lignes) a partir du rep regular_regions. se placer ds le rep contenant les fichiers _opr_regions
	#opr_like_regions(sys.argv[1])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py .
	#*********************Analyses proteomes
	#rangement des fichiers csv (etre dans le rep en question)
	#results_stockage(sys.argv[1])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py /data/mortaza/Lucifer/oprvsgenome/results/list/proteome_list
	#enlever les fichiers vides ==> pas de fichiers vides
	#frame_fichiers_non_vide(sys.argv[1])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py .
	#remplacer premiere colonne des fichiers NCL par le nom de NCL
	#sed -i 's/^1/NCL/' nom_du_fichier
	#replace_col1_name()
	#analyse contigs des proteomes
	#analysis_all_files(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py list/proteome_list proteome
	#compter le nombre de regions opr-like (nb de lignes) a partir du rep regular_regions. se placer ds le rep contenant les fichiers _opr_regions
	#opr_like_regions(sys.argv[1])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py .
	#selection des meilleurs hits
	#prot_result_Evalue_selection(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py . 1e-06
	#analyse des regions ds /data/mortaza/Lucifer/oprvsgenome/results
	#analysis_all_files(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py list/genome_list genome
	#calculer la taille de chaque proteine des proteomes 
	#contigs_length(sys.argv[1])
	#ajout colonne contig (y compris les proteines) length (a lancer dans resume)
	#add_contig_len_column(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py . /data/mortaza/Data/genomes/genome_len


#Notes & Commentaires
"""
- /home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py
- makeblastdb des OPRs.fasta
mortaza@athena:~$ scp -r mortaza@172.26.16.29:/home/mortaza/Documents/16S_PhylogenyTree/New_genomes .
- test commande tblastn ds athena
tblastn -query Chlre_OPR9.fasta -db /home/mortaza/New_genomes/PGGS01.1.fsa_nt -outfmt \
"6 qseqid qseq qstart qend sseqid sseq sstart send sstrand sframe length pident nident \
ppos gaps gapopen evalue bitscore qcovs qcovhsp qcovus" -evalue 10e-09 -out output_test.txt -num_threads 96
- repertoires 
/workdir/ibpc_team/lafontaine_team/mortaza/ingrid2shogofa/oprvsgenomes/fasta_files/Regular1
/workdir/ibpc_team/lafontaine_team/mortaza/ingrid2shogofa/oprvsgenomes/fasta_files/Regular2
/workdir/ibpc_team/lafontaine_team/mortaza/ingrid2shogofa/oprvsgenomes/fasta_files/Regular3
- distribution de la taille des genomes (/data/mortaza/Data/genomes)
- type des OPRs (/data/mortaza/Lucifer/oprvsopr/genes)
- modif d'une chaine de car par une autre
sed -ie 's/Chlre\_OPR67/Chlre\_OPR67\=Cre14\.g633150\.t1\.3/g' OPR_qlen.txt
- dps rep /data/.../genome, lancer
/home/mortaza/Documents/src/OPRvsgenome/oprvsgenome.py ../genome_list.txt .
- construction d'un multifasta de NCL ds rep results
cat /data/mortaza/Data/OPR/opr_byclass_chlamy/NCL/*.fasta >> NCL.fasta

"""