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
from result_analysis_part2 import *
from result_T2_analysis import *
from result_G2_analysis import *
sys.path.append('/data/mortaza/src/OPRvsOPR_rep/')
from cluster_quality import *

#Fonctions

def return_regions_non_chevauchantes(lhc_homologous_csv_file, type_ome, th) :
	"""
	Fonction permettant de retourner des regions non chevauchantes. 
	input :
		- lhc_homologous_csv_file : fichier csv résumé des résultats de blast
		obtenu avec la fonction dico_blast_results_pt()
		- type_ome : P ou T ou G
		- th : seuil pr selectionner les proteines et les transcrits
	output : 
		- retourne un dico de regions non chevauchantes
	"""
	#initialisation dico
	dico =  {}
	#lecture du fichier
	filin = open(lhc_homologous_csv_file, "r")
	#pr chaque ligne du fichier
	for line in filin :
		#recolte des infos nec
		ome = line.split(" ")[0]
		contig = line.split(" ")[1]
		sstart = int(line.split(" ")[2])
		send = int(line.split(" ")[3])
		strand = line.split(" ")[4]
		query = line.split(" ")[-3]
		qstart = int(line.split(" ")[-2])
		qend =  int(line.split(" ")[-1][:-1])
		#enlever les resultats de cette proteine
		if query != "Cre10.g454734" :
			#initialisation des elmts du dico
			if ome not in dico :
				dico[ome] = {}
			if contig not in dico[ome] :
				dico[ome][contig] = []
			#enregistrement des infos ds le dico
			dico[ome][contig].append((sstart, send, strand, qstart, qend))
	filin.close()
	#pr chaque esp
	for sp in dico :
		#pr chaque contig 
		for ctg in dico[sp] :
			#print "***"
			#liste des regions
			old_list = sorted(dico[sp][ctg], key=lambda data: data[0])
			#print old_list
			new_list = resume_regions_with_query(old_list)
			#print new_list
			#for i in new_list :
				#print sp, ctg, i
			#pr les PT
			if type_ome == "T" or type_ome == "P" :
				#initialisation des variables de positions
				true_sstart = true_send = "X"
				true_qstart = true_qend = "X"
				true_strand = "X"
				#mise en commun de regions
				#ss : sstart, se : send, st : stand, qs : qstart, qe : qend
				#pr chaque region
				for (ss, se, st, qs, qe) in new_list : 
					#print ss, se, st, qs, qe
					if true_sstart == "X" :
						true_sstart = ss
						true_send = se
						true_strand = st
						true_qstart = qs
						true_qend = qe
					else :
						if ss < true_sstart :
							true_sstart = ss
						if se > true_send :
							true_send = se
						if qs < true_qstart :
							true_qstart = qs
						if qe > true_qend :
							true_qend = qe
			if type_ome == "G" :
				G_new_list = genomic_prot_hits(new_list)
				for infos in G_new_list :
					#print infos
					true_sstart = infos[0]
					true_send = infos[1]
					true_strand = infos[2]
					true_qstart = infos[3]
					true_qend = infos[4]
			#couverture du query
			query_coverage = (true_qend - true_qstart)
			#limite pour selectionner les proteines ou transcrits
			#correspondant à une proteine LHCa
			th_coverage = 250 * int(th) / 100
			#si nb de res de la proteine sup ou egale à la limite, alors ok
			if query_coverage >= th_coverage :
				print sp, ctg, true_sstart, true_send, true_strand, query_coverage, th_coverage


def genomic_prot_hits(liste_regions) :
	"""
	Fonction permettant de donner les regions protéiques LHCa d'un contig
	input :
		- liste_regions : liste des regions retrouvés dans le contig génomique
	output :
		- retourne une liste de regions proteiques identifiées dans le contig
	"""
	#nb de prot identifiées dans le contig
	nb_prot_p = 0
	nb_prot_m = 0
	#initialisation des valeurs 
	r_sstart = r_send = r_qstart = r_qend = "X"
	#initialisation dico contenant les differentes infos des proteines
	#p pour brin plus et m pour brin moins
	dico_prot_p = {}
	dico_prot_m = {} 
	dico_prot_p[nb_prot_p] = (r_sstart, r_send, "+", r_qstart, r_qend)
	dico_prot_m[nb_prot_m] = (r_sstart, r_send, "-", r_qstart, r_qend)
	#pr chaque region non chevauchante
	for (sstart, send, strand, qstart, qend) in liste_regions :
		sstart = int(sstart)
		send = int(send)
		qstart = int(qstart)
		qend = int(qend)
		flag = 0
		if strand == "+" :
			r_sstart = dico_prot_p[nb_prot_p][0]
			r_send = dico_prot_p[nb_prot_p][1]
			r_qstart = dico_prot_p[nb_prot_p][3]
			r_qend = dico_prot_p[nb_prot_p][4]
			#initialisation avec les vraies valeurs
			if r_sstart == "X" :
				r_sstart = sstart
				r_send = send
				r_qstart = qstart
				r_qend = qend
				flag = 1
			#sinon tester les differents cas
			if r_qstart <= qstart :
				r_send = send
				r_qend = qend
				flag = 1
			if flag == 1 :
				dico_prot_p[nb_prot_p] = (r_sstart, r_send, "+", r_qstart, r_qend)
			if flag == 0 :
				nb_prot_p += 1
				dico_prot_p[nb_prot_p] = (sstart, send, "+", qstart, qend)
		if strand == "-" :
			r_sstart = dico_prot_m[nb_prot_m][0]
			r_send = dico_prot_m[nb_prot_m][1]
			r_qstart = dico_prot_m[nb_prot_m][3]
			r_qend = dico_prot_m[nb_prot_m][4]
			#initialisation avec les vraies valeurs
			if r_sstart == "X" :
				r_sstart = sstart
				r_send = send
				r_qstart = qstart
				r_qend = qend
				flag = 1
			#sinon tester les differents cas
			if r_qstart >= qstart :
				r_qstart = qstart
				r_send = send
				if r_qend < qend :
					r_qend = qend
				flag = 1
			if flag == 1 :
				dico_prot_m[nb_prot_m] = (r_sstart, r_send, "-", r_qstart, r_qend)
			if flag == 0 :
				nb_prot_m += 1
				dico_prot_m[nb_prot_m] = (sstart, send, "-", qstart, qend)
	#lecture des deux dicos non vides et stockage dans une liste :
	LHCa_proteins_list = [] 
	if dico_prot_p[0][0] != "X" :
		for i in dico_prot_p :
			#print dico_prot_p[i]
			LHCa_proteins_list.append(dico_prot_p[i])
	if dico_prot_m[0][0] != "X" :
		for i in dico_prot_m :
			#print dico_prot_m[i]
			LHCa_proteins_list.append(dico_prot_m[i])
	return LHCa_proteins_list

			
def moy_prot(fasta_len) :
	"""
	Fonction calculant la taille moyenne des protéines.
	input :
		- fasta_len : fichier de deux colonnes (prot - taille)
	output :
		- affiche la taille moyenne des protéines arrondie à l'unité
	"""
	#lecture du fichier
	filin = open(fasta_len, "r")
	#initialisation de la liste contenant les differentes tailles
	L_taille = []
	#pour chaque ligne du fichier
	for line in filin :
		#recolte des infos
		protein = line.split("\t")[0]
		prot_len = int(line.split("\t")[1][:-1])
		#ne pas prendre en compte la taille de cette proteine
		if protein != "Cre10_g454734" :
			L_taille.append(prot_len)
	filin.close()
	#affiche la moyenne arrondie à l'unité
	print round(sum(L_taille)/float(len(L_taille)),0)


def resume_regions_with_query(list_opr_region) :
	"""
	Fonction permettant de retourner des coord de regions non chevauchantes en prenant 
	en compte les start et end des query
	input :
		- list_opr_region : liste de regions sous forme de tuples(chevauchantes ou pas)
	output :
		- retourne une liste contenant des regions non chevauchantes
	"""
	rmin = 'X'
	rmax = 'X'
	rstrand = 'X'
	qmin = 'X'
	qmax = 'X'
	#les gdes regions dites opr
	real_regions = {}
	#initialisation nb de region
	nb_region = 1
	real_regions[nb_region] = (rmin, rmax, rstrand, qmin, qmax)
	#m = start, M = end, S = strand
	for (m,M,strand,qm,qM) in list_opr_region :
		#min et max
		m = int(m); M = int(M); qm = int(qm); qM = int(qM)
		flag = 0
		#pour les regions dejà recensées
		for n in real_regions :
			rmin = real_regions[n][0]
			rmax = real_regions[n][1]
			rstrand = real_regions[n][2]
			qmin = real_regions[n][3]
			qmax = real_regions[n][4]
			#initialisation du rmin et rmax et rstrand avec les vraies valeurs
			if rmin == 'X' and rmax == 'X' and rstrand == 'X' :
				rmin = m
				rmax = M
				rstrand = strand
				qmin = qm
				qmax = qM
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
				if qm < qmin : 
					qmin = qm
				if qM > qmax :
					qmax = qM
				real_regions[n] = (rmin, rmax, rstrand, qmin, qmax)
			#print nb_region, m, M, rmin, rmax
		#sinon nouvelle region non chevauchante
		if flag == 0 :
				nb_region += 1
				qmin = qm
				qmax = qM
				#print nb_region, m, M, rmin, rmax
				real_regions[nb_region] = (rmin, rmax, rstrand, qmin, qmax)
	#liste de tuples de coord de regions
	list_R = []
	for i in real_regions :
		list_R.append(real_regions[i]) 
	return list_R


def fasta_creation(ome, ctg, sstart, send, rep_fasta_files) :
	"""
	Fonction permettant de creer un fichier fasta
	input : 
		- ctg
		- sstart
		- send
		- rep_fasta_files : repertoire contenant les fichiers fasta des omes
	output :
		- creation de fichier fasta pr chaque contig voulu
	"""
	#dico contenant les infos des fichiers fasta d'un ome
	fastas = dico_info_fasta(rep_fasta_files)
	print ">" + ctg + "_" + str(sstart) + "_" + str(send)
	output = ctg + "_" + str(sstart) + "_" + str(send) + ".fasta"
	filout = open(output, "w")
	filout.write(">" + ctg + "_" + str(sstart) + "_" + str(send) + "\n")
	filout.write(fastas[ome][ctg][int(sstart):int(send)] + "\n")
	filout.close()


def prot_fasta_creation_from_data(file, rep_fasta_files) :
	"""
	Fonction permettant de créer les fichiers fastas voulus
	input :
		- file
		- rep_fasta_files
	output :
		- creation de fichiers fasta
	"""
	filin = open(file, "r")
	for line in filin :
		ome = line.split(" ")[0]
		ome_rep = ome
		ctg = line.split(" ")[1]
		sstart = line.split(" ")[2]
		send = line.split(" ")[3]
		#print ome, ctg, sstart, send
		if not os.path.exists(ome_rep):
			os.makedirs(ome_rep)
		else :
			os.chdir(ome_rep)
			fasta_creation(ome, ctg, sstart, send, rep_fasta_files)
			os.chdir("..")
	filin.close()


def exonerate_file_reader(exonerate_file_rep) :
	"""
	Fonction permettant de lire des fichiers obtenus avec exonerate 
	input :
		- exonerate_file_rep : rep des fichiers obtenus avec la commande 
		exonerate - modele protein2genome - et --ryo ">%ti %td %pi\n%tcs"
		et --showalignment FALSE
	output :
		- affiche les infos voulues
	"""
	#creation dico contenant les infos du fichier
	dico = {}
	#ens des fichiers resultats de exonerate
	files = glob.glob(exonerate_file_rep + "/*")
	#pr chaque fichier
	for file in files :
		ome = file.split("/")[-1].split(".fasta")[0]
		if ome not in dico :
			dico[ome] = {}
		#lecture du fichier
		filin = open(file, "r")
		#ne pas lire les headers
		line = filin.readline()
		line = filin.readline()
		line = filin.readline()
		seq = ""
		flag = 0
		#pr chaque ligne du fichier
		for line in filin :
			#si fin de la ligne ou nouvelle seq
			if ("vulgar" in line and seq != "") or ("exonerate" in line) :
				flag +=1
			#affichage des resultats complets
			if flag == 3:
				#print query, qstart, qend, p_identity
				#print subject, sstart, send, strand
				#print seq
				flag = 0
				#enregistrement des resultats
				if subject not in dico[ome] :
					dico[ome][subject] = []
				dico[ome][subject].append((sstart, send, strand, query, qstart, qend, p_identity))
				#print subject, query, p_identity
				#print (sstart, send, strand, query, qstart, qend, p_identity)
			#ligne infos vulgar
			if "vulgar" in line :
				query = line.split(" ")[1] 
				qstart = line.split(" ")[2]
				qend = line.split(" ")[3]
				subject =  line.split(" ")[5]
				sstart = line.split(" ")[6]
				send = line.split(" ")[7]
				strand = line.split(" ")[8]
				flag = 1
				#print subject, query,
				#print query, qstart, qend
				#print subject, sstart, send, strand
			#ligne header de la seq contenant le pourcentage d'identité
			if ">" in line :
				p_identity = line.split(" ")[2][:-1]
				flag += 1
				#print p_identity
				#initialisation de la variable contenant la séq 
				seq = ""
			#ligne de la seq
			else :
				seq += line[:-1]
		filin.close()
	#lecture du dico
	#pr chaque esp 
	for sp in dico :
		#pr chaque subject
		for ctg in dico[sp] :
			#print ctg
			maxi_pi = 0
			maxi_sstart = "X"
			maxi_send = "X"
			maxi_strand = "X"
			maxi_query = "X"
			maxi_qstart = "X"
			maxi_qend = "X"
			#pr chaque query
			for i in dico[sp][ctg] :
				#print i
				#print i
				if float(i[-1]) > maxi_pi :
					maxi_pi = float(i[-1])
					maxi_sstart = int(i[0])
					maxi_send = int(i[1])
					maxi_strand = i[2]
					maxi_query = i[3]
					maxi_qstart = int(i[4])
					maxi_qend = int(i[5])
			print ctg, maxi_sstart, maxi_send, maxi_strand, maxi_query, maxi_qstart, maxi_qend, maxi_pi


def resume_regions_with_query_and_Evalue(list_opr_region) :
	"""
	? fct pas terminée ...
	Fonction permettant de retourner des coord de regions non chevauchantes en prenant 
	en compte les start et end des query
	input :
		- list_opr_region : liste de regions sous forme de tuples(chevauchantes ou pas)
	output :
		- retourne une liste contenant des regions non chevauchantes
	"""
	rmin = 'X'
	rmax = 'X'
	rstrand = 'X'
	qmin = 'X'
	qmax = 'X'
	rEvalue = 'X'
	#les gdes regions dites opr
	real_regions = {}
	#initialisation nb de region
	nb_region = 1
	real_regions[nb_region] = (rmin, rmax, rstrand, qmin, qmax, rEvalue)
	#m = start, M = end, S = strand
	for (m,M,strand,qm,qM,Ev) in list_opr_region :
		#min et max
		m = int(m); M = int(M); qm = int(qm); qM = int(qM)
		Ev = float(Ev)
		flag = 0
		#pour les regions dejà recensées
		for n in real_regions :
			rmin = real_regions[n][0]
			rmax = real_regions[n][1]
			rstrand = real_regions[n][2]
			qmin = real_regions[n][3]
			qmax = real_regions[n][4]
			rEvalue = real_regions[n][5]
			#initialisation du rmin et rmax et rstrand avec les vraies valeurs
			if rmin == 'X' and rmax == 'X' and rstrand == 'X' :
				rmin = m
				rmax = M
				rstrand = strand
				qmin = qm
				qmax = qM
				rEvalue = Ev
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
			else : 
				rmin = 1
			#print nb_region, rmin, rmax
			#si region chevauchante
			if flag == 1 :
				if qm < qmin : 
					qmin = qm
				if qM > qmax :
					qmax = qM
				if float(Ev) < float(rEvalue) :
					rEvalue = Ev
				real_regions[n] = (rmin, rmax, rstrand, qmin, qmax, rEvalue)
			#print nb_region, m, M, rmin, rmax
		#sinon nouvelle region non chevauchante
		if flag == 0 :
				nb_region += 1
				qmin = qm
				qmax = qM
				#print nb_region, m, M, rmin, rmax
				real_regions[nb_region] = (rmin, rmax, rstrand, qmin, qmax, rEvalue)
	#liste de tuples de coord de regions
	list_R = []
	for i in real_regions :
		list_R.append(real_regions[i]) 
	return list_R


def genomic_prot_hits_Evalue(liste_regions) :
	"""
	Fonction permettant de donner les regions protéiques LHCa d'un contig
	input :
		- liste_regions : liste des regions retrouvés dans le contig génomique
	output :
		- retourne une liste de regions proteiques identifiées dans le contig
	"""
	#nb de prot identifiées dans le contig
	nb_prot_p = 0
	nb_prot_m = 0
	#initialisation des valeurs 
	r_sstart = r_send = r_qstart = r_qend = r_evalue = "X"
	#initialisation dico contenant les differentes infos des proteines
	#p pour brin plus et m pour brin moins
	dico_prot_p = {}
	dico_prot_m = {} 
	dico_prot_p[nb_prot_p] = (r_sstart, r_send, "+", r_qstart, r_qend, r_evalue)
	dico_prot_m[nb_prot_m] = (r_sstart, r_send, "-", r_qstart, r_qend, r_evalue)
	#pr chaque region non chevauchante
	for (sstart, send, strand, qstart, qend, Evalue) in liste_regions :
		sstart = int(sstart)
		send = int(send)
		qstart = int(qstart)
		qend = int(qend)
		Evalue = float(Evalue)
		flag = 0
		if strand == "+" :
			r_sstart = dico_prot_p[nb_prot_p][0]
			r_send = dico_prot_p[nb_prot_p][1]
			r_qstart = dico_prot_p[nb_prot_p][3]
			r_qend = dico_prot_p[nb_prot_p][4]
			r_evalue = dico_prot_p[nb_prot_p][5]
			#initialisation avec les vraies valeurs
			if r_sstart == "X" :
				r_sstart = sstart
				r_send = send
				r_qstart = qstart
				r_qend = qend
				r_evalue = Evalue
				flag = 1
			#sinon tester les differents cas
			if r_qstart <= qstart :
				r_send = send
				r_qend = qend
				flag = 1
			if flag == 1 :
				if Evalue < r_evalue :
					r_evalue = Evalue
				dico_prot_p[nb_prot_p] = (r_sstart, r_send, "+", r_qstart, r_qend, r_evalue)
			if flag == 0 :
				nb_prot_p += 1
				dico_prot_p[nb_prot_p] = (sstart, send, "+", qstart, qend, Evalue)
		if strand == "-" :
			r_sstart = dico_prot_m[nb_prot_m][0]
			r_send = dico_prot_m[nb_prot_m][1]
			r_qstart = dico_prot_m[nb_prot_m][3]
			r_qend = dico_prot_m[nb_prot_m][4]
			r_evalue = dico_prot_m[nb_prot_m][5]
			#initialisation avec les vraies valeurs
			if r_sstart == "X" :
				r_sstart = sstart
				r_send = send
				r_qstart = qstart
				r_qend = qend
				r_evalue = Evalue
				flag = 1
			#sinon tester les differents cas
			if r_qstart >= qstart :
				r_qstart = qstart
				r_send = send
				if r_qend < qend :
					r_qend = qend
				flag = 1
			if flag == 1 :
				if Evalue < r_evalue :
					r_evalue = Evalue
				dico_prot_m[nb_prot_m] = (r_sstart, r_send, "-", r_qstart, r_qend, r_evalue)
			if flag == 0 :
				nb_prot_m += 1
				dico_prot_m[nb_prot_m] = (sstart, send, "-", qstart, qend, Evalue)
	#lecture des deux dicos non vides et stockage dans une liste :
	LHCa_proteins_list = [] 
	if dico_prot_p[0][0] != "X" :
		for i in dico_prot_p :
			#print dico_prot_p[i]
			LHCa_proteins_list.append(dico_prot_p[i])
	if dico_prot_m[0][0] != "X" :
		for i in dico_prot_m :
			#print dico_prot_m[i]
			LHCa_proteins_list.append(dico_prot_m[i])
	return LHCa_proteins_list


def return_regions_non_chevauchantes_Evalue(lhc_homologous_csv_file, type_ome, th, Evalue_th) :
	"""
	Fonction permettant de retourner des regions non chevauchantes. 
	input :
		- lhc_homologous_csv_file : fichier csv résumé des résultats de blast
		obtenu avec la fonction dico_blast_results_pt()
		- type_ome : P ou T ou G
		- th : seuil pr selectionner les proteines et les transcrits
		- Evalue_th : seuil de selection pour les Evalues
	output : 
		- retourne un dico de regions non chevauchantes
	"""
	#initialisation dico
	dico =  {}
	#lecture du fichier
	filin = open(lhc_homologous_csv_file, "r")
	#pr chaque ligne du fichier
	for line in filin :
		#recolte des infos nec
		ome = line.split(" ")[0]
		contig = line.split(" ")[1]
		sstart = int(line.split(" ")[2])
		send = int(line.split(" ")[3])
		strand = line.split(" ")[4]
		query = line.split(" ")[-3]
		Evalue = line.split(" ")[5]
		qstart = int(line.split(" ")[-2])
		qend =  int(line.split(" ")[-1][:-1])
		#enlever les resultats de cette proteine
		if query != "Cre10.g454734" :
			#initialisation des elmts du dico
			if ome not in dico :
				dico[ome] = {}
			if contig not in dico[ome] :
				dico[ome][contig] = []
			#enregistrement des infos ds le dico
			dico[ome][contig].append((sstart, send, strand, qstart, qend, Evalue))
	filin.close()
	#pr chaque esp
	for sp in dico :
		#pr chaque contig 
		for ctg in dico[sp] :
			#print "***"
			#liste des regions
			old_list = sorted(dico[sp][ctg], key=lambda data: data[0])
			#print old_list
			new_list = resume_regions_with_query_and_Evalue(old_list)
			#print new_list
			#for i in new_list :
				#print sp, ctg, i
			
			#pr les PT
			if type_ome == "T" or type_ome == "P" :
				#initialisation des variables de positions
				true_sstart = true_send = "X"
				true_qstart = true_qend = "X"
				true_strand = "X"
				#mise en commun de regions
				#ss : sstart, se : send, st : stand, qs : qstart, qe : qend
				#pr chaque region
				for (ss, se, st, qs, qe, ev) in new_list : 
					#print ss, se, st, qs, qe
					if true_sstart == "X" :
						true_sstart = ss
						true_send = se
						true_strand = st
						true_qstart = qs
						true_qend = qe
						true_evalue = ev
					else :
						if ss < true_sstart :
							true_sstart = ss
						if se > true_send :
							true_send = se
						if qs < true_qstart :
							true_qstart = qs
						if qe > true_qend :
							true_qend = qe
						if ev < true_evalue : 
							true_evalue = ev
			if type_ome == "G" :
				G_new_list = genomic_prot_hits_Evalue(new_list)
				for infos in G_new_list :
					#print infos
					true_sstart = infos[0]
					true_send = infos[1]
					true_strand = infos[2]
					true_qstart = infos[3]
					true_qend = infos[4]
					true_evalue = infos[5]
			#couverture du query
			query_coverage = (true_qend - true_qstart)
			#limite pour selectionner les proteines ou transcrits
			#correspondant à une proteine LHCa
			th_coverage = 250 * int(th) / 100
			#si nb de res de la proteine sup ou egale à la limite, alors ok
			#selection par la taille de la couverture d'une proteine LHCa
			if query_coverage >= th_coverage :
				#selection par la Evalue
				if float(true_evalue) <= float(Evalue_th) :
					print sp, ctg, true_sstart, true_send, true_strand, query_coverage, th_coverage, true_evalue
			

def Chre_LHCa_P_th(list_file, csv_files) :
	"""
	Fonction permettant de determiner le seuil pour avoir seulement les 
	protéines de LHCa (et non les LHCb).
	input :
		- list_file : 
		- csv_files : liste de l'ens des rep contenant les fichiers csv 
		et obtenue avec glob.glob() 
	"""
	#enregistrement des infos utiles du blastp dans un dico
	dico = dico_blast_results_pt(list_file, csv_files)
	#lecture du dico 
	#pr chaque esp
	for ome in dico :
		#si y a des infos (dico du ome non vide)
		if dico[ome] != {} :
			#pr chaque contig
			for contig in dico[ome] :
				#print ">", contig,
				mini_Evalue = 2000.0
				mini_query = "X"
				mini_seq = "X"
				#pr chaque region retrouvée
				for region in dico[ome][contig] :
					#print region
					if float(region[3]) < mini_Evalue :
						mini_Evalue = float(region[3])
						mini_query = region[4]
						mini_seq = region[-1]
				if mini_query in contig:
					print ">", contig,
					#print mini_seq
					print mini_Evalue, mini_query
				#affichage des seq pr avoir le fichier multifasta de prot
				#print mini_seq

#Programme principal
if __name__ == "__main__":
	#1°) Transformation du fichier multifasta en plusieurs fichiers fasta 
	#(bdd query)
	#2°) Suppression des fichiers vides issus de blast
	#frame_fichiers_non_vide(sys.argv[1])
	#3°) Stockage des resultats par espece
	#results_stockage(sys.argv[1])
	#4°) Avoir le fichier de resume des infos de blast
	csv_files_g = glob.glob("/data/mortaza/Lucifer/evolution/LHCa_evolution/genomes/*")
	csv_files_t = glob.glob("/data/mortaza/Lucifer/evolution/LHCa_evolution/transcriptomes/*")
	csv_files_p = glob.glob("/data/mortaza/Lucifer/evolution/LHCa_evolution/proteomes/*")
	#dico_blast_results_pt(sys.argv[1], csv_files_p)
	#5°) Avoir les regions non chevauchantes
	#return_regions_non_chevauchantes(sys.argv[1])
	#6°) Taille des proteines LHCa
	#seq_length(sys.argv[1])
	#7°) Calcul de la taille moyenne des proteines LHCa
	#moy_prot(sys.argv[1])
	#8°) Regions non chevauchantes en fct start et end du subject et du query
	#return_regions_non_chevauchantes(sys.argv[1], "G", 80)
	#9°) Proteines obtenues avec Chre
	#prot_fasta_creation_from_data(sys.argv[1], sys.argv[2])
	#10°) Lecture de fichiers de resultats exonerate
	#exonerate_file_reader(sys.argv[1])
	#11°) 
	#return_regions_non_chevauchantes_Evalue(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	#12°) Analyse de blastp de Chre pour avoir le fichier resumé
	csv_files_p = glob.glob("/data/mortaza/Lucifer/evolution/LHCa_evolution/analysis/Chre_th/P_Chre_th/blastp_result/*")
	Chre_LHCa_P_th(sys.argv[1], csv_files_p)

"""
Notes & Commentaires :
- /data/mortaza/src/LHCa/LHCa_evolution.py
- multifasta en fasta
/home/mortaza/Documents/src/OPRvsgenome/multifasta_to_fasta.bash OPRs.fasta
- rep des listes des omes
/data/mortaza/Lucifer/evolution/OPR_evolution/oprvsgenome/results/list
- ligne de commande exonerate
"""