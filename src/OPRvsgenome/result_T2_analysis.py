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

#Fonctions

def dico_hit_T(transcriptome2_csv_file) :
	"""
	Fonction permettant de retourner un dico contenant les infos nec du 
	fichier csv donné en entrée
	input : 
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output :
		- retourne un dico [ome][transcrit][sframe][(sstart, send, sseq)]
	"""
	#initialisation dico 
	dico = {}
	#lecture du fichier
	filin = open(transcriptome2_csv_file, "r")
	#pr chaque ligne du fichier
	for line in filin :
		#recup des infos imp
		sep = line.split(" ")
		ome = sep[0]
		transcrit = sep[1]
		sstart = int(sep[2])
		send = int(sep[3])
		sframe = sep[4]
		sseq = sep[-1][:-1]
		#si ome pas ds dico, l'initialiser (dico)
		if ome not in dico :
			dico[ome] = {}
		#si transcrit pas ds dico, l'initialiser (dico)
		if transcrit not in dico[ome] :
			dico[ome][transcrit] = {}
		#si frame pas ds dico, l'initialiser (liste)
		if sframe not in dico[ome][transcrit] :
			dico[ome][transcrit][sframe] = []
		#enregistrement des resultats ds le dico
		dico[ome][transcrit][sframe].append((sstart, send, sseq))
		dico[ome][transcrit][sframe] = list(set(dico[ome][transcrit][sframe]))
	filin.close()
	#verif du dico
	#print "ome\ttranscrit\tsframe\tsstart\tsend\tsseq"
	for i in dico :
		#print i
		for j in dico[i] :
			#print j
			for k in dico[i][j] :
				#print k
				for l in sorted(dico[i][j][k], key=lambda data: data[0]) :
					#print i, j, k, 
					for m in l :
						#print m,
						a = 1 #qd print enlevé
					#print
	return dico

def dico_1transcrit_1seq(transcriptome2_csv_file) :
	"""
	Fonction permettant de faire correspondre 1 seq à 1 frame. 
	input : 
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output :
		- retourne un dico où 1 seq <=> 1 frame.
	"""
	dico_to_return_1transcrit_1seq = {}
	#recup du dico de la fct dico_hit_T()
	dico_hit = dico_hit_T(transcriptome2_csv_file)
	#pr chaque esp
	for i in dico_hit :
		dico_to_return_1transcrit_1seq[i] = {}
		#pr chaque transcrit
		for j in dico_hit[i] :
			dico_to_return_1transcrit_1seq[i][j] = {}
			#print i, j
			#creation liste contenant l'ens des 1seq-1frame
			L_frames = []
			#pr chaque frame
			for k in dico_hit[i][j] :
				#print "==> frame", k, 
				#sorted par rapport au sstart
				all_frames = sorted(dico_hit[i][j][k], key=lambda data: data[0])
				#print ", before : ", len(all_frames), ", after : ", 
				#rassembler les regions chevauchantes
				##si sur brin "+"
				if int(k) > 0 : 
					new_frames = regions_chevauchantes_plus(all_frames)
					#print len(new_frames)
				##si sur brin "-"
				else :
					new_frames = regions_chevauchantes_moins(all_frames)
					#print len(new_frames)
				#print new_frames
				#for l in new_frames :	
				#	print l
					#for m in l :
					#	print m,
					#print
				#mettre ens les regions non chevauchantes pour construire le
				#transcrit
				region_frame = regions_transcrit_per_frame(new_frames)
				#si pas de stop dans la seq, l'enregistrer
				if "*" not in region_frame[-1] :
					#print k, region_frame[-1]
					L_frames.append((k, region_frame[-1]))
				#print len(region_frame[-1])
			#initialisation des variables pour la seq qui sera choisie
			f_chosen = "X"
			s_chosen = "X"
			#si infos, choisir la meilleure
			if len(L_frames) != 0 :
				for infos in L_frames :
					if len(infos[-1]) > len(s_chosen) :
						s_chosen = infos[-1]
						f_chosen = infos[0]
			if f_chosen == "X" and s_chosen == "X" :
				#print i, j
				dico_to_return_1transcrit_1seq[i][j] = "X"
			else : 
				#print i, j, f_chosen, s_chosen
				dico_to_return_1transcrit_1seq[i][j] = s_chosen
	return dico_to_return_1transcrit_1seq

def regions_chevauchantes_plus(list_regions) :
	"""
	Fonction permettant de mettre ens les regions qui se chevauchent.
	input : 
		- list_regions : liste des regions -pr m frame- [(sstart, send, sseq)]
	output :
		- retourne la meme liste sans regions chevauchantes.
	"""
	rmin = "X"; rmax = "X"; rseq = "X"
	dico = {}
	nb_region = 1
	dico[nb_region] = (rmin, rmax, rseq)
	for (m, M, s) in list_regions :
		m = int(m); M = int(M)
		for n in dico :
			rmin = dico[n][0]
			rmax = dico[n][1]
			rseq = dico[n][2]
			flag = 0
			#initialisation
			if rmin == "X" :
				rmin = m; rmax = M; rseq = s
				flag = 1
			#chevauchement des hits
			##cas 1 : la region se trouve dans la grande region
			elif rmin <= m <= rmax and rmin <= M <= rmax :
				flag = 1
			##cas 2 : la region touche un bout de la grande region
			elif rmin <= m <= rmax and M > rmax :
				#print rmin, rmax
				#print rseq
				#print m, M
				#print s
				chevauchement_codon = (rmax - m + 1) / 3
				#print chevauchement_codon
				#verif que la partie qui se chevauche correspond à la meme seq
				if rseq[-chevauchement_codon:] == s[:chevauchement_codon] :
					new_rseq = rseq[:-chevauchement_codon] + s
					rmax = M
					rseq = new_rseq
					flag = 1
				else : 
					print "Error : les parties chevauchantes ne sont pas identiques.\n"
			if flag == 1 :
				dico[n] = (rmin, rmax, rseq)
		if flag == 0 :
			nb_region += 1
			dico[nb_region] = (m, M, s)
	#verif du dico et transformation en liste
	new_list_regions = []
	for i in dico :
		#print i, dico[i]
		new_list_regions.append(dico[i])
	#print new_list_regions
	return new_list_regions

def regions_chevauchantes_moins(list_regions) :
	"""
	Fonction permettant de mettre ens les regions qui se chevauchent.
	input : 
		- list_regions : liste des regions -pr m frame- [(sstart, send, sseq)]
	output :
		- retourne la meme liste sans regions chevauchantes.
	"""
	rmin = "X"; rmax = "X"; rseq = "X"
	dico = {}
	nb_region = 1
	dico[nb_region] = (rmin, rmax, rseq)
	for (M, m, s) in list_regions :
		m = int(m); M = int(M)
		for n in dico :
			rmin = dico[n][0]
			rmax = dico[n][1]
			rseq = dico[n][2]
			flag = 0
			#initialisation
			if rmin == "X" :
				rmin = m; rmax = M; rseq = s
				flag = 1
			#chevauchement des hits
			##cas 1 : la grande region se trouve dans la region
			elif M <= rmax <= m and M <= rmin <= m :
				rmin = m
				rmax = M
				rseq = s
				flag = 1
			##cas 2 : la region se trouve dans la grande region
			elif rmax <= m <= rmin and rmax <= M <= rmin :
				flag = 1
			##cas 3 : la region touche un bout de la grande region
			elif m >= rmin and rmax <= M <= rmin :
				#print rmin, rmax, m, M
				#print rmin, rmax
				#print rseq
				#print m, M
				#print s
				chevauchement_codon = (rmin - M + 1) / 3
				#print chevauchement_codon
				#print rseq[:chevauchement_codon]
				#print s[-chevauchement_codon:]
				#print s[:-chevauchement_codon]
				#verif que la partie qui se chevauche correspond à la meme seq
				if rseq[:chevauchement_codon] == s[-chevauchement_codon:] :
					new_rseq = s[:-chevauchement_codon] + rseq
					rmin = m
					rseq = new_rseq
					flag = 1
				else : 
					print "Error : les parties chevauchantes ne sont pas identiques.\n"
			if flag == 1 :
				dico[n] = (rmin, rmax, rseq)
		if flag == 0 :
			nb_region += 1
			dico[nb_region] = (m, M, s)
	#verif du dico et transformation en liste
	new_list_regions = []
	for i in dico :
		#print i, (dico[i][1], dico[i][0], dico[i][2])
		new_list_regions.append((dico[i][1], dico[i][0], dico[i][2]))
	#print new_list_regions
	return new_list_regions

def regions_transcrit_per_frame(list_regions_non_chevauchantes) :
	"""
	Fonction permettant de mettre ens les regions d'un transcrit. Ces regions 
	sont non chevauchantes et correspondent au meme frame.
	input : 
		- list_regions_non_chevauchantes : liste de tuples contenant pour un
		frame les differentes regions qui ont été retrouvées avec blast.
	output : 
		- retourne un tuple correspondant à une seq pr un frame.
	"""
	#raccourcissement du nom de la variable
	L = list_regions_non_chevauchantes
	#print L
	h_min = "X"
	h_max = "X"
	h_seq = ""
	for hit in L :
		#print hit
		if h_min == "X" and h_max == "X" :
			h_min = int(hit[0])
			h_max = int(hit[1])
		elif h_max < int(hit[1]) :
			h_max = int(hit[1])
		h_seq += hit[2]
	#print h_min, h_max
	#print h_seq
	return (h_min, h_max, h_seq)

def fasta_files_for_FTRep(transcriptome2_csv_file) :
	"""
	Fonction permettant de creer les fichiers multifasta de chaque esp
	avec les seq sans codant stop obtenues avec la methode approximative.
	selection avec la taille de la seq <=> au min du nb de repet OPR 
	(76 = 38*2).
	input :
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output : 
		- fichier multifasta par esp (transcriptome)
	"""
	#recup du dico
	dico_transcrit = dico_1transcrit_1seq(transcriptome2_csv_file)
	#pr chaque esp
	for ome in dico_transcrit :
		print ome
		#obtention de proteines
		filout = open("FTRep_" + ome + ".faa", "w")
		#pr chaque transcrit
		for transcrit in dico_transcrit[ome] :
			#que les seq obtenues
			if dico_transcrit[ome][transcrit] != "X" :
				#selection par la taille 
				if len(dico_transcrit[ome][transcrit]) >= 76 :
					#print transcrit,
					filout.write(">" + transcrit + "\n")
					filout.write(dico_transcrit[ome][transcrit] + "\n")
		filout.close()

def transcrit_for_exonerate(transcriptome2_csv_file) :
 	"""
 	Fonction permettant de donner les omes et transcrits pr exonerate
 	input : 
 		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output : 
		- retourne un dico qui contient les transcrit à analyser pour chaque 
		ome
 	"""
	#recup du dico
	dico_transcrit = dico_1transcrit_1seq(transcriptome2_csv_file)
	#initialisation dico contenant les transcrits de chaque ome pr exonerate
	dico_exonerate = {}
	#pr chaque esp
	for ome in dico_transcrit :
		#print ome
		dico_exonerate[ome] = []
		#pr chaque transcrit
		for transcrit in dico_transcrit[ome] :
			#print dico_transcrit[ome][transcrit]
			#affichage des omes et transcrits pr exonerate
			if "X" == dico_transcrit[ome][transcrit][-1] :
				#print ome, transcrit
				dico_exonerate[ome].append(transcrit)
	return dico_exonerate

def fasta_creation(transcriptome2_csv_file, rep_fasta_files) :
	"""
	Fonction permettant de creer un fichier fasta
	input : 
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
		- rep_fasta_files : repertoire contenant les fichiers fasta des omes
	output :
		- creation de fichier fasta pr chaque contig voulu
	"""
	#dico contenant les infos des fichiers fasta d'un ome
	fastas = dico_info_fasta(rep_fasta_files)
	#verif du dico
	#for ome in fastas:
	#	print ome
	#	for contig in fastas[ome] :
	#		print contig
	T_exonerate = transcrit_for_exonerate(transcriptome2_csv_file)
	for ome in T_exonerate : 
		#print ome
		for contig in T_exonerate[ome] :
			if ome in fastas :
				print ome, contig
				#print fastas[ome][contig]
				output = ome + "_" + contig + ".fasta"
				filout = open(output, "w")
				filout.write(">" + contig + "\n")
				filout.write(fastas[ome][contig] + "\n")
				filout.close()
				#print fastas[ome][contig]

def return_non_chevauchantes_regions(m_M_list) :
	"""
	Fonction permettant de retourner une liste de regions non chevauchantes.
	input : 
		- m_M_list : liste des start et end des regions
	output :
		- retourne une liste
	"""
	rmin = "X"; rmax = "X"
	dico = {}
	nb_region = 1
	dico[nb_region] = (rmin, rmax)
	for (m, M) in m_M_list :
		m = int(m); M = int(M)
		for n in dico :
			rmin = dico[n][0]
			rmax = dico[n][1]
			flag = 0
			#initialisation
			if rmin == "X" :
				rmin = m; rmax = M
				flag = 1
			#chevauchement des hits
			##cas 1 : la region se trouve dans la grande region
			elif rmin <= m <= rmax and rmin <= M <= rmax :
				flag = 1
			##cas 2 : la region touche un bout de la grande region
			elif rmin <= m <= rmax and M > rmax :
				rmax = M
				flag = 1
			if flag == 1 :
				dico[n] = (rmin, rmax)
		if flag == 0 :
			nb_region += 1
			dico[nb_region] = (m, M)
	#verif du dico et transformation en liste
	new_list_regions = []
	for i in dico :
		#print i, dico[i]
		new_list_regions.append(dico[i])
	#print new_list_regions
	return new_list_regions

def query_for_exonerate(transcriptome2_csv_file) :
	"""
	Fonction permettant de retourner les query des diff transcrits utilisés
	pr exonerate.
	input :
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output :
		- retourne les query. Affiche aussi certaines infos. 
	"""
	#####Enregistrement des données de query
	#transcrits à utiliser pr exonerate
	dico_exonerate = transcrit_for_exonerate(transcriptome2_csv_file)
	#initialisation du dico contenant les infos sur les query pr exonerate
	exonerate_query = {}
	#lecture du fichier resumé de blast
	filin = open(transcriptome2_csv_file, "r")
	#pr chaque ligne du fichier 
	for l in filin :
		#releve des infos
		sep = l.split(" ")
		ome = sep[0]; contig = sep[1]; sstart = int(sep[2]); send = int(sep[3])
		query = sep[6] 
		#si esp ds exonerate
		if ome in dico_exonerate :
			#initialisation de ome si pas ds le dico
			if ome not in exonerate_query :
				exonerate_query[ome] = {}
			#si contig parmi la liste des transcrits pr exonerate
			#region = contig
			#if "|" in region :
			#	region = region.split("|")[-2]
			#contig = region
			if contig in dico_exonerate[ome] :
				#initialisation du contig si pas ds le dico
				if contig not in exonerate_query[ome] :
					exonerate_query[ome][contig] = {}
				#initialisatio du query si pas ds le dico
				if query not in exonerate_query[ome][contig] :
					exonerate_query[ome][contig][query] = []
				#enregistrement des infos utiles
				exonerate_query[ome][contig][query].append((sstart, send))
	filin.close()
	#print exonerate_query

	#####Analyse des données de query pr prendre le meilleur pr 1 transcrit
	#liste contenant tous les query 
	L = []
	#pr chaque esp
	for sp in exonerate_query :
		#pr chaque transcrit de cette esp
		for ctg in exonerate_query[sp] :
			query_chosen = "X"
			max_len_query = 0
			#pr chaque query
			#print "*****"
			#print sp + "_" + ctg, 
			for q in exonerate_query[sp][ctg] :
				#print ">>>", q
				sorted_list = sorted(exonerate_query[sp][ctg][q], key=lambda data: data[0])
				#print sorted_list
				new_sorted_list = return_non_chevauchantes_regions(sorted_list)
				#print new_sorted_list
				s = 0
				#pr chaque elmt de ma liste q, calcul de la couverture du transcrit
				for i in new_sorted_list :
					s += abs(int(i[1])-int(i[0]))
				#print s
				if s > max_len_query :
					max_len_query = s
					query_chosen = q
			#print query_chosen
			if query_chosen not in L :
				L.append(query_chosen)
	return L

def fasta_query_creation(transcriptome2_csv_file, rep_fasta_files) :
	"""
	Fonction permettant de creer un fichier fasta
	input : 
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
		- rep_fasta_files : repertoire contenant les fichiers fasta des omes
	output :
		- creation de fichier fasta pr chaque contig voulu
	"""
	#dico contenant les infos des fichiers fasta d'un ome
	fastas = dico_info_fasta(rep_fasta_files)
	#verif du dico
	#for ome in fastas:
	#	print ome
	#	for contig in fastas[ome] :
	#		print contig
	Q_exonerate = query_for_exonerate(transcriptome2_csv_file)
	for ome in fastas : 
		for contig in fastas[ome] :
			if contig in Q_exonerate :
				print ">" + contig
				#print fastas[ome][contig]
				output = contig + ".fasta"
				filout = open(output, "w")
				filout.write(">" + contig + "\n")
				filout.write(fastas[ome][contig] + "\n")
				filout.close()

def exonerate_cmd(subject_query_file, subject_rep, query_rep) :
	"""
	Fonction permettant de lancer la commande exonerate du modèle coding2coding.
	permet la production de fichiers d'alignements comme resultat. avoir créé 
	un rep results qui contiendra les fichiers créés.
	input :
		- subject_query_file : fichier de deux colonnes avec la premiere 
		contenant le nom des subject et le second le nom des query.
		- subject_rep : nom du rep contenant les fichiers fasta des subjects
		- query_rep : nom du rep contenant les fichiers fasta des querys.
	output :
		- creation des fichiers de'alignements (resultats de exonerate)
	"""
	filin = open(subject_query_file, "r")
	for line in filin : 
		subject = line.split(" ")[0].replace("|", "_")
		query = line.split(" ")[1][:-1].replace("|", "_")
		Sfile = subject_rep + subject + ".fasta"
		Qfile = query_rep + query + ".fasta"
		output = "./results/" + subject + "_" + query + ".res"
		cmd = "/home/mortaza/exonerate-2.2.0-x86_64/bin/exonerate --model coding2coding "
		os.system(cmd + Qfile + " " + Sfile + " > " + output)
		print subject, query
	filin.close()

def comparison_T2_FTRep_Ingrid(transcriptome2_csv_file, rep_FTRep) :
	"""
	Fonction permettant de comparer les transcrits issus du blast2 et ceux
	validés avec FTRep (les ORFs pris en compte sont supposés ne pas avoir 
	d'introns).
	input : 
		- transcriptome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
		- rep_FTRep : rep contenant les resultats de FTRep d'Ingrid
	output :
		- affiche le nom des espces et des transcrits validés par FTRep
	"""
	#recup du dico contenant les infos nec
	dico_hit = dico_hit_T(transcriptome2_csv_file)
	#Liste contenant les transcrits validés par FTRep
	T_from_FTRep_list = []
	files = glob.glob(rep_FTRep + "/*res") + [rep_FTRep + "/Dunaliella-tertiolecta_res2"]
	for file in files :
		#print file
		filin = open(file, "r")
		for line in filin :
			T = line.split("\t")[0][:-2]
			if T not in T_from_FTRep_list :
				T_from_FTRep_list.append(T)
		filin.close()
	print T_from_FTRep_list
	#pr chaque esp
	for ome in dico_hit :
		#pr chaque transcrit
		n = 0
		total = 0
		for contig in dico_hit[ome] :
			#print contig
			total += 1
			#s'il a été validé par FTRep, l'afficher
			if contig in T_from_FTRep_list :
				#print ome, contig
				n += 1
		print ome, n, "/", total


#Programme principal
if __name__ == "__main__":
	#analyse transcriptomes
	csv_files_t_t = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryT/SubjectT/*')
	csv_files_t_g = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryG/SubjectT/*')
	rep_resume1_files = "/data/mortaza/Lucifer/evolution/oprvsgenome/results/genome/regions/resume/part1"
	#lignes de commandes
	l = [(2180, 2314, 'LIEVAQAIQHLHNMKLIHCDIKPENVLLKSDPSKAIGFVTKLSDF'), (2369, 2440, 'GTVTHLAPELFQVGSKLTTAVDTF'), (2444, 2497, 'FGIMMWELYTGQRAYGGL')]
	#regions_chevauchantes_moins(l)
	#fasta_files_for_FTRep(sys.argv[1]) #fasta_files_for_FTRep
	#regions_transcrit_per_frame(l)
	blast1_file = '/data/mortaza/Lucifer/evolution/hitvsome/QueryT/src/transcriptomic_regions.txt'
	#blast2_file = FTRep_OPRcandidates_rep
	#dico_region_blast1(blast1_file)
	#blast1_blast2_comparison(blast1_file, sys.argv[1], sys.argv[2])
	blast2_file = '/data/mortaza/Lucifer/evolution/hitvsome/T_analysis/exonerate/ome_transcrit_for_exonerate.csv'
	#dico_region_blast2(blast2_file)
	#transcrit_for_exonerate(sys.argv[1])
	#fasta_creation(sys.argv[1], sys.argv[2])
	#fasta_query_creation(sys.argv[1], sys.argv[2])
	#exonerate_cmd(sys.argv[1], sys.argv[2], sys.argv[3])
	#comparison_T2_FTRep_Ingrid(sys.argv[1], sys.argv[2])
	exonerate_file_reader(sys.argv[1])

"""
Notes & Commentaires :
- /data/mortaza/src/OPRvsgenome/result_T2_analysis.py
- liste de notre db perso (g, t, p) ds le rep
/data/mortaza/Lucifer/evolution/OPR_evolution/oprvsgenome/results/list
- chemin pour avoir les rep des tailles des diff contigs
/data/mortaza/Data
- cmd pr suppression d'une ligne contenant un certain car ds un fichier
sed -i '/query/d' transcriptome2.csv
- pr utiliser l'executable exonerate :
/home/mortaza/exonerate-2.2.0-x86_64/bin/exonerate
"""