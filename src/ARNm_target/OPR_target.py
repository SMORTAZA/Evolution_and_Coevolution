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
sys.path.append('/data/mortaza/src/OPRvsOPR_rep/')
from cluster_quality import *

#Fonctions

def print_seq(ome, start, end) :
	"""
	Fonction permettant de retourner la région voulue du ome
	input :
		- ome : fichier fasta du ome
		- start : position start de la region concernée dans le ome
		- end : position end de la region concernée dans le ome
	output :
		- retourne la region voulue
	"""
	seq = ''
	filin = open(ome, 'r')
	for l in filin :
		if '>' not in l :
			#print l[:-1]
			seq = seq + l[:-1]
	filin.close()
	#print len(seq)
	#print seq[int(start):int(end)]
	return seq[int(start) - 1 :int(end)]

def psaA_exons(psaA_fasta = "/data/mortaza/Lucifer/coevolution/Chre_genes/psaA.fasta") :
	"""
	Fonction permettant de créer les fichiers fasta des trois exons de psaA à partir du fichier fasta de psaA from NCBI.
	exon = true_position <=> position_fasta
	exon1 = 32737-32825 <=> 1-89
	exon2 = 69520-71506 <=> 36784-38770
	exon3 = 174205-174384 <=> 141469-141648
	input :
		- psaA_fasta : fichier fasta du gene psaA
	output : 
		- trois fichiers fasta des exons de psaA
	"""
	#chaine de car initialisée qui contiendra la seq complete du gene psaA
	seq = ""
	filin = open(psaA_fasta, "r")
	for line in filin :
		if '>' not in line :
			seq = seq + line[:-1]
	filin.close()
	exon1 = open('psaA1.fasta', 'w')
	exon1.write('> psaA1\n' + seq[0:89])
	exon1.close()
	exon2 = open('psaA2.fasta', 'w')
	exon2.write('> psaA2\n' + seq[36784:38770])
	exon2.close()
	exon3 = open('psaA3.fasta', 'w')
	exon3.write('> psaA3\n' + seq[141469:])
	exon3.close()

def print_species_target_cds(rep_files, th_Evalue = 1e-03) :
	"""
	Fonction permettant de donner les cds des genes OPR_target.
	input :
		- rep_files : nom rep contenant les fichiers resultats csv
		- th_Evalue : seuil de Evalue
	output :
		- affichage des résultats voulues à l'écran
	"""
	dico_cds = {}
	#liste des fichiers de resultats contenu ds le rep en question
	files = glob.glob(rep_files + "/*csv")
	#pour chaque fichier (par esp)
	for file in files :
		filin = open(file, 'r')
		for line in filin :
			if float(line.split('\t')[16]) <= th_Evalue :
				#print float(line.split('\t')[16])
				cds = file.split('/')[-1].split('vs')[0].replace("_prot","")
				#print cds
				species = file.split('/')[-1].split('vs')[1].split('.')[0]
				#print species
				if species not in dico_cds : 
					dico_cds[species] = [cds]
				elif cds not in dico_cds[species] :
					dico_cds[species].append(cds)
		filin.close()
		
	#verif dico
	for espece in dico_cds :
		print espece, dico_cds[espece]

def dico_species_target_cds_genome(rep_files, prot_len, th_Evalue = 1e-03) :
	"""
	Fonction permettant de donner les regions cds des genes OPR_target. specifiée pour genome.
	input :
		- rep_files : nom rep contenant les fichiers resultats csv
		- th_Evalue : seuil de Evalue
	output :
		- retourn 
	"""
	#dico de taille des prot
	dico_prot_len = {}
	filin = open(prot_len, 'r')
	for line in filin : 
		if "gene" in line :
			prot = line.split("_gene")[0]
		else :
			prot = line.split("_prot")[0]
		prot_len = float(line.split("\t")[1])
		dico_prot_len[prot] = prot_len
	filin.close()
	#verif du dico
	#for i in dico_prot_len :
	#	print i, dico_prot_len[i]
	#dico contenant les infos
	dico_cds = {}
	#liste des fichiers de resultats contenu ds le rep en question
	files = glob.glob(rep_files + "/*csv")
	#pour chaque fichier (par esp), remplissage du dico
	for file in files :
		filin = open(file, 'r')
		for line in filin :
			Evalue = float(line.split('\t')[16]) 
			if Evalue <= th_Evalue :
				#print Evalue
				if "gene" in file.split('/')[-1].split('vs')[0] :
					cds = file.split('/')[-1].split('vs')[0].replace("_gene","")
				else :
					cds = file.split('/')[-1].split('vs')[0].replace("_prot","")
				#print cds
				species = file.split('/')[-1].split('vs')[1].split('.')[0]
				#print species
				qstart = int(line.split('\t')[2])
				qend = int(line.split('\t')[3])
				sstart = int(line.split('\t')[6])
				send = int(line.split('\t')[7])
				sframe =  int(line.split('\t')[9])
				contig = line.split('\t')[4]
				cds_percentage = round(100 * (qend - qstart) / dico_prot_len[cds], 1)
				#print cds, species, Evalue, sstart, send, sframe, contig
				if species not in dico_cds : 
					dico_cds[species] = {}
				if cds not in dico_cds[species] :
					dico_cds[species][cds] = []
				if send < sstart :
					tmp = sstart
					sstart = send
					send = tmp
				#print sframe, sstart, send
				#dico_cds[species][cds].append((sstart, send, sframe))
				dico_cds[species][cds].append((qstart, qend, sframe, sstart, send, Evalue, cds_percentage))
		filin.close()	
	"""
	#verif dico + resume regions
	print 'species\tcds\tnb\tqstart\tqend\tsframe\tsstart\tsend\tEvalue\t%CDS'
	for espece in dico_cds :
		#print 
		#print espece
		for c in dico_cds[espece] :
			#print c, 
			#print sorted(dico_cds[espece][c], key=lambda data: data[0])
			n = -1
			for i in sorted(dico_cds[espece][c], key=lambda data: data[0]) :
				n += 1 
				print espece + '\t' + c + '\t' + str(n), 
				print i
			#regions = resume_regions(sorted(dico_cds[espece][c], key=lambda data: data[0]))
			#print regions
			#for region in dico_cds[espece][c] :
				#print region
	"""
	return dico_cds
	
def dico_species_target_cds_contig(rep_files, prot_len, th_Evalue = 1e-03) :
	"""
	Fonction permettant de donner les regions cds des genes OPR_target. 
	quand plusieurs contigs (transcriptomes, wgs)
	input :
		- rep_files : nom rep contenant les fichiers resultats csv
		- th_Evalue : seuil de Evalue
	output :
		- retourn 
	"""
		#dico de taille des prot
	dico_prot_len = {}
	filin = open(prot_len, 'r')
	for line in filin : 
		if "gene" in line :
			prot = line.split("_gene")[0]
		else :
			prot = line.split("_prot")[0]
		prot_len = float(line.split("\t")[1])
		dico_prot_len[prot] = prot_len
	filin.close()
	#dico contenant les infos
	dico_cds = {}
	#liste des fichiers de resultats contenu ds le rep en question
	files = glob.glob(rep_files + "/*csv")
	#pour chaque fichier (par esp), remplissage du dico
	for file in files :
		filin = open(file, 'r')
		for line in filin :
			Evalue = float(line.split('\t')[16]) 
			if Evalue <= th_Evalue :
				#print Evalue
				if "gene" in file.split('/')[-1].split('vs')[0] :
					cds = file.split('/')[-1].split('vs')[0].replace("_gene","")
				else :
					cds = file.split('/')[-1].split('vs')[0].replace("_prot","")
				#print cds
				species = file.split('/')[-1].split('vs')[1].split('.')[0]
				#print species
				qstart = int(line.split('\t')[2])
				qend = int(line.split('\t')[3])
				sstart = int(line.split('\t')[6])
				send = int(line.split('\t')[7])
				sframe =  int(line.split('\t')[9])
				contig = line.split('\t')[4]
				#print cds, species, Evalue, sstart, send, sframe, contig
				if species not in dico_cds : 
					dico_cds[species] = {}
				if cds not in dico_cds[species] :
					dico_cds[species][cds] = []
				if send < sstart :
					tmp = sstart
					sstart = send
					send = tmp
				#print sframe, sstart, send
				cds_percentage = round(100 * (qend - qstart) / dico_prot_len[cds], 1)
				#print (qstart, qend, sframe, sstart, send, Evalue, cds_percentage, contig)
				dico_cds[species][cds].append((qstart, qend, sframe, sstart, send, Evalue, cds_percentage, contig))
		filin.close()
	
	#ttt dico
	#Dico contenant la meilleure region pour le cds
	CDS = {}
	#pour chaque esp
	for espece in dico_cds :
		#print espece
		CDS[espece] = {}
		#pour chaque cds
		for c in dico_cds[espece] :
			#print c
			#print dico_cds[espece][c]
			#print "***", len(dico_cds[espece][c])
			#si taille de la liste != de 1
			if len(dico_cds[espece][c]) != 1 :
				#initialisation des valeurs_region
				m = 0
				M = 0
				f = 'X'
				ctg = 'X'
				CDS[espece][c] = (m, M, f, ctg)
				#pour chaque elmt / region de la liste
				for r in dico_cds[espece][c] :
					#print r,
					#si taille (send-sstart) > taille (M-m)
					if int(r[1]) - int(r[0]) > M - m :
						m = r[3]; M = r[4]; f = r[2]; ctg = r[-1]
						qs = r[0]; qe = r[1]; ev = r[5]; cp = r[6]					
						#best region
						CDS[espece][c] = [(qs, qe, f, m, M, ev, cp, ctg)]
			else :
				CDS[espece][c] = dico_cds[espece][c]
			#print "==>", CDS[espece][c]
	print 'species\tcds\tqstart\tqend\tsframe\tsstart\tsend\tEvalue\t%CDS\tcontig'
	#verif dico modifié
	for espece in CDS :
		#print espece
		for c in CDS[espece] :
			print espece + '\t' + c, CDS[espece][c][0]
	
def cds_regions(rep_files, prot_len, th_Evalue = 1e-03) :
	"""
	Fonction permettant de retourner les vraies regions cds (1/cds)
	"""
	#target_region = {}
	print 'species\tcds\tqstart\tqend\tsframe\tsstart\tsend\tEvalue\t%CDS'
	dico_cds = dico_species_target_cds_genome(rep_files, prot_len, th_Evalue)
	for espece in dico_cds :
		#target_region[espece] = {}
		for c in dico_cds[espece] :
			#si resultat unique
			if len(dico_cds[espece][c]) == 1 :
				#sauf pour ces deux cibles
				if not((c == "atpA" or c == "atpB") and dico_cds[espece][c][0][-1] < 95) :
					print espece + '\t' + c,
					print dico_cds[espece][c][0]
			#si plusieurs resultats
			else :
				flag = 0
				#si un des resultats est le meilleur
				for i in sorted(dico_cds[espece][c], key=lambda data: data[0]) :
					if (float(i[-2]) == 0 or float(i[-1]) >= 95.0) and flag == 0 :
						flag = 1
						print espece + '\t' + c,
						print i
				#si pas de choix
				if flag == 0 :
					rep = cds_choice(espece, c, sorted(dico_cds[espece][c], key=lambda data: data[0]), prot_len) 
					print espece + '\t' + c, 
					print rep

					"""
					n = 0
					for i in sorted(dico_cds[espece][c], key=lambda data: data[0]) :
						n += 1
						print espece + '\t' + c + '\t' + str(n), 
						print i
					"""

def cds_choice(species, cds, cds_liste, prot_len) :
	"""
	Fonction permettant de choisir la meilleure region cds parmi plusieurs regions
	input :
		- species : nom de l'esp
		- cds : nom du cds
		- cds_liste : liste contenant les infos sur les differentes regions
	output :
		- retourne les coordonnées de la meilleure région
	"""
	#dico de taille des prot
	dico_prot_len = {}
	filin = open(prot_len, 'r')
	for line in filin : 
		if "gene" in line :
			prot = line.split("_gene")[0]
		else :
			prot = line.split("_prot")[0]
		prot_len = float(line.split("\t")[1])
		dico_prot_len[prot] = prot_len
	filin.close()
	#valeurs brins +
	pqstart = pqend = psstart = psend = pevalue = "X"
	#valeurs brins -
	mqstart = mqend = msstart = msend = mevalue = "X"
	for i in sorted(cds_liste, key=lambda data: data[0]) :
		#print "***\n",i, "\n***"
		qstart = int(i[0]); qend = int(i[1]); sstart = int(i[3]); send = int(i[4]); evalue = float(i[5])
		#si region sur brin - 
		if int(i[2]) < 0 :
			#initialisation des variables avec les vraies valeurs
			if mqstart == "X" :
				mqstart = qstart
				mqend = qend
				msstart = sstart
				msend = send
				mevalue = evalue
			if sstart < msstart : msstart = sstart
			if send > msend : msend = send
			if qstart < mqstart : mqstart = qstart
			if qend > mqend : mqend = qend
			if evalue < mevalue : mevalue = evalue
		#si region sur brin +
		else :
			#initialisation des variables avec les vraies valeurs
			if pqstart == "X" :
				pqstart = qstart
				pqend = qend
				psstart = sstart
				psend = send
				pevalue = evalue
			if sstart < psstart : psstart = sstart
			if send > psend : psend = send
			if qstart < pqstart : mqstart = qstart
			if qend > pqend : pqend = qend
			if evalue < pevalue : pevalue = evalue
	#affichage à l'écran
	if mqstart != 'X' and pqstart != 'X' :
		if mevalue < pevalue :
			cds_percentage = round(100 * (mqend - mqstart) / dico_prot_len[cds], 1)
			return (mqstart, mqend, '-', msstart, msend, mevalue, cds_percentage)
		else :
			cds_percentage = round(100 * (pqend - pqstart) / dico_prot_len[cds], 1)
			return (pqstart, pqend, '+', psstart, psend, pevalue, cds_percentage)
	elif mqstart != 'X' and pqstart == 'X' :
		cds_percentage = round(100 * (mqend - mqstart) / dico_prot_len[cds], 1)
		return (mqstart, mqend, '-', msstart, msend, mevalue, cds_percentage)
	elif pqstart != 'X' and mqstart == 'X' :
		cds_percentage = round(100 * (pqend - pqstart) / dico_prot_len[cds], 1)
		return (pqstart, pqend, '+', psstart, psend, pevalue, cds_percentage)


#Programme principal
if __name__ == "__main__":
	#print_seq(sys.argv[1], sys.argv[2], sys.argv[3])
	#*****Analyse*****
	#a lancer ds le rep result
	#dico_species_target_cds_genome(sys.argv[1], sys.argv[2])
	cds_regions(sys.argv[1], sys.argv[2])
	#dico_species_target_cds_contig(sys.argv[1], sys.argv[2])
	#/home/mortaza/Documents/src/ARNm_target/OPR_target.py .
	#seq_length(sys.argv[1]) # calcul de la taille des prot



"""
Notes & Commentaires :
- /data/mortaza/src/ARNm_target/OPR_target.py
- remplacer la ligne header de chaque gene par son vrai nom 
sed -ir 's/^>.*c161640-160165.*$/>atpB/g' Chre_target_CDS.fasta
- /data/mortaza/Lucifer/coevolution/Chre_prot/chre_prot_len.txt
"""