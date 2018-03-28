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
import re
sys.path.append('/home/mortaza/Documents/src/OPRvsgenome/')
from oprvsgenome import *
from OPR_target import *

#Fonctions

def return_dico_cibles_wgs(cibles_cds_file) :
	"""
	Fonction permettant de retourner un dico contenant les infos cds des cibles
	input :
		- cibles_cds_file : fichier de resultats cds des cibles
	output :
		- retourne un dico
	"""
	#initialisation du dico
	dico_cible = {}
	filin = open(cibles_cds_file, "r")
	#ne pas lire le header
	line = filin.readline()
	#pour chaque ligne du fichier
	for line in filin :
		species = line.split(",")[0].replace(" ","")
		cds = line.split(",")[1].replace(" ","")
		frame = line.split(",")[4].replace(" ","")
		if frame == "-1" or frame == "-2" or frame == "-3" or frame == "'-'" :
			strand = "-"
		else :
			strand = "+"
		start = line.split(",")[5]
		end = line.split(",")[6]
		if species not in dico_cible :
			if species != "Casymmetrica_atpA_contig" or \
			species != "Squadricauda_atpA_contig" :
				#pour ne prendre en compte que les resultats de WGS
				if "contig" in species :
					dico_cible[species] = {}
		if species in dico_cible :
			dico_cible[species][cds] = (strand, int(start), int(end))
	filin.close()
	#verif dico
	#for cle in dico_cible : 
	#	print cle 
	#	for cle2 in dico_cible[cle] :
	#		print cle2, dico_cible[cle][cle2]
	return dico_cible

def intergenic_regions_wgs(cibles_cds_file, amonts_cds_file) :
	"""
	Fonction permettant donner les coordonnées des regions intergeniques
	cibles - amonts.
	input :
		- cibles_cds_file : fichier de resultats cds des cibles
		- amonts_cds_file : fichier de resultats cds des amonts
	output :
		- retourne les coordonnées des regions intergeniques
	"""
	#Dico de correspondance entre les cibles et les amonts
	dico_corr_amont_cible = {"atpB":"ORF1995","atpA":"rbcL","ccsA":"trnL2", \
	"petD":"petA","petG":"psbL","psaA":"rps8","psaB":"trnQ","psbC":"ORF2971",\
	"psbI":"atpA"}
	#Dico cibles
	dico_cible = return_dico_cibles_wgs(cibles_cds_file)
	#Dico cibles
	dico_amont = return_dico_cibles_wgs(amonts_cds_file)
	#pour chaque espece de la cible
	print "species gene1 gene2 end1 start2 intergenic_region strand"
	for species in dico_cible :
		#si l'espece existe dans le dico amont
		if species in dico_amont :
			#pour les cds psts ds dico  cibles
			for cds in dico_cible[species] :
				#si le cds amont (acds) est present dans le dico amont
				acds = dico_corr_amont_cible[cds]
				if acds in dico_amont[species] :
					#print species
					#print cds, acds,
					#comparaison des strands avec s1 pr cible et s2 pr amont
					s1 = dico_cible[species][cds][0]
					s2 = dico_amont[species][acds][0]
					#print s1, s2, 
					#si les cds sur meme strand
					if s1 == s2 :
						#start et end des cibles (1) et amonts (2)
						info1 = (dico_cible[species][cds][1], \
						dico_cible[species][cds][2])
						info2 = (dico_amont[species][acds][1], \
						dico_amont[species][acds][2])
						print species, cds, acds,
						if s1 == "+" :
							print info1[0], info2[1],
							print info1[0] - info2[1],
						else :
							print info1[1], info2[0],
							print info2[0] - info1[1],
						print s1

def dico_gene_infos_from_genbank(genbank_file) :
	"""
	Fonction permettant de retourner un dico contenant les infos sur les genes
	de l'organisme en question.
	input :
		- genbank_file : fichier genbank
	output :
		- retourne un dico avec les infos des genes
	"""
	#initialisation de listes de genes sur le brin + ou -
	list_plus = []
	list_moins = []
	filin = open(genbank_file, "r")
	#type de ligne qu'on veut avoir
	regex_pos = re.compile(r"\s+gene(.)+")
	regex_name = re.compile(r"\s+/gene=(.)+")
	#initialisation du flag 
	flag = "X"
	#pour chaque ligne de ce fichier Genbank
	for line in filin :
		if regex_pos.match(line) : 
			positions = line.split("gene")[1].replace(" ","")[:-1]
			#pas encore eu le nom du gene
			flag = 0
		if regex_name.match(line) :
			if flag == 0 :
				gene = line.split('"')[1]
				#eu le nom du gene (eviter d'avoir des repetitions)
				flag = 1
				if "complement" in positions :
					positions = positions.split("(")[1].split(")")[0]
					#print gene, positions
					list_moins.append((gene, positions))
				else :
					#print gene, positions
					list_plus.append((gene, positions))
	filin.close()
	#initialisation d'un dico contenant les deux listes
	dico_list = {"+" : list_plus, "-" : list_moins}
	#Liste des cibles
	list_cible = ["atpB", "atpA", "ccsA", "petD", "petG", "psaB", "psbC", \
	"psbI"]
	#dico pour ne retourner que les resultats nous concernant
	true_dico = {"+" : [], "-" : []} 
	#ttt du dico pour ne prendre que les cibles et les amonts
	for i in range(len(dico_list["+"])) :
		if dico_list["+"][i][0] in list_cible :
			#print dico_list["+"][i-1], dico_list["+"][i]
			true_dico["+"].append((dico_list["+"][i-1], dico_list["+"][i]))
	for i in range(len(dico_list["-"])) :
		if dico_list["-"][i][0] in list_cible :
			#print dico_list["-"][i], dico_list["-"][i+1]
			true_dico["-"].append((dico_list["-"][i], dico_list["-"][i+1]))
	#verif du dico
	#for i in true_dico :
	#	print i
	#	print true_dico[i]
	return true_dico

def intergenic_regions_genome(rep_gb_files) :
	"""
	Fonction permettant donner les coordonnées des regions intergeniques
	cibles - amonts.
	input :
		- rep_gb_files : repertoire des fichiers Genbank des especes 
		ayant genome chloroplastique complet
	output :
		- retourne les coordonnées des regions intergeniques
	"""
	print "species gene1 gene2 end1 start2 intergenic_region strand"
	files = glob.glob(rep_gb_files + "/*.gb")
	for file in files :
		dico_gene = dico_gene_infos_from_genbank(file)
		species = file.split("/")[-1][:-3]
		for sign in dico_gene :
			#print sign
			for genes in dico_gene[sign] :
				#print genes
				#noms de genes
				gene1 = genes[0][0]
				gene2 = genes[1][0]
				#start et end de la region intergenique
				end1 = genes[0][1].split("..")
				#enlever les join
				if "join" not in end1 :
					end1 = end1[1]
					start2 = genes[1][1].split("..")[0]
					intergenic_region = int(start2) - int(end1)
					#affichage du resultat
					print species, gene1, gene2, end1, start2, intergenic_region, sign

def utr_regions_transcriptome(cds_transcriptome_file, rep_Tcontig_length) :
	"""
	Fonction permettant de retourner la partie utr des cds transcriptomiques
	input :
		- cds_transcriptome_file : fichier resumé des resultat transcriptomi-
		ques blast
		- rep_Tcontig_length : nom du rep contenant les fichiers de taille des
		contigs pour chaque transcriptome
	output :
		- retourne les regions utr
	"""
	print "species contig cds start end strand"
	#dico de taille des contigs
	dico_contigs = dico_len_contigs(rep_Tcontig_length)
	#for i in dico_contigs :
	#	print i
	#lecture du fichier resumé des cds des cibles
	filin = open(cds_transcriptome_file, "r")
	#ne pas lire le header
	line = filin.readline()
	#pour chaque ligne
	for line in filin :
		#print line[:-1]
		#recuperer les infos imp
		species = line.split(",")[0].replace(" ", "")
		contig =  line.split(",")[-1][2:-2]
		cds = line.split(",")[1].replace(" ", "")
		frame = int(line.split(",")[4])
		if frame < 0 :
			strand = "-"
		else :
			strand = "+"
		start = line.split(",")[5].replace(" ", "")
		end = line.split(",")[6].replace(" ", "")
		#ttt selon le strand
		#print species, contig, cds,
		if strand == "+" :
			print species, contig, cds, "1", start, strand
		else :
			print species, contig, cds, end, \
			dico_contigs[species][contig], strand
	filin.close()

def region_fasta_creation(ome, start, end, region) :
	"""
	Fonction permettant de créer un fichier fasta d'une region d'un ome.
	input :
		- genome : fichier fasta du genome / transcriptome / contig wgs
		- start : position start de la region concernée dans le ome
		- end : position end de la region concernée dans le ome
		- region : nom de la region (correspondant au nom du gene par ex)
	output :
		- creation d'un fichier fasta
	"""
	ome_name = ome.split('/')[-1].split(".")[0]
	header = ">" + region + "_" + start + "_" + end + "_" + ome_name
	output = region + "_" + ome_name + ".fasta"
	filout = open(output, "w")
	#print header
	#print print_seq(ome, start, end)
	filout.write(header + "\n")
	filout.write(print_seq(ome, start, end))
	filout.close()

def utr_fasta(rep_ome, region_result, type_ome) :
	"""
	Fonction permettant de creer les fichiers fasta des regions utr a partir
	des resultats region_start_end_strand des omes.
	input :
		- rep_ome : repertoire des omes concernés
		- region_result : fichier resultat csv des regions utr
		- type_ome : G (genome) ou T (transcriptome) ou W (wgs)
	output :
	"""
	omes = glob.glob(rep_ome + "/*.fasta")
	for ome in omes :
		print ome


#Programme principal
if __name__ == "__main__":
	utr_fasta(sys.argv[1], sys.argv[2], sys.argv[3])



"""
Notes & Commentaires :
- /home/mortaza/Documents/src/ARNm_target/cible_UTR.py
- emplacement des fichiers de taille des contigs pr les transcriptomes :
/data/mortaza/Data/transcriptomes/transcriptome_len
"""