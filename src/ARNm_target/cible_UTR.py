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
sys.path.append('/data/mortaza/src/OPRvsgenome/')
from oprvsgenome import *
from db_building_part2 import *
from OPR_target import *
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import random
#sys.path.append('/data/mortaza/Clotilde_result/')
#from Align_Global import *

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
	filout.write(print_seq(ome, start, end) + "\n")
	filout.close()

def utr_fasta_GW(rep_ome, region_result, type_ome) :
	"""
	Fonction permettant de creer les fichiers fasta des regions utr a partir
	des resultats region_start_end_strand des omes.
	input :
		- rep_ome : repertoire des omes concernés
		- region_result : fichier resultat csv des regions utr
		- type_ome : G (genome) ou W (wgs)
	output :
		- creation de fichier fasta par espece et par utr des genes/cds
	"""
	#liste des fichiers fasta des omes des esp
	omes = glob.glob(rep_ome + "/*.fasta")
	#lecture du fichier resultat
	filin = open(region_result, "r")
	#ne pas lire le header
	line = filin.readline()
	#pour chaque ligne du fichier
	for line in filin :
		#print "***"
		#print line[:-1]
		#recup des infos
		#Si G ou W
		species = line.split(" ")[0]
		strand = line.split(" ")[-1][:-1] 
		if type_ome == "G" :
			if strand == "+" :
				gene = line.split(" ")[2] + "_" + strand
				#print gene
			else :
				gene = line.split(" ")[1] + "_" + strand
				#print gene
		else :
			gene = line.split(" ")[1]
		start = line.split(" ")[3]
		end = line.split(" ")[4]
		#si le start est bien avant le end (verif car parfois des erreurs)
		if int(start) < int(end) :
			#print line[:-1]
			#print gene
			#print species, gene, start, end
			for ome in omes :
				if species in ome :
					region_fasta_creation(ome, start, end, gene)
	filin.close()

def utr_fasta_T(rep_ome, region_result) :
	"""
		Fonction permettant de creer les fichiers fasta des regions utr a partir
	des resultats region_start_end_strand des omes.
	input :
		- rep_ome : repertoire des omes concernés
		- region_result : fichier resultat csv des regions utr transcriptomiques
	output :
		- creation de fichier fasta par espece et par utr des genes/cds
	"""
	#dico contenant les infos sur contigs des transcriptomes
	dico = dico_info_fasta(rep_ome)
	#verif du dico
	#for ome in dico :
	#	print ome
	#	for ctg in dico[ome] :
	#		print ctg
	#		print dico[ome][ctg]
	#lecture du fichier resultat
	filin = open(region_result, "r")
	#ne pas lire le header
	line = filin.readline()
	#pour chaque ligne du fichier
	for line in filin :
		#print "***"
		#print line[:-1]
		#recup des infos
		species = line.split(" ")[0]
		contig = line.split(" ")[1]
		gene = line.split(" ")[2]
		start = line.split(" ")[3]
		end = line.split(" ")[4]
		strand = line.split(" ")[-1][:-1]
		header = ">" + gene + "_" + strand + "_" + start + "_" + end + "_" + \
		species
		for ome in dico :
			if species in ome :
				filout = open(header[1:] + ".fasta", "w")
				filout.write(header + "\n")
				#print dico[ome][contig][int(start)-1:int(end)]
				filout.write(dico[ome][contig][int(start)-1:int(end)] + "\n")
				filout.close()
	filin.close()

def footprint_dico(footprint_seq_rep) :
	"""
	Fonction permettant de retourner un dico cle = cible et valeur = seq
	input :
		- footprint_seq_rep : rep des fichiers fasta des seq du footprint 
	output :
		- retourne dico
	"""
	#initialisation du dico
	dico_f = {}
	#reĉup de tous les fichiers fasta des footprints
	footprint_seq = glob.glob(footprint_seq_rep + "/*fasta")
	#pour chaque fichier fasta de footprint
	for footprint in footprint_seq :
		filin = open(footprint, "r")
		for line in filin : 
			#si header, prendre le nom de la cible
			if '>' in line :
				f = line[2:-1]
				#print f
			#sinon prendre la seq de la cible
			else :
				f_seq = line[:-1]
				#print f_seq
		#enregistrement des données dans le dico
		dico_f[f] = f_seq
		filin.close()
	return dico_f

def utr_dico(species_utr_sequence_rep) :
	"""
	Fonction permettant de retourner un dico cle = cible et valeur = seq
	input :
		- species_utr_sequence_rep : rep des fichiers fasta des seq 5'UTR 
		d'une esp
	output :
		- retourne dico
	"""
	#initialisation du dico
	dico_utr = {}
	#reĉup de tous les fichiers fasta des utr des especes
	utr_species = glob.glob(species_utr_sequence_rep + "/*fasta")
	#pour chaque fichier fasta d'especes
	for utr in utr_species :
		filin = open(utr, "r")
		for line in filin : 
			#si header, prendre le nom de la cible
			if '>' in line :
				utr_name = line[1:-1]
				#print utr_name
			#sinon prendre la seq de la cible
			else :
				utr_seq = line[:-1]
				#print utr_seq
		#enregistrement des données dans le dico
		dico_utr[utr_name] = utr_seq
		filin.close()
	return dico_utr

def Align_Global_Z_I (seq1,seq2,I=100):
	"""
	A partir de deux séquences proteique retourne le z-score et le pourcentage
	 identité de leur alignement global. Par défaut 100 séquences aléatoire. 
	Paramètres : BLOSUM 50 /  gap-open = -12 / extend-gap = -2 / free gap end
	"""
	matrix = matlist.blosum50
	S=[] # Initialisation liste des scores des alignements seq1 et seq2 mélangé
	#aléatoirement
	for i in range(I):
		S.append(pairwise2.align.globalds(seq1, \
			''.join(random.sample(seq2,len(seq2))),matrix,-12,-2, \
			penalize_end_gaps=(False, False))[0][2])
	Align = pairwise2.align.globalds(seq1, seq2, matrix, -12, -2)
	s = Align[0][2]
	if np.std(S) == 0:
		print seq1,seq2,S
	#else :
		#print "ok"
	z = (s-np.mean(S))/np.std(S)
	print z
	#i = Identity(Align[0])
	#return z
	return ' '

def global_alignment(footprint_seq_rep, species_utr_sequence_rep) :
	"""
	Fonction permettant de faire un alignement global par programmation 
	dynamique
	input :
		- footprint_seq_rep : rep des fichiers fasta des seq du footprint
		- species_utr_sequence_rep : rep des fichiers fasta des seq 5'UTR 
		d'une esp
		- cible : nom du cible (mRNA)
	output :
		retourne l'alignement global
	"""
	cible = species_utr_sequence_rep.split("/")[-1] + "_f"
	dico_f = footprint_dico(footprint_seq_rep)
	dico_utr = utr_dico(species_utr_sequence_rep)
	footprint = dico_f[cible]
	#print cible,
	for utr in dico_utr :
		utr_species = dico_utr[utr]
		#print utr
		#print Align_Global_Z_I (utr_species,footprint,I=100)
		alignments = pairwise2.align.globalxx(utr_species,footprint)
		print alignments[0]
		#print format_alignment(*alignments[0])

#Programme principal
if __name__ == "__main__":
	global_alignment(sys.argv[1], sys.argv[2]) 


"""
Notes & Commentaires :
- /data/mortaza/src/ARNm_target/cible_UTR.py
- emplacement des fichiers de taille des contigs pr les transcriptomes :
/data/mortaza/Data/transcriptomes/transcriptome_len
"""