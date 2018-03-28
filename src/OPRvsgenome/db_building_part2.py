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
#from matplotlib import pyplot as plt
#import numpy as np

#Fonctions

def dico_info_resume1_file(rep_resume1_files) :
	"""
	Fonction permettant de resumer ds un dico les infos necessaires pour d'autres fct pour chaque espece
	input : 
		- rep_resume1_files : repertoire contenant les fichiers resume1 des especes qu'on etudie (genomes)
	output :
		- retourne un dico
	"""
	files = glob.glob(rep_resume1_files + "/*.csv")
	genomic_hit = {}
	for file in files :
		#print file
		genome = file.split("resume1_")[1]
		#print genome
		if genome not in genomic_hit :
			genomic_hit[genome] = {}
		filin = open(file, "r")
		#ne pas lire le header des fichiers
		line = filin.readline()
		for line in filin :
			contig = line.split("\t")[0]
			contig_length = int(line.split("\t")[1])
			region = line.split("\t")[2]
			region_length = float(line.split("\t")[3][:-2])
			#print contig, contig_length, region, region_length
			if contig not in genomic_hit[genome] :
				genomic_hit[genome][contig] = (contig_length, region, region_length)
		filin.close()
	return genomic_hit


def G_contig(rep_resume1_files) :
	"""
	Fonction permettant de relever l'ens des contigs pst ds resultat de blast contre OPR_Chre
	input : 
		- rep_resume1_files : repertoire contenant les fichiers resume1 des especes qu'on etudie (genomes)
	output :
		- retourne l'ens des contigs retenus
	"""
	dico_region = dico_info_resume1_file(rep_resume1_files)
	db_from_blast = {}
	for genome in dico_region :
		db_from_blast[genome] = []
		for contig in dico_region[genome] :
			#print dico_region[genome][contig]
			contig_length = int(dico_region[genome][contig][0])
			sstart = int(dico_region[genome][contig][1].split("-")[0])
			send = int(dico_region[genome][contig][1].split("-")[1])
			region_length = float(dico_region[genome][contig][2])
			#si la region trouvée est inf à 50kb, alors on considere que c'est une region d'un gene
			if region_length < 50 :
				#print contig_length, sstart, send, region_length
				#start et end du "gene"
				start = sstart - 50000
				end = send + 50000
				if start < 1 :
					start = 1
				if end > contig_length :
					end = contig_length
				#print contig, start, end
				db_from_blast[genome].append((contig,start,end))
	#verif du dico
	for g in db_from_blast :
		print '>' + g[:-4] #sans le csv
		for r in db_from_blast[g] :
			#contig, start et end
			print str(r[0]) + '\t' + str(r[1]) + '\t' + str(r[2])
	return db_from_blast


def PandT_contig(rep_resume1_files) :
	"""
	Fonction permettant de relever l'ens des contigs / proteines pst dans resultat de blast contre OPR_Chre
	input :
		- rep_resume1_files : rep contenant les fichiers resume1 des especes qu'on etudie (proteomes et transcriptomes)
	output :
		- retourne l'ens des contigs retenus
	"""
	files = glob.glob(rep_resume1_files + "/*.csv")
	db_from_blast = {}
	for file in files :
		#print file
		ome = file.split("resume1_")[1]
		#print genome
		db_from_blast[ome] = []
		filin = open(file, "r")
		#ne pas lire le header des fichiers
		line = filin.readline()
		for line in filin :
			contig = line.split("\t")[0]
			if contig not in db_from_blast[ome] :
				db_from_blast[ome].append(contig)
		filin.close()
	#verif dico
	for i in db_from_blast :
		print '>', i[:-4]
		for j in db_from_blast[i] :
			print j
	return db_from_blast


def dico_info_fasta(rep_fasta_files) :
	"""
	Fonction permettant de retourner un dictionnaire avec comme cle le nom du fichier fasta et comme valeur le contig et sa seq
	input :
		- rep_fasta_files : repertoire contenant l'ens des fichiers fasta d'un ome 
	output :
		- retourne un dico contenant les infos des fichiers fasta
	"""
	dico_fasta = {}
	fasta_files = glob.glob(rep_fasta_files + "/*.faa") + \
	glob.glob(rep_fasta_files + "/*.fsa_nt") + \
	glob.glob(rep_fasta_files + "/*.fasta") + \
	glob.glob(rep_fasta_files + "/*.nt.fa") + \
	glob.glob(rep_fasta_files + "/*.fna")
	for file in fasta_files :
		#print file
		species = file.split("/")[-1]
		print species
		dico_fasta[species] = {}
		fasta = open(file, "r")
		contig = ""
		#initialisation 
		flag = "X"
		for line in fasta :
			if ">" in line :
				#initialisation
				if flag == "X" :
					contig_name = line.split()[0][1:]
					flag = 0
				if contig != "" : # pour ne pas mettre la seq vide
					#print contig_name
					#print contig
					dico_fasta[species][contig_name] = contig
					contig_name = line.split()[0][1:]
					contig = ""
					#print contig_name
					#print contig
			if ">" not in line :
				contig = contig + line[:-1]
		#print contig_name
		#print contig
		dico_fasta[species][contig_name] = contig
		fasta.close()
	#verif dico
	#for s in dico_fasta :
		#print s
		#for c_n in dico_fasta[s] :
			#print c_n
			#print dico_fasta[s][c_n]
	#print dico_fasta["GCA_002335675.1_C.eustigma_genome_v1.0_protein.faa"]["GAX86705.1"]
	#print dico_fasta["GCA_002891735.1_TetSoc1_protein.faa"]["PNG51489.1"]
	return dico_fasta


def dico_G_regions(G_regions_file) :
	"""
	Fonction permettant de retourner un dictionnaire
	input : 
		- G_regions_file : fichier txt contenant les regions genomiques
		> nom fichier de l'espece
		contig sstart send
	output :
		- retourne un dico
	"""
	dico_G = {}
	filin = open(G_regions_file, 'r')
	for line in filin :
		#print line[:-1]
		if '>' in line :
			species = line[1:-1].replace(" ","")
			dico_G[species] = {}
		else : 
			#print line.split('\t')[0], line.split('\t')[1], line.split('\t')[2]
			region_name = line.split('\t')[0]
			dico_G[species][region_name] = (int(line.split('\t')[1]), int(line.split('\t')[2][:-1]))
	filin.close()
	#verif dico
	#pour chaque fichier / espece
	for s in dico_G :
		print s, 
		print len(dico_G[s])
		#pour chaque region
		#for r in dico_G[s] :
			#print r, dico_G[s][r]
	return dico_G


def dico_TandP_regions(TorP_regions_file) :
	"""
	Fonction permettant de retourner un dictionnaire
	input : 
		- TorP_regions_file : fichier txt contenant les regions transcriptomiques ou proteomiques
		> nom fichier de l'espece
		contig 
	output :
		- retourne un dico
	"""
	dico_TorP = {}
	filin = open(TorP_regions_file, 'r')
	for line in filin :
		#print line[:-1]
		if '>' in line :
			species = line[1:-1].replace(" ","")
			dico_TorP[species] = []
		else : 
			dico_TorP[species].append(line[:-1])
	#verif dico
	#pour chaque fichier / espece
	#for s in dico_TorP :
	#	print s
	#	#pour chaque element de la liste
	#	for r in dico_TorP[s] :
	#		#l'afficher
	#		print r
	return dico_TorP


def blast_part2(cmd_blast, type_ome, regions_file, rep_fasta_files, rep_db, ecut = "1e-09") :
	"""
	Fonction permettant de lancer cmd blast
	input : 
		- cmd_blast : type de cmd blast à lancer (tblastx ou blastx)
		- type_ome : G ou T ou P
		- regions_file : fichier des regions
		- rep_fasta_files : rep contenant l'ens des fichiers fasta des diff esp
		- rep_db : rep de la banque de données (genomes, transcriptomes ou proteomes)
		- ecut : seuil pour la evalue (par defaut 1e-09 pr G et T)
	output :
		- fichier de resultat csv
	"""
	#output des fichiers blast
	outopt=" -outfmt \"6 qseqid qseq qstart qend sseqid sseq sstart send sstrand sframe length pident nident ppos gaps gapopen evalue bitscore qcovs qcovhsp qcovus\" "
	#omes des diff especes
	dbase = glob.glob(rep_db + "/*.faa") + glob.glob(rep_db + "/*.fsa_nt") + \
	glob.glob(rep_db + "/*.fasta") + glob.glob(rep_db + "/*.nt.fa") + \
	glob.glob(rep_db + "/*.fna")
	#creation du fichier fasta
	#dico contenant les infos fasta
	fastas = dico_info_fasta(rep_fasta_files)
	#pour G
	if type_ome == "G" :
		regions = dico_G_regions(regions_file)
		for species in regions :
			#print species
			for contig in regions[species] :
				filout = open(contig + ".fasta", "w")
				#print ">" + contig
				sstart = int(regions[species][contig][0])
				send = int(regions[species][contig][1])
				#print fastas[species][contig][sstart:send]
				filout.write(">" + contig + "\n")
				filout.write(fastas[species][contig][sstart:send] + "\n")
				filout.close()
				#lancement de la ligne de commande "bash" */* blast
				quer = contig + ".fasta"
				for ome in dbase :
					print quer, ome
					outfile = contig + "vs" + str(ome.split('/')[-1]) + "." + cmd_blast + ".csv"
					cmd = cmd_blast + " -query " + quer + " -db " + ome + " " + outopt + " -evalue " + ecut + " -out " + outfile + " -num_threads 32" 
					os.system(cmd)
				os.system("rm " + quer)
	#pour T ou P
	if type_ome == "T" or type_ome == "P" :
		if cmd_blast == "blastp" :
			ecut = "5"
		regions = dico_TandP_regions(regions_file)
		for species in regions :
			#print species
			for contig in regions[species] :
				filout = open(contig.replace("|","_") + ".fasta", "w")
				#print ">" + contig
				#print fastas[species][contig]
				filout.write(">" + contig + "\n")
				filout.write(fastas[species][contig] + "\n")
				filout.close()
				#lancement de la ligne de commande "bash" */* blast
				quer = contig.replace("|","_") + ".fasta"
				for ome in dbase :
					#print quer, ome
					outfile = contig.replace("|","_") + "vs" + str(ome.split('/')[-1]) + "." + cmd_blast + ".csv"
					#cmd = cmd_blast + " -query " + quer + " -db " + ome + " " + outopt + " -out " + outfile + " -num_threads 96"
					cmd = cmd_blast + " -query " + quer + " -db " + ome + " " + outopt + " -evalue " + ecut + " -out " + outfile + " -num_threads 32" 
					os.system(cmd)
				os.system("rm " + quer)


def genomic_regions_split(genomic_regions_file) :
	"""
	Fonction permettant de creer des fichiers de 21 regions max à partir du fichier de départ contenant toutes les regions genomiques
	input :
		- genomic_regions_file : fichier txt contenant toutes les regions genomiques issues du blast de la premiere partie de l'analyse
	output :
		_ fichiers txt 
	"""
	filin = open(genomic_regions_file, 'r')
	#dico contenant des listes de 21 regions par espece
	dico = {}
	for line in filin :
		if ">" in line :
			espece = line[:-1]
			dico[espece] = {}
			n = 0
			dico[espece][n] = []
		else : 
			if len(dico[espece][n]) < 21 :
				dico[espece][n].append(line[:-1])
			else :
				n += 1
				dico[espece][n] = []
				dico[espece][n].append(line[:-1])
	filin.close()
	#pour chaque espece ou couple d'especes, mettre ensemble les listes non pleines
	l = [('>GCA_002891735.1_TetSoc1_genomic.fna', '>GCA_000147415.1_v_1.0_genomic.fna'),
	('>GCA_002317545.1_ASM231754v1_Squadricauda_genomic.fna','>GCA_001662365.1_Cap_assembly01_genomic.fna'),
	('>GCA_001584585.1_ASM158458v1_genomic.fna','>MPUG01.1.fsa_nt'),
	('>GCA_002588565.1_ASM258856v1_genomic.fna','>GCA_002284615.1_Dsalina_v1.0_genomic.fna')]
	N = 0
	for esp in l :
		N += 1
		filout = open("genomic_region_" + str(N) + ".txt", "w")
		l0 = dico[esp[0]][len(dico[esp[0]]) - 1]
		l1 = dico[esp[1]][len(dico[esp[1]]) - 1]
		#print len(l0) + len(l1)
		#print l0, l1
		filout.write(esp[0] + "\n")
		#print esp[0]
		for region0 in l0 :
			filout.write(region0 + "\n")
			#print region0
		filout.write(esp[1] + "\n")
		#print esp[1] 
		for region1 in l1 :
			filout.write(region1 + "\n")
			#print region1
		filout.close()
		#supprimer ces listes du dico
		del dico[esp[0]][len(dico[esp[0]]) - 1]
		del dico[esp[1]][len(dico[esp[1]]) - 1]
	#pour les listes manquantes du dico (pleines), creer les fichiers
	
	for e in dico :
		#print e
		#print len(dico[e])
		for num in dico[e] :
			#print num
			N += 1
			filout = open("genomic_region_" + str(N) + ".txt", "w")
			filout.write(e + "\n")
			#print e
			for region in dico[e][num] :
				filout.write(region + "\n")
				#print region 


#Programme principal
if __name__ == "__main__":
	#regions imp
	rep_resume1_G = "/data/mortaza/Lucifer/evolution/oprvsgenome/results/genome/regions/resume/part1"
	rep_resume1_T = "/data/mortaza/Lucifer/evolution/oprvsgenome/results/transcriptome/regions/resume/part1"
	rep_resume1_P = "/data/mortaza/Lucifer/evolution/oprvsgenome/results/proteome/regions/resume/part1"
	#G_contig(rep_resume1_G) # >> genomic_regions.txt
	#PandT_contig(rep_resume1_T) # >> transcriptomic_regions.txt
	#PandT_contig(rep_resume1_P) # >> proteomic_regions.txt
	#infos fichiers fasta
	#dico_info_fasta(sys.argv[1])
	#dico_G_regions(sys.argv[1])
	#dico_TandP_regions(sys.argv[1])
	blast_part2(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
	#genomic_regions_split(sys.argv[1])


#Notes & Commentaires
"""
- /home/mortaza/Documents/src/OPRvsgenome/db_building_part2.py

"""