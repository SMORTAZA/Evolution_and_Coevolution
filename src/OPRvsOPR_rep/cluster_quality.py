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
sys.path.append('/home/mortaza/Documents/src/16S_rep/')
from blast_analysis import *
sys.path.append('/home/mortaza/Documents/src/OPRvsOPR_rep/')
from fonctions import *

#Fonctions

def true_cluster(output_MCL) :
	"""
	Fonction permettant de retourner les vrais clusters (au moins deux composants)
	input : 
		- output_MCL : fichier de sortie MCL avec une certaine valeur de I
	output :
		- retourne un dico avec comme clé le num du cluster et comme valeur la liste des proteines OPR
	"""
	dico_cluster = {}
	#initialisation du numero du cluster
	num_cluster = 0
	filin = open(output_MCL, 'r')
	#pour chaque ligne du fichier
	for l in filin :
		num_cluster += 1
		#si le nombre de composants du cluster est different de 1
		if l.count('\t') != 0 :
			if num_cluster not in dico_cluster :
				dico_cluster[num_cluster] = []
			for sep in range(l.count('\t') + 1) :
				#ajout du nom de la proteine dans la liste du cluster
				#print l.split()[sep]
				dico_cluster[num_cluster].append(l.split()[sep])
	filin.close()
	#verif du contenu du dico
	#for i in dico_cluster :
		#print i, dico_cluster[i]
	return dico_cluster

def make_fasta_clusters(output_MCL, fasta_rep_name) :
	"""
	Fonction permettant de creer un fichier fasta par cluster
	input : 
		- output_MCL : fichier de sortie MCL avec une certaine valeur de I
		- fasta_rep_name : nom du rep contenant l'ensemble des fichiers fasta
	output :
		- fichier fasta par cluster
	"""
	#recup des infos cluster-prot avec la fonction true_cluster
	dico_cluster = true_cluster(output_MCL)
	#division du cluster 1 en deux 
	#dico_cluster = {'NCL_1' : ['Chlre_OPR75', 'Chlre_OPR72', 'Chlre_OPR73', 'Chlre_OPR70', \
	#'Chlre_OPR71', 'Chlre_OPR76', 'Chlre_OPR77', 'Chlre_OPR79', 'Chlre_OPR94', 'Chlre_OPR74', \
	#'Chlre_OPR82', 'Chlre_OPR106', 'Chlre_OPR92', 'Chlre_OPR93', 'Chlre_OPR90', 'Chlre_OPR91', \
	#'Chlre_OPR96', 'Chlre_OPR97', 'Chlre_OPR95', 'Chlre_OPR98', 'Chlre_OPR99', 'Chlre_OPR87', \
	#'Chlre_OPR84', 'Chlre_OPR80', 'Chlre_OPR20', 'Chlre_OPR112', 'Chlre_OPR110', 'Chlre_OPR111', \
	#'Chlre_OPR88', 'Chlre_OPR83', 'Chlre_OPR81', 'Chlre_OPR86', 'Chlre_OPR69', 'Chlre_OPR78', 'Chlre_OPR21'], 
	#'autre_1' : ['Chlre_OPR107', 'Chlre_OPR59', 'Chlre_OPR44', 'Chlre_OPR109', 'Chlre_OPR6', \
	#'Chlre_OPR54', 'Chlre_OPR4', 'Chlre_OPR46', 'Chlre_OPR5', 'Chlre_OPR17', 'Chlre_OPR52', \
	#'Chlre_OPR57', 'Chlre_OPR58', 'Chlre_OPR39', 'Chlre_OPR68']}
	#pour chaque cluster
	for cluster in dico_cluster :
		#print cluster
		#creation fastafile pour cluster en question
		filout = open("cluster_" + str(cluster) + ".fasta", "w")
		#pour chaque prot
		for prot in dico_cluster[cluster] :
			#print prot
			#print fasta_rep_name + "/" + prot + ".fasta"
			#lecture du fichier fasta correspondant à cette prot
			fasta_file = open(fasta_rep_name + "/" + prot + ".fasta", "r")
			for line in fasta_file :
				#print line[:-1]
				filout.write(line)
			fasta_file.close()
		filout.close()

def MAFFT_alignment(fasta_rep) :
	"""
	Fonction permettant de réaliser un alignement multiple MAFFT de fichiers fasta (avec des param par defaut)
	input :
		- fasta_rep : rep contenant les fichiers multifasta donnés en entrée
	output :
		- fichiers d'alignements multiple de ces multifasta au format fasta
	"""
	fasta_rep_files = glob.glob(fasta_rep + "/*.fasta")
	for fasta_file in fasta_rep_files :
		#commande bash à lancer
		cmd = "/home/mortaza/mafft-linux64/mafft.bat " + str(fasta_file) + " >> MAFFT_" + str(fasta_file.split('/')[-1])
		os.system(cmd)

def seq_length(fasta_file_rep, output = 'fasta_len.txt') :
	"""
	Fonction permettant de donner la longueur de la séquence (utile pour avoir le qlen)
	input :
		- fasta_file_rep : rep des fichiers fasta pour lesquels on veut avoir la longueur de la séq
	output :
		- un fichier texte résumant tout les resultats
	"""
	filout = open(output, 'w')
	fasta_file_list = glob.glob(fasta_file_rep + '/*f*')
	for fasta_file in fasta_file_list : 
		#identifiant de la sequence
		qid = fasta_file.split('/')[-1].split('.')[0]
		file = open(fasta_file, 'r')
		#ne pas lire le header du fasta (>)
		line = file.readline()
		#lecture de la sequence sur une ou plusieurs lignes
		qlen = 0
		for line in file :
			qlen = qlen + len(line[:-1])
		#print qid, qlen
		filout.write(str(qid) + '\t' + str(qlen) + '\n')
		file.close()
	filout.close()

def general_add_qcov_blast(blast_file, query_length) : 
	"""
	Fonction permettant d'ajouter la valeur du qcov dans le fichier de resultat blast.
	input : 
		- blast_file : fichier resultat de blastn (outfmt 6 - tabular)
		- query_length : fichier texte contenant en deux colonnes, le qid et le qlen
	output : 
		- ce meme fichier avec une colonne en plus du qcov
	"""
	#dico contenant qlen des diff qid
	qlen = {}
	qlen_file = open(query_length, 'r')
	#pour chaque ligne du fichier 
	for line in qlen_file :
		#si qid pas dans le dico qlen
		if line.split()[0] not in qlen :
			if line.split()[0] == "Chlre_OPR67=Cre14" :
				qlen[line.split()[0] + ".g633150.t1.3"] = line.split()[1]
			else :
				#l'ajouter avec la valeur de son qlen
				qlen[line.split()[0]] = line.split()[1]
	qlen_file.close()
	#verif de la composition du dico
	#for i in qlen :
		#print i, qlen[i]
	#nom du fichier de sortie
	output_file = blast_file + ".qcov"
	#print output_file
	#ouverture des fichiers
	filin = open(blast_file, 'r')
	filout = open(output_file, 'w')
	#pour ne pas lire le header
	l = filin.readline()
	#print l[:-1]+" qcov"
	#ajout de header dans le nouveau fichier
	#filout.write(l[:-1]+" qcov\n")
	#pour chaque ligne, calcul du qcov et ajout de la valeur dans le fichier 
	for l in filin :
		#print l.split()[2], l.split()[3], qlen[l.split()[0]]
		#print query_coverage(l.split()[2], l.split()[3], qlen[l.split()[0]])
		filout.write(l[:-1]+ " " + str(query_coverage(l.split()[2], l.split()[3])) + "\n")
	#fermeture des fichiers
	filin.close()
	filout.close()


#Programme principal
if __name__ == "__main__":
	print 
	#make_fasta_clusters(sys.argv[1], sys.argv[2])
	#MAFFT_alignment(sys.argv[1])
	#****************************
	#calcul qlen pour chaque chaque OPR
	#fasta_file_list = glob.glob(sys.argv[1] + '/*fasta')
	#seq_length(sys.argv[1])
	#Ajout qcov dans fichier resultat brut de blast
	#general_add_qcov_blast(sys.argv[1], sys.argv[2])
	#awk '{if ($22 > 30 && $17 < 0.001) print $1, $5, $17}' opr_vs_opr.blastp.csv.qcov > opr_vs_opr_selected_30_qcov_0.001_evalue.txt
	#prendre les meilleurs hits
	#convert_blast_to_cytoscape(sys.argv[1], sys.argv[2], sys.argv[3])
	#commande bash MCL avec I = 5.0
	#/home/mortaza/local/bin/mcl opr_vs_opr_selected_30_qcov_0.001_evalue.txt --abc --abc-neg-log10 -I 5.0
	#analyse type proteines avec les resultats MCL
	#mcl_analysis(sys.argv[1], sys.argv[2]) 
	#creation multifasta des clusters
	#make_fasta_clusters(sys.argv[1], sys.argv[2])
	#alignement MAFFT
	#MAFFT_alignment(sys.argv[1])

#Notes & Commentaires
"""
- Les fichiers fasta des prot OPRs se trouvent dans le rép /data/mortaza/Lucifer/oprvsgenome/OPRs_fasta/fasta_files
- ligne de commande à lancer dans le rep /data/mortaza/Lucifer/oprvsopr/genes/out1
/home/mortaza/Documents/src/OPRvsOPR_rep/cluster_quality.py out.opr_vs_opr_0.001.txt.I50 /data/mortaza/Lucifer/oprvsgenome/OPRs_fasta/fasta_files
- installation de MAFFT pour l'alignement multiple
-- telechargement en ligne de commande
wget https://mafft.cbrc.jp/alignment/software/mafft-7.222-with-extensions-src.tgz ==> a besoin de administrateur qd meme
-- autre telechargement (All-in-one package for Linux) suivre le site
wget https://mafft.cbrc.jp/alignment/software/mafft-7.380-linux.tgz ==> ok
- utilisation de MAFFT en lancant 
/home/mortaza/mafft-linux64/mafft.bat fichier.fasta (affiche a l'ecran l'alignement obtenu)
*******************
- avec qcov
- creation un fichier contenant qlen pour chaque OPR.
- selection des qcov > 30 % et evalue < 0.001
awk '{if ($22 > 30 && $17 < 0.001) print $5, $17, $22}' opr_vs_opr.blastp.csv.qcov 
awk '{if ($22 > 30 && $17 < 0.001) print $1, $5, $17}' opr_vs_opr.blastp.csv.qcov > opr_vs_opr_selected_30_qcov_0.001_evalue.txt
"""