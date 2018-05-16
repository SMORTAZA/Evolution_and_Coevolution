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

#Fonctions

def dico_regions_genomiques_non_chevauchantes(genome2_csv_file) :
	"""
	Fonction permettant de faire correspondre 1 seq à 1 frame. 
	Cette fct a egalement servie à retrouver les regions à considérer 
	pour l'outil exonerate (avec quelques lignes à (de)-commenter) et 
	en utilisant la commande uniq de bash pr enlever les doublons
	input : 
		- genome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output :
		- retourne un dico
		- affiche les resultats voulus pr les enregistrer dans un 
		fichier ensuite.
	"""
	#recup du dico de la fct dico_hit_T()
	dico_hit = dico_hit_T(genome2_csv_file)
	#nouveau dico pour mettre les regions non chevauchantes de chaque frame
	new_dico = {}
	#pr chaque esp
	for i in dico_hit :
		new_dico[i] = {}
		#pr chaque contig du genome
		for j in dico_hit[i] :
			##pr resultats transcriptomiques + G à traiter dans exonerate
			print i, j, 
			##initialisation des min et max pr resultats génomiques à traiter 
			#ds exonerate.
			mini = "X"
			maxi = "X"
			new_dico[i][j] = []
			#pr chaque frame
			for k in dico_hit[i][j] :
				#print "==> frame", k
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
				new_dico[i][j].append((k,new_frames))
				##pr resultats génomiques à traiter dans exonerate
				for (smin, smax, sseq) in new_frames :
					#si pas encore de valeurs numeriques
					if mini == "X" and maxi == "X" :
						mini = smin
						maxi = smax
					#sinon comparaison
					elif mini > smin :
						mini = smin
					elif maxi < smax :
						maxi = smax
			print mini, maxi
	return new_dico

def region_sans_codon_stop(genome2_csv_file) :
	"""
	Fonction permettant de retourner pour chaque transcrit la "meilleure region"
	proteique sans codon stop pour une sequence nucleotidique.  
	input :
		- genome2_csv_file : fichier csv obtenu avec les fct 
		dico_blast_results_pt() et dico_g_query_result_analysis().
		header non indiqué ds le fichier :
		ome	contig	sstart	send	strand	Evalue	query	sseq
	output :
		- 
	"""
	#recup du dico {esp : {contig : [(sframe, (sstart, send, sseq))]}}
	dico_frames = dico_regions_genomiques_non_chevauchantes(genome2_csv_file)
	#pr chaque esp
	for sp in dico_frames :
		#pr chaque transcrit
		for ctg in dico_frames[sp] :
			#initialisation de liste pour les regions sur brin plus ou 
			#brin moins
			a = 0; b = 0
			strand_plus = []
			strand_moins = []
			#info = [(sframe, (sstart, send, sseq))]
			for (frame, info) in dico_frames[sp][ctg] :
				#print sp, ctg, frame		
				#afficher chaque elmt de la liste
				for region in info :
					#si pas de codon stop dans la region, l'enregistrer
					if "*" not in region[-1] :
						if int(frame) > 0 :
							strand_plus.append(region)
							a += 1
						else :
							strand_moins.append(region)
							b += 1
			print "*"
			for i in sorted(strand_plus, key=lambda data: data[0]) :
				print sp, ctg, i
			print "*"
			for i in sorted(strand_moins, key=lambda data: data[0]) :
				print sp, ctg, i
			#ttt des strands 
			## si y a rien ds les listes, longueur nulle pr la prot predite
			#if len(strand_plus) == 0 :
			#	seq_plus = ""
			#elif len(strand_moins) == 0 :
			#	seq_moins = ""
			##sinon, retrouver la region predite dans chacun des strands
			#seq_plus = fct(strand_plus)
			#seq_moins = fct(strand_moins)
			##comparer ensuite ces regions pour prendre la plus grande
			#if len(seq_plus) >= len(seq_moins) :
			#	seq_prot = seq_plus
			#else :
			#	seq_prot = seq_moins

def list_of_regions_to_1region(list_of_regions) :
	"""
	Fonction permettant de donner la seq de la region issue d'une liste de 
	regions.
	input :
		- list_of_regions : liste des regions [(sstart, send, sseq)]
	output :
		- retourne une séquence (chaine de car)
	"""
	#print list_of_regions
	#voir s'il y a des regions chevauchantes :
	#pr chaque region :
	smin = "X"
	smax = "X"
	for (sstart, send, sseq) in list_of_regions :
		if smin == "X" and smax == "X" :
			smin = sstart
			smax = send

	#en cours d'écriture

def list_regions_per_species_T(ome_exonerate_file, rep_fasta_files) :
	"""
	Fonction permettant de créer pour chaque espèce un fichier contenant
	les regions du blast2.
	input :
		- ome_exonerate_file : fichier obtenu avec la fct
		dico_regions_genomiques_non_chevauchantes() avec header :
		ome contig (start end)
		- rep_fasta_files : rep contenant les fichiers fasta des transcriptomes
	output :
		- crée des fichiers fasta par esp
	"""
	#dico contenant les infos fasta
	fastas = dico_info_fasta(rep_fasta_files)
	#dico contenant les infos à afficher et à enregistrer dans fichiers
	dico = {}
	filin = open(ome_exonerate_file, "r")
	for line in filin :
		species = line.split(" ")[0]
		contig = line.split(" ")[1][:-1]
		if species not in dico :
			dico[species] = []
		dico[species].append(contig)
	filin.close()
	for i in dico :
		#creation du fichier fasta par espece
		output = i + "_exonerate.fasta"
		print output
		filout = open(output, "w")
		for j in dico[i] :
			filout.write(">" + j + "\n")
			filout.write(fastas[i][j] + "\n")
			#print ">" + j
			#print fastas[i][j]
		filout.close()

def list_regions_per_species_G(ome_exonerate_file, rep_fasta_files) :
	"""
	Fonction permettant de créer pour chaque espèce un fichier contenant
	les regions du blast2.
	input :
		- ome_exonerate_file : fichier obtenu avec la fct
		dico_regions_genomiques_non_chevauchantes() avec header :
		ome contig (start end)
		- rep_fasta_files : rep contenant les fichiers fasta des transcriptomes
	output :
		- crée des fichiers fasta par esp
	"""
	#dico contenant les infos fasta
	fastas = dico_info_fasta(rep_fasta_files)
	#dico contenant les infos à afficher et à enregistrer dans fichiers
	dico = {}
	filin = open(ome_exonerate_file, "r")
	for line in filin :
		species = line.split(" ")[0]
		contig = line.split(" ")[1]
		if species != "Chlamydomonasreinhardtii.fasta" :
			region = contig
			if "|" in region :
				region = region.split("|")[-2]
			contig = region
		sstart = int(line.split(" ")[2])
		send = int(line.split(" ")[3][:-1])
		if species not in dico :
			dico[species] = []
		dico[species].append((contig, sstart, send))
	filin.close()
	for i in dico :
		#creation du fichier fasta par espece
		output = i + "_exonerate.fasta"
		print output
		filout = open(output, "w")
		#print fastas[i]
		for (j, st, se) in dico[i] :
			#print j
			filout.write(">" + j + " " + i + "\n")
			filout.write(fastas[i][j][st:se] + "\n")
			#print ">" + j, i,
			#print fastas[i][j][st:se]
		#print fastas["Chlamydomonasreinhardtii.fasta"]["gi|154749079|gb|ABCN01007600.1|"][2524:3405]
		filout.close()

def protein_2_genome_exonerate(type_ome, query_fasta, subject_fasta) :
	"""
	Lancement de la commande exonerate P2G
	"""
	#pour T ou P
	if type_ome == "T" :
		cmd = '/home/mortaza/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome '\
		 + query_fasta + ' ' + subject_fasta + \
		 ' --ryo ">%ti %td %pi\n%tcs" --showalignment FALSE > ' + "OPR_" + subject_fasta.split("/")[-1]
		os.system(cmd)

def exonerate_file_reader(exonerate_file_rep) :
	"""
	Fonction permettant de lire des fichiers obtenus avec exonerate 
	input :
		- exonerate_file_rep : rep des fichiers obtenus avec la commande 
		exonerate - modele protein2genome - et --ryo ">%ti %td %pi\n%tcs"
		et --showalignment FALSE
	output :
		- affiche les infos voulues et retourne un dico 1 ctg = 1 seq
	"""
	#creation dico contenant les infos du fichier
	dico = {}
	#ens des fichiers resultats de exonerate
	files = glob.glob(exonerate_file_rep + "/*")
	#pr chaque fichier
	for file in files :
		#print file
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
				dico[ome][subject].append((sstart, send, strand, query, qstart, qend, score, p_identity, seq))
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
				score = line.split(" ")[9]
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
	#lecture du dico et enregistrement ds un autre dico de l'info 1 contig = 1 seq
	dico_1ctg_1seq = {}
	#pr chaque esp 
	for sp in dico :
		#initialisation nvo dico avec comme cle le nom de l'esp
		dico_1ctg_1seq[sp] = {}
		#pr chaque subject
		for ctg in dico[sp] :
			#print ctg
			maxi_pi = 0
			maxi_score = 0
			maxi_sstart = "X"
			maxi_send = "X"
			maxi_strand = "X"
			maxi_query = "X"
			maxi_qstart = "X"
			maxi_qend = "X"
			maxi_seq = "X"
			#pr chaque query
			for i in dico[sp][ctg] :
				#print i
				#print i
				#selection du hit ayant le meilleur score 
				if float(i[-2]) > maxi_score :
					maxi_score = float(i[-3])
					maxi_pi = float(i[-2])
					maxi_sstart = int(i[0])
					maxi_send = int(i[1])
					maxi_strand = i[2]
					maxi_query = i[3]
					maxi_qstart = int(i[4])
					maxi_qend = int(i[5])
					maxi_seq = i[-1]
			#print ctg, maxi_sstart, maxi_send, maxi_strand, maxi_query, \
			#maxi_qstart, maxi_qend, maxi_score, maxi_pi, maxi_seq
			#ajout contig et ses infos si pas dans le dico
			dico_1ctg_1seq[sp][ctg] = (maxi_sstart, maxi_send, maxi_strand,\
			maxi_query, maxi_qstart, maxi_qend, maxi_score, maxi_pi, maxi_seq)
	return dico_1ctg_1seq

def exonerate_multifasta_prot(exonerate_file_rep) :
	"""
	Fonction permettant de donner un affichage type multifasta avec les cds
	traduits.
	input : 
		- exonerate_file_rep :  rep des fichiers obtenus avec la commande 
		exonerate - modele protein2genome - et --ryo ">%ti %td %pi\n%tcs"
		et --showalignment FALSE
	output :
		- 
	"""
	#recup des infos
	dico = exonerate_file_reader(exonerate_file_rep)
	#pr chaque esp
	for sp in dico :
		#pr chaque contig
		for ctg in dico[sp] :
			#print sp, ctg,
			info = dico[sp][ctg]
			strand = info[2]
			#print strand
			region_name = sp[4:-10] + "_" + ctg + "_" + str(info[0]) + "_" + str(info[1])
			print ">" + region_name
			seq = info[-1]
			#print seq
			print translate_seq(seq)

def translate_seq(seq, strand = "+") :
	"""
	Fonction permettant de traduire une sequence nucléotidique.
	input :
		- seq : chaine de car. seq nucleotidique.
		- strand : "+" ou "-" . par defaut, "+".
	output : 
		- retourne la seq protéique
	"""
	#dico des codons - aa
	code_genetique = {

  'TCA':'S', # Sérine

  'TCC':'S', # Sérine

  'TCG':'S', # Sérine

  'TCT':'S', # Sérine

  'TTC':'F', # Phénylalanine

  'TTT':'F', # Phénylalanine

  'TTA':'L', # Leucine

  'TTG':'L', # Leucine

  'TAC':'Y', # Tyrosine

  'TAT':'Y', # Tyrosine

  'TAA':'*', # Codon Stop

  'TAG':'*', # Codon Stop

  'TGC':'C', # Cystéine

  'TGT':'C', # Cystéine

  'TGA':'*', # Codon Stop

  'TGG':'W', # Tryptophane

  'CTA':'L', # Leucine

  'CTC':'L', # Leucine

  'CTG':'L', # Leucine

  'CTT':'L', # Leucine

  'CCA':'P', # Proline

  'CCC':'P', # Proline

  'CCG':'P', # Proline

  'CCT':'P', # Proline

  'CAC':'H', # Histidine

  'CAT':'H', # Histidine

  'CAA':'Q', # Glutamine

  'CAG':'Q', # Glutamine

  'CGA':'R', # Arginine

  'CGC':'R', # Arginine

  'CGG':'R', # Arginine

  'CGT':'R', # Arginine

  'ATA':'I', # Isoleucine

  'ATC':'I', # Isoleucine

  'ATT':'I', # Isoleucine

  'ATG':'M', # Méthionine

  'ACA':'T', # Thréonine

  'ACC':'T', # Thréonine

  'ACG':'T', # Thréonine

  'ACT':'T', # Thréonine

  'AAC':'N', # Asparagine

  'AAT':'N', # Asparagine

  'AAA':'K', # Lysine

  'AAG':'K', # Lysine

  'AGC':'S', # Sérine

  'AGT':'S', # Sérine

  'AGA':'R', # Arginine

  'AGG':'R', # Arginine

  'GTA':'V', # Valine

  'GTC':'V', # Valine

  'GTG':'V', # Valine

  'GTT':'V', # Valine

  'GCA':'A', # Alanine

  'GCC':'A', # Alanine

  'GCG':'A', # Alanine

  'GCT':'A', # Alanine

  'GAC':'D', # Acide Aspartique

  'GAT':'D', # Acide Aspartique

  'GAA':'E', # Acide Glutamique

  'GAG':'E', # Acide Glutamique

  'GGA':'G', # Glycine

  'GGC':'G', # Glycine

  'GGG':'G', # Glycine

  'GGT':'G', # Glycine 

  }
	#seq proteique initialisation
	protein_seq = ""
	if strand == "+" :
		#lecture de la seq et traduction au fur et a mesure
		for i in range(0, len(seq), 3) :
			if seq[i:i+3].upper() in code_genetique :
				protein_seq += code_genetique[seq[i:i+3].upper()]
			else :
				protein_seq += "X"
	elif strand == "-" :
		#print seq
		#print "*"
		c_seq = complementaire_seq(seq)
		#print c_seq
		#lire la seq en sens inverse (-1 pr Frame 6, -2 pr Frame 5)
		for i in range(len(c_seq) - 1, -1, -3) :
			r_seq = ""
			for j in reversed(c_seq[i:i+3].upper()) :
				r_seq += j
			if len(r_seq) == 3 :
				#print r_seq, code_genetique[r_seq]
				if r_seq in code_genetique :
					protein_seq += code_genetique[r_seq]
				else :
					protein_seq += "X"
	#print c_seq 
	return protein_seq
	#print len(seq) 
	#print len(protein_seq) * 3

def complementaire_seq(seq) :
	"""
	Fonction permettant de donner la sequence complementaire d'une sequence
	input :
		- seq : une chaine de car correspondant à la seq
	output :
		- retourne une chaine de car, le complementaire de la seq donnée
	"""
	#dico indiquant le complementaire de chaque base
	dico = {"A" : "T", "T" : "A", "U" : "A", "C" : "G", "G" : "C", \
	"N" : "N", "S" : "N", "Y" : "N"}
	complementaire_seq = ""
	for nt in seq :
		complementaire_seq += dico[nt]
	return complementaire_seq




#Programme principal
if __name__ == "__main__":
	#analyse régions génomiques
	csv_files_g_t = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryT/SubjectG/*')
	#dico_blast_results_pt(sys.argv[1], csv_files_g_t)
	rep_resume1_files = "/data/mortaza/Lucifer/evolution/oprvsgenome/results/genome/regions/resume/part1"
	csv_files_g_g = glob.glob('/data/mortaza/Lucifer/evolution/hitvsome/QueryG/SubjectG/*')
	#dico_g_query_result_analysis(sys.argv[1], csv_files_g_g, rep_resume1_files, "g") 
	#region_sans_codon_stop(sys.argv[1])
	l = [(1,	69,	'RKAPTRALVPAAVRGASALPLSG'),
	(24,	74,	'RPCCCPRCERSSPLWPS'),
	(147,	212,	'CLTLSATSQARQTGFPPGSASP'),
	(163,	249,	'PHRKLGRQVSRRAVLLQRGVGAVELACGP'),
	(304,	474,	'DGAGRAGSALPGAGGPTGVAAACVVQSVPRGVSPAASQRTPRALGSHDLVSQPLSPC'),
	(1025,	1171,	'SPTHYPCRRERVWSGGIDAAMHSMLVQLWCMRTAASYKLQRRSLVVPGD')]
	#list_of_regions_to_1region(l)
	#dico_regions_genomiques_non_chevauchantes(sys.argv[1])
	#list_regions_per_species_G(sys.argv[1], sys.argv[2])
	#protein_2_genome_exonerate("T", sys.argv[1], sys.argv[2])
	# /data/mortaza/Data/OPR/Chlre.fasta
	exonerate_multifasta_prot(sys.argv[1])
	#seq = "GCGTCGGCGGCGAACAGTCGGCATGCGCCGACACGACCGTCGCGTGCGCCGCCGTCTGGTTGACCTGCTCCACCACCAGCTCCCTGAACCCCGCATCTGCGTCCTGCTGCATCACCACCGGCACCAGCGTTTCCATCGAGTCCGGGCCTGCCAGCAGCCCCGCCTCGCCCGGCCGCGCCTCTGCCGCCAGCACCGCAACGCCGTCGCTCTGCGTCAGCTCGGCAAACTCGAGCGAGAGGCTCACGCTGCGTCGCCGCAGGCACAGCCGCGGCGCCACCGTCCGCACGTCCGATAGCAGCGAGCGCAGCG" 
	#ranslate_seq(seq, "-")
 
"""
Notes & Commentaires :
- /data/mortaza/src/OPRvsgenome/result_G2_analysis.py
- liste de notre db perso (g, t, p) ds le rep
/data/mortaza/Lucifer/evolution/OPR_evolution/oprvsgenome/results/list
- chemin pour avoir les rep des tailles des diff contigs
/data/mortaza/Data
- pr utiliser l'executable exonerate :
/home/mortaza/exonerate-2.2.0-x86_64/bin/exonerate


/data/mortaza/src/OPRvsgenome/result_G2_analysis.py /data/mortaza/Data/OPR/Chlre.fasta ../
"""



