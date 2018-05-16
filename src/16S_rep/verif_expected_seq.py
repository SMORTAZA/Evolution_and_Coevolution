#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author : Shogofa MORTAZA
IBPC Internship
"""

#chemin de ce code 
#/data/mortaza/src/16S_rep/verif_expected_seq.py

#import des modules
import os 
import sys
import glob

#fonctions

#fichiers fasta genome + transcriptome
#genome_list = glob.glob("/workdir/ibpc_team/lafontaine_team/shared/databases/BlastDB/genomes/*.fna") + glob.glob("/workdir/ibpc_team/lafontaine_team/shared/databases/BlastDB/genomes/*.fsa_nt")
#transcriptome_list = glob.glob("/workdir/ibpc_team/lafontaine_team/shared/databases/BlastDB/transcriptomes/*.fa")

#resultats de blastn
#genome_S16_seq = glob.glob("/workdir/ibpc_team/lafontaine_team/mortaza/ingrid2shogofa/rna16s/genomes_16S/*.blastn.csv")
#transcriptome_S16_seq = glob.glob("/workdir/ibpc_team/lafontaine_team/mortaza/ingrid2shogofa/rna16s/transcriptomes_16S/*.blastn.csv")

#exemple
#GCA_000143455.1_v1.0_genomic.fna ==> genome de l'espece
#G_16s_vs_GCA_000143455.1_v1.0_genomic.fna.blastn.csv ==> resultats blast de l'espece etudiée
#G_16S_GCA_000143455.1_v1.0_genomic.fna.blastn.csv ==> qques resultats imp sur la seq de l'esp selectionnés à partir de blast (1 ligne = 1 resultat)

def s_from_blast(fichier_blast) :
	"""
	Fonction permettant de prendre les col sseqid, sstart, send et sseq du fichier blast
	input :
		- fichier_blast : fichier de resultat blast avec toutes les données en colonnes
	output : 
		- fichier de sortie contenant les infos imp de blast
		- retourne le nom du fichier de sortie
	"""
	#nom du fichier de sortie
	name = 'G_16S' + fichier_blast.split('/')[-1][8:] #NOTE : selon nom debut de fichier, modif de [8:]
	#print name
	#elimine le fichier s'il existe deja
	cmd_rm = "rm " + name
	os.system(cmd_rm)
	# commande awk
	cmd = "awk '{print $5,$7,$8,$6}' " + fichier_blast + " >> " + name
	#print cmd
	os.system(cmd)
	return name

def dico_blast_16S(blast_16S, type_ome) :
	"""
	Fonction permettant de mettre dans un dictionnaire toutes les infos du fichier blast_16S
	input :
		- blast_16S : fichier contenant l'id du contig, le start, le end et la seq du 16s donné par blast
		- type_ome : G pour genome ou T pour transcriptome
	output :
		- retourne un dico[id du contig] = (start, end, seq)
	"""
	#dico qui contiendra les infos de la seq
	seq = {}
	#lecture des lignes du fichier
	filin = open(blast_16S, 'r')
	lines = filin.readlines()
	filin.close()
	#pour chaque ligne du fichier, mettre dans le dico les infos sous forme de tuples
	for l in lines :
		s_id, sstart, send, sseq = l.split()
		if type_ome == 'G' :
			s_id = s_id.split('|')[1]
		if s_id not in seq :
			seq[s_id] = (sstart, send, sseq.replace('-',''))
		else : 
			print 'Error, cet id a déjà été utilisé ...'
			break
	return seq

def comparison_with_ome_16S(fichier_ome, s_id, sstart, send, sseq) :
	"""
	Fonction permettant de trouver la seq correspondant à l'id et aux start et end
	input :
		- fichier_ome : fichier fasta contenant le genome ou le transcriptome
		- s_id : identifiant du contig
		- sstart : debut de la seq voulue
		- send : fin de la seq voulue

	output :
		- retourne une sequence sous forme de chaine de car
	"""
	filin = open(fichier_ome, 'r')
	for l in filin :
		if l[0] == ">":
			flag = 0
			#print l[:-1]
			id_contig = l.split()[0][1:]
			#print id_contig
			if s_id == id_contig :
				#print l[:-1]
				flag = 1
				seq = ''
		if l[0] != ">" and flag == 1 :
			#print l[:-1]
			seq = seq + l[:-1]
	if int(sstart) < int(send) :
		if seq[int(sstart)-1:int(send)] == sseq :
			print 'ok'
	else :
		print seq[int(send):int(sstart)+1]
		print 'pas ok'

def seq16S_expected(fichier_ome, type_ome, fichier_blast) :
	"""
	Fonction permettant de verifier pour un ome si bonne seq 16S donnée ou pas avec blast
	input : 
		- fichier_ome : fichier fasta contenant le genome ou le transcriptome
		- type_ome : G pour genome ou T pour transcriptome
		- fichier_blast : fichier de resultat blast avec toutes les données en colonnes
	output :
		- message donnant l'id du contig + 'ok' ou 'pas ok'
	"""

	#1°) a partir du fichier blast, ne prendre que les colonnnes les plus importantes
	blast_16S = s_from_blast(fichier_blast)
	#print blast_16S

	#2°) Mettre dans un dico les infos du fichier blast_16S : cle = id de la seq et valeur = (start, end, seq)
	dico_seq_16S = dico_blast_16S(blast_16S, type_ome)

	#3°) Comparaison avec la seq du genome
	for i in dico_seq_16S:
		print i 
		print dico_seq_16S[i][0], dico_seq_16S[i][1]
		print dico_seq_16S[i][2]
		comparison_with_ome_16S(fichier_ome, i, dico_seq_16S[i][0], dico_seq_16S[i][1], dico_seq_16S[i][2])



#programme principal
if __name__ == "__main__":
	#s_from_blast('/home/mortaza/Documents/16S_PhylogenyTree/G_16s_vs_GCA_000143455.1_v1.0_genomic.fna.blastn.csv')
	#dico_blast_16S('/home/mortaza/Documents/16S_PhylogenyTree/G_16S_GCA_000143455.1_v1.0_genomic.fna.blastn.csv', 'G')
	#sequence = "TCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTCGTAAAGTCTAATGTCAAATACCAGGGCTCAACCTTGGATCGGCATTAGAGCACTCACGAGCTTGAGTACGGTAGAGGCAGAGGGAACTCCAAGTGGAGCGGTGAAATGCGTAGAGATTTGGAAGAACACCAGAGGCGAAGGCGCTCTGCTGGGCCGAAACTGACACTGAGAGACGAAAGCTGGGGGAGCGAATAGGATTAGATACCCTAGTAGTCCCAGCCGTAAACTATGGAGACTAAGTGCTGCCGCAAGCAGTGCTGTAGCTAACGCGTTAAGTCTCCCGCCTGGGGAGTATGCTCGCAAGAGTAAACTCAAAGGAATTGACGGGGACCCGCACAAGAGGTGGATTATGTGGATTAATTCGATGCAACGCGAAGAACCTTACCAGGGTTTGACATGTCAAGAACTTCCCAGAAATGGGAGGGTGCCCTAACGG"
	#comparison_with_ome_16S('/home/mortaza/Documents/16S_PhylogenyTree/GCA_000143455.1_v1.0_genomic.fna', 'GL378338.1', '937889', '938362', sequence)
	seq16S_expected('/home/mortaza/Documents/16S_PhylogenyTree/GCA_000143455.1_v1.0_genomic.fna', 'G', '/home/mortaza/Documents/16S_PhylogenyTree/G_16s_vs_GCA_000143455.1_v1.0_genomic.fna.blastn.csv')



"""
#Aide à l'ecriture du code

liste contenant l'ens des fichiers de resultats blast 16S à partir des genomes.
liste contenant l'ens des genomes

GL378338.1
937889 938362
TCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTCGTAAAGTCTAATGTCAAATACCAGGGCTCAACCTTGGATCGGCATTAGAGCACTCACGAGCTTGAGTACGGTAGAGGCAGAGGGAACTCCAAGTGGAGCGGTGAAATGCGTAGAGATTTGGAAGAACACCAGAGGCGAAGGCGCTCTGCTGGGCCGAAACTGACACTGAGAGACGAAAGCTGGGGGAGCGAATAGGATTAGATACCCTAGTAGTCCCAGCCGTAAACTATGGAGACTAAGTGCTGCCGCAAGCAGTGCTGTAGCTAACGCGTTAAGTCTCCCGCCTGGGGAGTATGCTCGCAAGAGTAAACTCAAAGGAATTGACGGGGACCCGCACAAGAGGTGGATTATGTGGATTAATTCGATGCAACGCGAAGAACCTTACCAGGGTTTGACATGTCAAGAACTTCCCAGAAATGGGAGGGTGCCCTAACGG
GL378454.1
53219 53482
GL378376.1
509603 509571
GL378351.1
181928 182326
"""
