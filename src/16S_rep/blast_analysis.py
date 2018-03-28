#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author : Shogofa MORTAZA
IBPC Internship
"""

"""
Ce fichier regroupe les différents traitements possibles après un resultat de Blastn.
Exemple avec :
/home/mortaza/Documents/16S_PhylogenyTree/G_16s_vs_GCA_000143455.1_v1.0_genomic.fna.blastn.csv
"""

#import de modules
import os
import sys
import glob

#Fonctions

def add_header_blast(blast_file) :
	"""
	Fonction permettant d'ajouter un header au fichier resultat de Blastn sans header.
	input : 
		- blast_file : fichier resultat de blastn (outfmt 6 - tabular)
	output : 
		- ce meme fichier avec un header
	"""
	header="qseqid qseq qstart qend sseqid sseq sstart send sstrand sframe length pident nident ppos gaps gapopen evalue bitscore qcovs qcovhsp qcovus"
	os.system("sed -i '1i"+header+"' "+blast_file)

def query_coverage(qstart, qend, qlen = 1474) :
	"""
	Fonction permettant de calculer le Query Coverage [qcov = abs((qend-qstart)/length(query))]
	input :
		- qstart : start du query
		- qend : end du query
		- qlen : longueur du query. par defaut, longueur de 16S de Chre.
	output :
		- retour du qcov (en pourcentage)
	"""
	return abs((float(qend) - float(qstart)) / float(qlen))*100

def add_qcov_blast(blast_file) : 
	"""
	Fonction permettant d'ajouter la valeur du qcov dans le fichier de resultat blast.
	input : 
		- blast_file : fichier resultat de blastn (outfmt 6 - tabular)
	output : 
		- ce meme fichier avec une colonne en plus du qcov
	"""
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
	filout.write(l[:-1]+" qcov\n")
	#pour chaque ligne, calcul du qcov et ajout de la valeur dans le fichier 
	for l in filin :
		#print l.split()[2], l.split()[3]
		#print query_coverage(l.split()[2], l.split()[3])
		filout.write(l[:-1]+ " " + str(query_coverage(l.split()[2], l.split()[3])) + "\n")
	#fermeture des fichiers
	filin.close()
	filout.close()

def selection_by_qcov(qcov_file_list, th_qcov = 70, output = "qcov_selection.csv") :
	"""
	Fonction permettant de selectionner les resultats de blast ayant un qcov inf au seuil_qcov
	input : 
		- qcov_file_list : liste des fichiers resultat blast + qcov (/!\ car y a des splits)
		- th_qcov : seuil pour le qcov. par defaut, à 70 %
		- output : nom donné pour le fichier 
	output :
		- un fichier avec l'ens des resultats 
	"""
	filout = open(output, 'w')
	#nb d'especes apres la selection
	nb = 0
	for file in qcov_file_list :
		#index = index + 1
		index = file.split('/')[7][9:-16]
		flag = 0
		reading_file = open(file, 'r')
		#pas de lecture du header
		reading_file.readline()
		for l in reading_file : 
			#selection des resultats de blastn avec un qcov sup au seuil
			if float(l.split()[21]) > th_qcov :
				if flag == 0 :
					nb = nb + 1
					flag = 1
				#print index + "\t" + l.split()[4] + "\t" + l.split()[5] + "\t" + l.split()[6] + "\t" + l.split()[7] + "\t" + l.split()[8] + "\t" + l.split()[16] + "\t" + l.split()[21] + "\n"
				filout.write(index + "\t" + l.split()[4] + "\t" + l.split()[5] + "\t" + l.split()[6] + "\t" + l.split()[7] + "\t" + l.split()[8] + "\t" + l.split()[16] + "\t" + l.split()[21] + "\n")
	filout.close()
	print "\nIndication : Sur les", len(qcov_file_list), "especes,", nb, "ont ete selectionnées."

def add_header_ncbi_blastn_csv(ncbi_blastn_file) :
	"""
	Fonction permettant d'ajouter un header au fichier resultat blastn NCBI.
	input :
		- ncbi_blastn_file : fichier resultat de blastn a partir de NCBI WEB, format csv
	output :
		- meme fichier avec un header
	"""
	header="qseqid,qvss,qseqid,sseqid,pident,length,mismatchs,gapopen,qstart,qend,sstart,send,evalue,score"
	os.system("sed -i '1i"+header+"' "+ncbi_blastn_file)

def from_ncbi_blastn_to_qcov_file_format(csv_files, output = "ncbi_qcov_file.csv") :
	"""
	Fonction permettant d'avoir un fichier de type qcov (voir fct selection_by_qcov) à partir d'un fichier csv de ncbi web
	qseqid,qvss,qseqid,sseqid,pident,length,mismatchs,gapopen,qstart,qend,sstart,send,evalue,score
	Rmq = il faut les fichiers fasta d'alignement ayant le meme nom et etant au meme endroit que les csv_files
	input :
		- csv_files : liste de fichiers resultat de blastn a partir de NCBI WEB, format csv avec un header
			(qseqid,qvss,qseqid,sseqid,pident,length,mismatchs,gapopen,qstart,qend,sstart,send,evalue,score)
		- output : nom du fichier de sortie. 
	output :
		- creation d'un fichier contenant le resultat voulu
	"""
	#liste pour avoir id = nom du fichier
	filout = open(output, "w")
	for i in range(len(csv_files)) :
		#id = nom du fichier
		id_file = csv_files[i].split('/')[8][:-4]
		#print id_file
		#traitement du fichier csv
		c = open(csv_files[i], 'r')
		l = c.readline()
		#liste contenant les infos nec sous forme de tuples pour chaque seq
		info_csv = []
		for l in c :
			#print l[:-1]
			sseqid = l.split(',')[3]
			sstart = l.split(',')[10]
			send = l.split(',')[11]
			evalue = float(l.split(',')[12])
			qcov = query_coverage(l.split(',')[8], l.split(',')[9])
			info_csv.append((sseqid, sstart, send, evalue, qcov))
		c.close()
		#traitement du fichier fasta
		f = open(csv_files[i][:-4]+'.fasta', 'r')
		#liste contenant toutes les seq du fichier fasta
		info_fasta = []
		#chaine de car pour la seq
		seq = ""
		#pour chaque ligne, 
		for l in f:
			#si debut d'une seq, ajout de la seq precedente
			if ">" in l:
				info_fasta.append(seq)
				seq = ""
			#sinon, completer la seq
			if ">" not in l :
				#print l[:-1]
				seq = seq + l[:-1]
		#ajout de la derniere seq du fichier
		info_fasta.append(seq)
		#suppression de l'element vide
		info_fasta.remove("")
		f.close()
		for n in range(len(info_csv)) :
			#print info_csv[n]
			#print info_fasta[n]
			#print str(id_file)+"\t"+str(info_csv[n][0])+"\t"+str(info_fasta[n])+"\t"+str(info_csv[n][1])+"\t"+str(info_csv[n][2])+"\t?\t"+str(info_csv[n][3])+"\t"+str(info_csv[n][4])+"\n"
			filout.write(str(id_file)+"\t"+str(info_csv[n][0])+"\t"+str(info_fasta[n])+"\t"+str(info_csv[n][1])+"\t"+str(info_csv[n][2])+"\t?\t"+str(info_csv[n][3])+"\t"+str(info_csv[n][4])+"\n")
	filout.close()

def seq_CCMP681() :
	"""
	Fonction permettant de retourner la sequence d'un contig particulier
	"""
	CCMP681 = open("/data/mortaza/Data/transcriptomes/MMETSP1180.nt.fa","r")
	flag = 0
	for l in CCMP681 :
		if '>CAMNT_0007085561' in l :
			flag = 1
			seq = ""
			#print l
		elif '>CAMNT_0007085561' not in l and '>' in l :
			flag = 0
		if flag == 1 and '>' not in l:
			seq = seq + l[:-1]
	CCMP681.close()
	return seq

def qcov_end_file_to_fasta_file(qcov_end_file, output = "multifasta.fasta") :
	"""
	Fonction permettant de creer un fichier multifasta
	input :
		- qcov_end_file : fichier base de données (resultats de qcov) dt le header est :
		(idfile sseqid sseq sstart send strand evalue qcov)
		- output : nom donné pour le fichier de sortie multifasta
	output :
		- fichier multifasta
	"""
	fin = open(qcov_end_file,'r')
	fout = open(output, 'w')
	#Ajout de Chre = ref	
	s = "ATCCATGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGAGCAAAGCAATTTGTGTAGTGGCGAACGGGTGCGTAACGCGTAAGAACCTACCTATCGGAGGGGGATAACATTGGGAAACTGTTGCTAATACCCCATACAGCTGAGGAGTGAAAGGTGAAAAACCGCCGATAGAGGGGCTTGCGTCTGATTAGCTAGTTGGTGGGGGTAACGGCCTCCCAAGGCCACGAGCAGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGAGGAATTTTTCGCAATGGGCGCAAGCGACGGAGCAATGCCGCGTGCAGGAAGAAGGCCTGTGGGTCGTAAACTGCTTTTCTCAGAGAAGAAGTTCTGACGGTATCTGAGGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTGTCCGCAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTCGTAAAGTCTAATGTCAAATACCAGGGCTCAACCTTGGACCGGCATTGGAGTACTCACGAGCTTGAGTACGGTAGGGGCAGAGGGAATTCCATGTGGAGCGGTGAAATGCGTAGAGATATGGAGGAACACCAGTGGCGAAGGCGCTCTGCTGGGCCGAAACTGACACTGAGAGACGAAAGCTGGGGGAGCGAATAGGATTAGATACCCTAGTAGTCCCAGCCGTAAACTATGGAGACTAAGTGCTGCCGCAAGCAGTGCTGTAGCTAACGCGTTAAGTCTCCCGCCTGGGGAGTATGCTCGCAAGAGTGAAACTCAAAGGAATTGACGGGACCGCACAAGCGGTGGATTATGTGGATTAATTCGATACAACGCGAAGAACCTTACCAGGGTTTGACATGTCAAGAACCTCTCAGAAATGGGAGGGTGCCCTAACGGACTTGAACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGCTGTGAAGTGTATAGTTAAGTCTCATAACGAGCGCAACCCTCGTCTTTAGTTGCCATTTGGTTCTCTAAAGAGACTGCCAGTGTAAGCTGGAGGAAGGTGAGGATGACGTCAAGTCAGCATGCCCCTTACATCCTGGGCTTCACACGTAATACAATGGTTGGGACAATCAGAAGCGACTCGTGAGAGCTAGCGGCTCTGTTAAACCCAACCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCCAGTCAGCTATATGGCGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCACACCATGGAAGCTGGTTCTGCTCCAAGTCGTTACCCTAACCTTCGGGAGGGGGGCGCCTAAAGCAGGGCTAGTGACTAGGGTGAAGTCGTAACAAGGTAGGGCTACTGGAAGGTGGCCCTGGCTCACCTCCTTC"	
	fout.write('>Chre\n'+s+'\n')
	#l = fin.readline() #a mettre si le fichier csv contient un header
	for l in fin :
		#/!\ il y a des seq qui contiennent des gaps !!!
		seq = l.split()[2]
		seq = seq.replace('-','')
		#print '>'+str(l.split()[1])+'\n'+str(seq)+'\n'
		fout.write('>'+str(l.split()[1])+'\n'+str(seq)+'\n')
	fout.close()
	fin.close()

def add_num_seqid(database_csv, output = "new_database.csv") :
	"""
	Fonction permettant d'ajouter un num à chaque seqid pour que ca soit unique
	input : 
		- database_csv : fichier csv contenant les meilleurs hits de blastn
		- output : fichier de sortie csv (meme fichier mais changement de seqid)
	output :
		- fichier de sortie csv
	"""
	#num de seq
	n = 0
	filin = open(database_csv, 'r')
	filout = open(output, 'w')
	l = filin.readline()
	#print l
	filout.write(l)
	for l in filin :
		#print l
		n = n + 1
		#print l.split()[0]+"\t"+str(n)+"."+l.split()[1]+"\t"+l.split()[2]+"\t"+l.split()[3]+"\t"+l.split()[4]+"\t"+l.split()[5]+"\t"+l.split()[6]+"\t"+l.split()[7]+"\n"
		filout.write(l.split()[0]+"\t"+str(n)+"."+l.split()[1]+"\t"+l.split()[2]+"\t"+l.split()[3]+"\t"+l.split()[4]+"\t"+l.split()[5]+"\t"+l.split()[6]+"\t"+l.split()[7]+"\n")
	filout.close()
	filin.close()


def replace_fasta_by_name(fasta_file, id_name, output = "fasta_named") :
	"""
	Fonction permettant de remplacer les identifiants des especes par leur vrai nom
	input : 
		- fasta_file : fichier fasta
		- id_name : fichier contenant la relation entre les id et les noms des especes
		- fasta_named : nom fichier sortie fasta avec identifiants différents
	output :
		- fichier fasta de sortie
	"""
	#dico avec cle = seqid et valeur = nom de l'espece
	dico_name = {}
	#Lecture du fichier pour completer le dico
	id_name_file = open(id_name, 'r')
	#ne pas lire le header
	l = id_name_file.readline()
	for l in id_name_file :
		dico_name[l.split(',')[1]] = l.split(',')[2]
	#ouverture du fichier fasta pour lecture
	filin = open(fasta_file, 'r')
	#ouverture d'un fichier pour ecriture
	filout = open(output, 'w')
	flag = 0 #flag pour dire si ecriture ou pas de la ligne sseqid et de la ligne sseq
	for l in filin :
		if ">" in l :
			if l[1:-1] in dico_name :
				#print ">"+str(dico_name[l[1:-1]])+"\n"
				filout.write(">"+str(dico_name[l[1:-1]]))
				flag = 0
			elif "sseqid" in l :
				flag = 1  #ne pas ecrire la ligne
			else : 
				#print l
				filout.write(l)
				flag = 0 #ecrire la ligne
		else :
			#print l
			if flag == 0 :
				filout.write(l)
	filout.close()
	filin.close()

def seq_length(fasta_file, output = "length_seq.txt") :
	"""
	Fonction permettant de donner la taille des seq pstes dans un fichier fasta
	input : 
		- fasta_file : fichier fasta
		- output : fichier texte contenant deux colonnes (nom_esp, taille de la seq)
	output :
		- fichier de sortie format texte
	"""
	#ouverture d'un fichier pour ecriture
	filout = open(output, 'w')
	#ajout du header
	filout.write("name\t16s_seq_length\n")
	#ouverture d'un fichier pour lecture
	filin = open(fasta_file, 'r')
	for l in filin :
		if ">" in l :
			name = l[1:-1]
		else :
			length = len(l[:-1])
			#print name+"\t"+str(length)
			filout.write(name+"\t"+str(length)+"\n")
	#fermeture des fichiers
	filout.close()
	filin.close()

def compare_2_seq(seq1, seq2) :
	"""
	Fonction permettant de comparer deux séquences
	input :
		- seq1 : sequence 1
		- seq2 : sequence 2
	output :
		- imprime resultat de la comparaison
	"""
	for i in range(len(seq1)) :
		if seq1[i] != seq2[i] :
			print "difference entre les deux sequences\n"
			break
	print "tout est pareil"

def remove_seq(fasta_file, header_list, output = "pruned.fasta") :
	"""
	Fonction qui permet d'enlever les seq indesirables.
	input :
		- fasta_file : fichier multifasta
		- header_list : liste contenant le header des seq que l'on souhaite enlever
		- output : fichier de sortie format fasta
	output :
		- fichier multifasta
	"""
	filin = open(fasta_file, 'r')
	filout = open(output, 'w')
	flag = 0
	for l in filin :
		if ">" in l :
			if l[1:-1] not in header_list :
				#print l
				filout.write(l)
				flag = 1
		if ">" not in l and flag == 1:
			#print l
			filout.write(l)
			flag = 0
	filout.close()
	filin.close()




#programme principal
if __name__ == "__main__":
	"""
	###############################################
	genomic_file = glob.glob("/data/mortaza/Lucifer/rna16s/genomes_16S/*.csv")
	transcriptomic_file = glob.glob("/data/mortaza/Lucifer/rna16s/transcriptomes_16S/*.csv")
	ome_list = [genomic_file,transcriptomic_file]
	for ome in ome_list :
		for csv_file in ome :
			#print csv_file
			#add_header_blast(csv_file)
			#add_qcov_blast(csv_file)


	#test pour calcul de qcov
	start = [662, 994, 1259]
	end = [1068, 1314, 1342]
	for i in range(len(start)) :
		print query_coverage(start[i],end[i])
	#add_qcov_blast("/home/mortaza/Documents/16S_PhylogenyTree/G_16s_vs_GCA_000143455.1_v1.0_genomic.fna.blastn.csv")
	#add_qcov_blast("/data/mortaza/Lucifer/rna16s/genomes_16S/G_16s_vs_GCA_000147415.1_v_1.0_genomic.fna.blastn.csv") #exemple de fichier vide (sans resultat)
	###############################################
	#Selection des meilleurs hits blast avec le qcov
	g_file = glob.glob("/data/mortaza/Lucifer/rna16s/genomes_16S/qcov/*.blastn.csv.qcov")
	t_file = glob.glob("/data/mortaza/Lucifer/rna16s/transcriptomes_16S/qcov/*.blastn.csv.qcov")
	print '------> resultats du genome\n'
	selection_by_qcov(g_file, output = "/data/mortaza/Lucifer/rna16s/G_qcov_selection.csv")
	print 
	print '------> resultat du transcriptome\n'
	selection_by_qcov(t_file, output = "/data/mortaza/Lucifer/rna16s/T_qcov_selection.csv")
	#Concatenation des 2 fichiers pr rassembler les resultats omiques
	os.system("cat /data/mortaza/Lucifer/rna16s/*.csv >> /data/mortaza/Lucifer/rna16s/qcov_16S_results.csv")
	#ajout de header
	os.system("sed -i '1iid sseqid sseq sstart send sstrand evalue qcov' /data/mortaza/Lucifer/rna16s/*.csv")
	###############################################
	ncbi_results = glob.glob("/home/mortaza/Documents/16S_PhylogenyTree/New_genomes/from_NCBIblastn/blastn_results/*.csv")
	#for i in range(len(ncbi_results)) :
		#print ncbi_results[i]
		#Ajout de header dans les fichiers de blastn ncbi
		#add_header_ncbi_blastn_csv(ncbi_file)
	from_ncbi_blastn_to_qcov_file_format(ncbi_results, '/data/mortaza/Lucifer/rna16s/ncbi_qcov_file.csv')
	#MMETSP1392
	selection_by_qcov(['/data/mortaza/Lucifer/rna16s/transcriptomes_16S/qcov/T_16s_vs_MMETSP1392.nt.fa.blastn.csv.qcov'], th_qcov = 58, output = "qcov_MMETSP1392.csv")
	os.system('mv /home/mortaza/Documents/src/16S_rep/qcov_MMETSP1392.csv /data/mortaza/Lucifer/rna16s/')
	#cat en_cours.csv qcov_MMETSP1392.csv >> en_cours1.csv #lancé directement dans le terminal
	#MMETSP1180
	filout = open("/data/mortaza/Lucifer/rna16s/qcov_MMETSP1180.csv","w")
	filout.write("MMETSP1180.nt.fa\tCAMNT_0007085561\t"+str(seq_CCMP681())+"\t1\t"+str(len(seq_CCMP681()))+"\t?\t?\t?\n")
	filout.close()
	#cat en_cours1.csv qcov_MMETSP1180.csv >> database.csv #lancé directement dans le terminal
	"""
	#Creation du fichier multifasta
	#qcov_end_file_to_fasta_file("/data/mortaza/Lucifer/rna16s/database.csv", "/data/mortaza/Lucifer/rna16s/database.fasta")
	#qcov_end_file_to_fasta_file("/data/mortaza/Lucifer/rna16s/database.csv")
	##############################################
	#remplacement des id par les noms des especes dans les fichiers fasta
	#replace_fasta_by_name("/data/mortaza/Lucifer/rna16s/database.fasta", "/data/mortaza/Lucifer/rna16s/PhyML/id_name.csv", "/data/mortaza/Lucifer/rna16s/PhyML/named_database.fasta")
	#seq_length("/data/mortaza/Lucifer/rna16s/PhyML/named_database.fasta", "/data/mortaza/Lucifer/rna16s/PhyML/length_seq.txt")
	#voir la difference entre deux sequences de T. obliquus /GCA_900108755.1_Tobliquus_sob1_genomic.fna,emb|FNXT01001328.1|,T.obliquus/
	#seq1 = "CATGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGTATGCCTAACACATGCAAGTCGTACGAGATCTACTTTTTTTCGGAAAGGTAGATTGCAGTGGCGGACGGGTGAGTAATGCGTAAGAACTTACCTTTTGGTGGGGGATAACAGCGGGAAACTGCTGCTAATACCCCATAATGCTGAGGAGCTAAAAGTGAAAACTGCCAAGAGAGAGGCTTGCGTCTGATTAGCTAGTTGGTGGGGGTAAAGGCCTCCCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCACAATGGGCGAAAGCCTGATGGAGCAATACCGCGTGAGGGATGAAGCATCGTGGTGTGTAAACCTCTTTTCTTAAGGAAGAATCAATGACGGTACTTAAGGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTGTCCGGAATGATTGGGCGTAAAGCATCTGTAGGTGGTTTCTAAAGTCAACTGTTAAATCCCAGAGCTCAACTTTGGCCAAGCAGTTGAGTACTTAGAGACTTGAGTACGGTAGGGGTAGAGGGAATTCCTTGTTTAGCGGTGAAATGCATAGAGATAAGGAAGAACACCAGTGGCGAAAGCGCTCTACTGGGCCGTAACTGACACTGAGAGATGAAAGTTAGGGGAGCGAAAAGGATTAGATACCCTTGTAGTCCTAACTGTAAACGATGGATACTAAGTGCTGTCAAAAACAGTGCTGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTATGCTCGCAAGAGTGAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGATTATGCGGTTTAATTCGATGATACACGAAGAACCTTACCAGGGATTGACATGCCAGGAACTTTCTAGAGATAGATTGGTGCCTTTTTTAGGAACCTGGACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTAGAGTTAAGTCTTCTAACGAGCGCAACCCTTGTCTTTAGTTAAAACCATTTGGAACTCTAAAGAGACTGCCGGTGTAAACCGGAGGAAGGAGAGGATGACGTCAAGTCAGCATGCCCCTTACATCCTGGGCTACACGCGTAATACAATGGCTAAGACAATAGGCAGCGAACCCGCGAGGGGAAGCGAATCTAGCAAACTTAGCCTCAGTTCAGATTGTAGGCTGCAACTCGCCTGCATGAAGCCGGAATCGCTAGTAATCGCCGGTCAGCTATACGGCGGTGAATACGTTCTCGGGTCTTGTACACACCGCCCGTCACATGCTGGAAGC"
	#seq2 = "CATGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGTATGCCTAACACATGCAAGTCGTACGAGATCTACTTTTTTTCGGAAAGGTAGATTGCAGTGGCGGACGGGTGAGTAATGCGTAAGAACTTACCTTTTGGTGGGGGATAACAGCGGGAAACTGCTGCTAATACCCCATAATGCTGAGGAGCTAAAAGTGAAAACTGCCAAGAGAGAGGCTTGCGTCTGATTAGCTAGTTGGTGGGGGTAAAGGCCTCCCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCACAATGGGCGAAAGCCTGATGGAGCAATACCGCGTGAGGGATGAAGCATCGTGGTGTGTAAACCTCTTTTCTTAAGGAAGAATCAATGACGGTACTTAAGGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTGTCCGGAATGATTGGGCGTAAAGCATCTGTAGGTGGTTTCTAAAGTCAACTGTTAAATCCCAGAGCTCAACTTTGGCCAAGCAGTTGAGTACTTAGAGACTTGAGTACGGTAGGGGTAGAGGGAATTCCTTGTTTAGCGGTGAAATGCATAGAGATAAGGAAGAACACCAGTGGCGAAAGCGCTCTACTGGGCCGTAACTGACACTGAGAGATGAAAGTTAGGGGAGCGAAAAGGATTAGATACCCTTGTAGTCCTAACTGTAAACGATGGATACTAAGTGCTGTCAAAAACAGTGCTGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTATGCTCGCAAGAGTGAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGATTATGCGGTTTAATTCGATGATACACGAAGAACCTTACCAGGGATTGACATGCCAGGAACTTTCTAGAGATAGATTGGTGCCTTTTTTAGGAACCTGGACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTAGAGTTAAGTCTTCTAACGAGCGCAACCCTTGTCTTTAGTTAAAACCATTTGGAACTCTAAAGAGACTGCCGGTGTAAACCGGAGGAAGGAGAGGATGACGTCAAGTCAGCATGCCCCTTACATCCTGGGCTACACGCGTAATACAATGGCTAAGACAATAGGCAGCGAACCCGCGAGGGGAAGCGAATCTAGCAAACTTAGCCTCAGTTCAGATTGTAGGCTGCAACTCGCCTGCATGAAGCCGGAATCGCTAGTAATCGCCGGTCAGCTATACGGCGGTGAATACGTTCTCGGGTCTTGTACACACCGCCCGTCACATGCTGGAAGC"
	#seq0 = "CATGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGTATGCCTAACACATGCAAGTCGTACGAGATCTACTTTTTTTCGGAAAGGTAGATTGCAGTGGCGGACGGGTGAGTAATGCGTAAGAACTTACCTTTTGGTGGGGGATAACAGCGGGAAACTGCTGCTAATACCCCATAATGCTGAGGAGCTAAAAGTGAAAACTGCCAAGAGAGAGGCTTGCGTCTGATTAGCTAGTTGGTGGGGGTAAAGGCCTCCCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCACAATGGGCGAAAGCCTGATGGAGCAATACCGCGTGAGGGATGAAGCATCGTGGTGTGTAAACCTCTTTTCTTAAGGAAGAATCAATGACGGTACTTAAGGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTGTCCGGAATGATTGGGCGTAAAGCATCTGTAGGTGGTTTCTAAAGTCAACTGTTAAATCCCAGAGCTCAACTTTGGCCAAGCAGTTGAGTACTTAGAGACTTGAGTACGGTAGGGGTAGAGGGAATTCCTTGTTTAGCGGTGAAATGCATAGAGATAAGGAAGAACACCAGTGGCGAAAGCGCTCTACTGGGCCGTAACTGACACTGAGAGATGAAAGTTAGGGGAGCGAAAAGGATTAGATACCCTTGTAGTCCTAACTGTAAACGATGGATACTAAGTGCTGTCAAAAACAGTGCTGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTATGCTCGCAAGAGTGAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGATTATGCGGTTTAATTCGATGATACACGAAGAACCTTACCAGGGATTGACATGCCAGGAACTTTCTAGAGATAGATTGGTGCCTTTTTTAGGAACCTGGACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTAGAGTTAAGTCTTCTAACGAGCGCAACCCTTGTCTTTAGTTAAAACCATTTGGAACTCTAAAGAGACTGCCGGTGTAAACCGGAGGAAGGAGAGGATGACGTCAAGTCAGCATGCCCCTTACATCCTGGGCTACACGCGTAATACAATGGCTAAGACAATAGGCAGCGAACCCGCGAGGGGAAGCGAATCTAGCAAACTTAGCCTCAGTTCAGATTGTAGGCTGCAACTCGCCTGCATGAAGCCGGAATCGCTAGTAATCGCCGGTCAGCTATACGGCGGTGAATACGTTCTCGGGTCTTGTACACACCGCCCGTCACATGCTGGAAGC"
	#compare_2_seq(seq1, seq0)
	##############################################
	#Comme des seq st id et que dans la phylogenie, on veut enlever les seq redondantes, faut que chaque seq soit identifiable. Dc changement des seqid
	#add_num_seqid("/data/mortaza/Lucifer/rna16s/database.csv", output = "/data/mortaza/Lucifer/rna16s/new_16S_database.csv")
	#qcov_end_file_to_fasta_file("/data/mortaza/Lucifer/rna16s/new_16S_database.csv", output = "/data/mortaza/Lucifer/rna16s/new_16S_database.fasta")
	#replace_fasta_by_name("/data/mortaza/Lucifer/rna16s/new_16S_database.fasta", "/data/mortaza/Lucifer/rna16s/PhyML/id_name.csv", "/data/mortaza/Lucifer/rna16s/PhyML/named_new_16S_database.fasta")
	#seq_length("/data/mortaza/Lucifer/rna16s/PhyML/results_with_diff_seqname_database/named_new_16S_database.fasta", "/data/mortaza/Lucifer/rna16s/PhyML/results_with_diff_seqname_database/length_seq.txt")
	##############################################
	#header_list = ['E.2oleoabundans', '2Monoraphidium', 'T.1obliquus', 'T.2obliquus', 'D.1tertiolecta', 'D.3tertiolecta', 'D.4tertiolecta','V.2carteri', 'G.1pectorale']
	#remove_seq("/data/mortaza/Lucifer/rna16s/PhyML/results_with_diff_seqname_database/complet_seq/named_new_16S_database.fasta", header_list, output = "/data/mortaza/Lucifer/rna16s/PhyML/results_with_diff_seqname_database/removed_seq/16S_db_removed_seq.fasta")
	##############################################
	#Traitement du nouveau génome : T. socialis
	#ncbi_result = ["/home/mortaza/Documents/16S_PhylogenyTree/New_genomes/from_NCBIblastn/blastn_results/Tetrabaena_socialis_chloroplast_partial_genome.csv"]
	#from_ncbi_blastn_to_qcov_file_format(ncbi_result, '/data/mortaza/Lucifer/rna16s/16S_csv_blastn/inter_results/ncbi_qcov_1file.csv')
	#add_num_seqid(sys.argv[1], output = sys.argv[2])
	#creation du fichier multifasta
	#qcov_end_file_to_fasta_file(sys.argv[1], sys.argv[2])
	#replace_fasta_by_name(sys.argv[1], "/data/mortaza/Lucifer/rna16s/PhyML/results_with_diff_seqname_database/id_name.csv", sys.argv[2])
	#header_list = ['E.2oleoabundans', '2Monoraphidium', 'T.1obliquus', 'T.2obliquus', 'D.1tertiolecta', 'D.3tertiolecta', 'D.4tertiolecta','V.2carteri', 'G.1pectorale']
	#remove_seq(sys.argv[1], header_list, output = sys.argv[2])


	Seaview_tree_building_MUSCLE_complete_16S_db_removed_seq.txt