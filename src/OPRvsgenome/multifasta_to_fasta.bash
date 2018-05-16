#!/usr/bin/env bash

#"""
#Author : Shogofa MORTAZA
#IBPC Internship
#Passer d'un fichier multifasta à plusieurs fichiers fasta
#exemple de commande lancée :
#/home/mortaza/Documents/src/OPRvsgenome/multifasta_to_fasta.bash OPRs.fasta
#"""

#multifasta file en fasta files
awk '/^>/{s=++d".fasta"} {print > s}' $1 #ligne de commande trouvée sur Internet
#cration d'un rep contenant l'ensemble des fichiers fasta créés
mkdir fasta_files
mv *.fasta fasta_files
cd fasta_files
mv $1 ..

#Création d'un fichier intermediaire contenant l'ensemble des noms de fichiers fasta
ls >> fasta_name.txt
sed -i '/fasta_name.txt/d' fasta_name.txt
#Modification des noms de fichiers par les headers des fastas
LIST='fasta_name.txt'
while read ligne
do
	#echo $ligne
	NAME=`grep ">" $ligne`
	NAME=${NAME:1}
	#echo $NAME.fasta
	mv $ligne $NAME.fasta
done < $LIST

#Effacer le fichier intermediaire
rm fasta_name.txt
