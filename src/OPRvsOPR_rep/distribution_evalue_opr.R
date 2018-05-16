#ouverture du fichier resultat blast (OPR1, OPR2, Evalue)
data <- read.table("opr_vs_opr_cytoscape.txt")
#recuperation des Evalues
Evalue <- data$V3
#recuperation du min des Evalues diff de 0
min <- min(Evalue[which(Evalue != 0)])

log.Evalue <- rep()
for (i in Evalue){
  if (i == 0){
    print (c(i, -log(min/10))) 
    log.Evalue <- rep(c(log.Evalue, -log(min/10)))
  }else{
    print(c(i,-log(i)))
    log.Evalue <- rep(c(log.Evalue, -log(i)))
  }
}

#verification de la taille des vecteurs
length(Evalue)
length(log.Evalue)

#ajout de la colonne à data
data <- cbind(data,log.Evalue)
#enregistrement de la nouvelle table sans les noms de colonnes.
write.table(data, "opr_vs_opr_logEvalue.txt", col.names = FALSE, row.names = FALSE)

#distribution des Evalues
hist(Evalue)
#distribution des -log(Evalues) diff de 0
hist(Evalue[which(Evalue != 0)])
