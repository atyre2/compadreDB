## predicting vital rates across phylogenies
load("data/COMADRE_v.1.0.0.RData")
library(Mage)

amniotes <- read.csv("workshop/Amniote_Database.csv")

ssdb <- subsetDB(db = comadre, sub = MatrixDimension == 2 & 
                   Class == "Aves" &
                   MatrixComposite != "Individual" & 
                   MatrixTreatment == "Unmanipulated")

ssdb$metadata$MatrixID <- row.names(ssdb$metadata)

#unique(ssdb$metadata$SpeciesAccepted) %in% with(amniotes,paste(genus, species, sep = "_"))

amniotes$SpeciesAccepted <- with(amniotes,paste(genus, species, sep = "_"))
# OK... so I should be able to match them up ... 

ssmd <- merge(ssdb$metadata, amniotes, all.x = TRUE)

library(popdemo)

ssmd$lambda <- sapply(ssdb$mat,function(X)popbio::lambda(X$matA))
ssmd$reactivity <- sapply(ssdb$mat,function(X)reactivity(X$matA))
ssmd$adult_survival <- sapply(ssdb$mat, function(X) X$matA[2,2])

with(ssmd,plot(adult_survival,adult_body_mass_g))
