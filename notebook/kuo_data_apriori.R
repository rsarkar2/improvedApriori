# Kuo et al setting
# 13 reports, 13 AEs and 5 drugs

library(readxl)
tab <- read_excel("/Users/reetikasarkar/Documents/dissertation_research/apriori_project/codes copy/data/kuo_data.xlsx")
drugs <- c()
aes <- c()

#make list of all drugs and aes in the data
for (row in 1:nrow(tab)) {
  drugs <- append(drugs, unlist(strsplit(as.character(tab[row,2]), ", ")))
  aes <- append(aes, unlist(strsplit(as.character(tab[row,3]), ", ")))
}

unique_drugs <- unique(drugs)
unique_aes <- unique(aes)

#restructure data into array of binary variables corresponding to drugs and aes
dat_drug <- matrix(0, nrow=nrow(tab), ncol=length(unique_drugs))
colnames(dat_drug) <- unique_drugs

dat_ae <- matrix(0, nrow=nrow(tab), ncol=length(unique_aes))
colnames(dat_ae) <- unique_aes

for (row in 1:nrow(tab)) {
  for (dr in 1:length(unique_drugs)) {
    if (unique_drugs[dr] %in% unlist(strsplit(as.character(tab[row,2]), ", "))) {
      dat_drug[row,dr] = 1
    }
  }
  for (ae in 1:length(unique_aes)) {
    if (unique_aes[ae] %in% unlist(strsplit(as.character(tab[row,3]), ", "))) {
      dat_ae[row,ae] = 1
    }
  }
}

kuo_data <- cbind.data.frame(dat_drug, dat_ae)
num_drugs <- length(unique_drugs)
num_ae <- length(unique_aes)

num_reports <- nrow(kuo_data)

