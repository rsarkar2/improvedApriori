#' R Package for implementing improvedApriori to identify frequently co-occurring drug/vaccine and AE pair
#'
#' This is an R package to identify frequently co-occurring drug/vaccine and adverse event (AE) pairs from safety datasets using either the traditional Apriori method ("Confidence" as second hyper parameter) or an improved Apriori method (using "ROR", "RR", or "PRR" as the second hyper parameter)
#' @docType _PACKAGE
#' @name improvedApriori
NULL
#' @param data The binarized drug/vaccine safety data. Data should be saved as a data frame.
#' @param numdrug The number of unique drugs/vaccines. Default value is 5.
#' @param numae The number of unique adverse events. Default value is 13.
#' @param numreports The number of unique report IDs. Default value is 13.
#' @param second_par The choice for second hyper parameter, i.e., ROR, PRR, RR, Confidence.
#' @param minsup The threshold for the first hyper parameter - Support. Default value is 3.
#' @param minthresh The threshold for the second hyper parameter. Default value is 1.
#' @return Data frame containing frequently co-occurring drugs/vaccines and AEs.
#' @export
#' 
#' @examples improvedApriori(data=kuo_data, numdrug=5, numae=13, numreports=13, second_par="ROR", minsup=3, minthresh=1)
#' 
# improved Apriori function
# data needs to be binarized before applying function
# each row represents a unique report
improvedApriori <- function(data=data, numdrug=5, numae=13, numreports=13, second_par="ROR", minsup=3, minthresh=1){
  
  # make a list of unique VAERS IDs, vaccine and ae codes in the data
  unique_vacs <- colnames(data)[1:numdrug]
  unique_aes <- colnames(data)[(numdrug+1):ncol(data)]
  
  #1 itemsets (include all drugs and AEs)
  c1 <- data.frame(item=character(), count=integer(), stringsAsFactors = FALSE) #1 item candidate sets
  for (col in 1:ncol(data)) {
    c1[col,1] <- colnames(data)[col]
    c1[col,2] <- sum(data[,col])
  }
  
  if (nrow(c1) == 0) {
    print("no frequent 1 itemsets")
  } else{
    l1 <- data.frame(item=character(), count=integer(), 
                     stringsAsFactors = FALSE) #1 item item sets
    for (item in 1:nrow(c1)) {
      if((as.integer(c1[item,2])) >= minsup) {
        l1 <- rbind(l1, c1[item,])
      }
    }
    rownames(l1) <- 1:nrow(l1)
  }
  
  # 2 itemsets (only pairs with one drug/vaccine and one AE to be considered)
  c2 <- data.frame(item1=character(), item2=character(), count=integer(), stringsAsFactors = FALSE) #2 item candidate sets
  for (i in 1:(nrow(l1)-1)) {
    for (j in (i+1):nrow(l1)) {
      if(((l1[i,1] %in% unique_vacs) & (l1[j,1] %in% unique_aes)) | ((l1[i,1] %in% unique_aes) & (l1[j,1] %in% unique_vacs))){
        com <- c(l1[i,1], l1[j,1], sum((data[,l1[i,1]]==1) & (data[,l1[j,1]]==1)))
        c2 <- rbind(c2, com)
      }
    }
  }
  names(c2) <- c("item1","item2","count")
  
  if (nrow(c2) == 0) {
    print("no frequent 2 itemsets")
  } else{
    l2 <- data.frame(item1=character(), item2=character(), count=integer(), 
                     stringsAsFactors = FALSE) #2 item item sets
    for (item in 1:nrow(c2)) {
      if(as.integer(c2[item,ncol(c2)]) >= minsup) {
        l2 <- rbind(l2, c2[item,])
      }
    }
    if (nrow(l2)!=0){
      names(l2) <- c("item1","item2","count")
      rownames(l2) <- 1:nrow(l2)
    }
  }
  
  freq2_selected <- data.frame(item1=character(), item2=character(), support=integer(), confidence=numeric())
  i1 <- 0
  cf <- 0
  
  if (nrow(l2)==0){
    print("No frequent 2 itemsets")
  } else {
    for (row in 1:nrow(l2)) {
      a <- sum((data[,l2[row,1]] == 1) & (data[,l2[row,2]] == 1))
      b <- sum((data[,l2[row,1]] == 1) & (data[,l2[row,2]] == 0))
      c <- sum((data[,l2[row,1]] == 0) & (data[,l2[row,2]] == 1))
      d <- sum((data[,l2[row,1]] == 0) & (data[,l2[row,2]] == 0))
      
      prr <- (a/(a + b + 0.00001))/((c+0.00001)/(c+d+0.00001))
      rr <- num_reports*a/((a + c)*(a + b)+0.00001)
      ror <- (a*d)/(c*b+0.00001)
      
      if(second_par == "ROR"){
        if(ror >= minthresh){
          freq2_selected <- rbind(freq2_selected, c(unlist(l2[row,]), ror))
        }
        colnames(freq2_selected) <- c("item1","item2","support","ror")
      } else if(second_par == "RR"){
        if(rr >= minthresh){
          freq2_selected <- rbind(freq2_selected, c(unlist(l2[row,]), rr))
        }
        colnames(freq2_selected) <- c("item1","item2","support","rr")
      } else if(second_par == "PRR"){
        if(prr >= minthresh){
          freq2_selected <- rbind(freq2_selected, c(unlist(l2[row,]), prr))
        }
        colnames(freq2_selected) <- c("item1","item2","support","prr")
      } else if(second_par == "Confidence"){
        i1 <- which(l1[,1] == l2[row,1])
        cf <- as.numeric(l2[row,3])/as.numeric(l1[i1,2])
        if (cf >= minthresh) {
          freq2_selected <- rbind(freq2_selected, c(unlist(l2[row,]), cf))
        }
        colnames(freq2_selected) <- c("item1","item2","support","confidence")
      } else print("Please select a valid second parameter.")
    } 
  }
  print(freq2_selected)
}

#' R function for binarizing safety data for implementing improvedApriori function
#' @rdname binarizeData
#' @export
binarizeData <- function(file=file, numdrug=5, numae=13, numreports=13){
tab <- read_excel(file)
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

binarydata <- cbind.data.frame(dat_drug, dat_ae)
return(binarydata)
}