#get start time
time.start = Sys.time()
#import and process gene expression data
raw = read.csv(
  "leukemia.txt",
  header = FALSE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE
)
raw = t(raw)
row.names(raw) = NULL
raw = as.data.frame(raw, stringsAsFactors = FALSE)
for (i in 2:ncol(raw)) {
  raw[, i] = as.integer(raw[, i])
}

#subset ALL and AML
raw.all = raw[raw$`gene/patient` == "ALL", 2:ncol(raw)]
raw.aml = raw[raw$`gene/patient` == "AML", 2:ncol(raw)]

#create a list of genes and a new data frame
genes = colnames(raw[2:ncol(raw)])
test = as.data.frame(genes, stringsAsFactors = FALSE)
colnames(test)[1] = "gene"

#use t-test too score and rank the differential gene expression
for (i in genes) {
  test$p.value[test$gene == i] = t.test(raw.all[, i], raw.aml[, i])$p.value
}
for (i in test$gene) {
  test$mean.all[test$gene == i] = mean(raw.all[, i])
}
for (i in test$gene) {
  test$mean.aml[test$gene == i] = mean(raw.aml[, i])
}
test$diff = test$mean.all - test$mean.aml
test$sign[test$diff > 0] = 1
test$sign[test$diff < 0] = -1
test$psign = test$p.value * test$sign
test$rank = rank(test$p.value)

#calculate phenotypic relationship
test$qvalue[test$sign > 0] = log10(test$p.value[test$sign > 0]) / log10(min(test$p.value[test$sign >
                                                                                           0])) * test$sign[test$sign > 0]
test$qvalue[test$sign < 0] = log10(test$p.value[test$sign < 0]) / log10(min(test$p.value[test$sign <
                                                                                           0])) * test$sign[test$sign < 0]
test$ranksign = rank(-test$qvalue)
barplot(test$qvalue[order(test$ranksign)])

#import pathways
pathway = read.csv(
  file = "pathways.txt",
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)
#get pathway with more than 15 genes
for (i in row.names(pathway)) {
  if (length(as.character(pathway[row.names(pathway) == i, 3:21])[as.character(pathway[row.names(pathway) ==
                                                                                       i, 3:21]) != ""]) < 15) {
    pathway = pathway[!row.names(pathway) == i, ]
  }
}
report = as.data.frame(cbind(geneset = pathway[, 1]), stringsAsFactors = FALSE)

#for(k in report$geneset[c(582,301,662,300,596)]){
for (k in report$geneset) {
  pathwaygene = c(as.character(pathway[pathway$V1 == k, 3:21]))
  
  #create list of hit genes
  hitgene = NULL
  for (i in genes) {
    if (i %in% pathwaygene == TRUE) {
      hitgene = c(hitgene, i)
    }
  }
  #calculate ES
  NR = 0
  p = 1
  for (i in hitgene) {
    NR = NR + abs(test$qvalue[test$gene == i]) ^ p
  }
  ranks = order(test$ranksign)
  phit = 0
  pmiss = 0
  ESs = NULL
  for (i in ranks) {
    if (test$gene[i] %in% hitgene == TRUE) {
      phit = phit + abs(test$qvalue[i]) ^ p / NR
    }
    if (test$gene[i] %in% hitgene == FALSE) {
      pmiss = pmiss + 1 / (length(genes) - length(hitgene))
    }
    ESs[i] = phit - pmiss
  }
  test$ES = ESs
  #plot(test$ranksign,test$ES)
  test$ESrank = rank(-test$ES)
  test$absES = abs(test$ES)
  test$absESrank = rank(-test$absES)
  ES = test$ES[test$absESrank == 1]
  report$ES[report$geneset == k] = ES
  #calculate randomized ES
  nper = 10
  per = as.data.frame(matrix(nrow = nper, ncol = length(ranks)))
  colnames(per) = ranks
  esper = NULL
  rankper = as.data.frame(rbind(ranks))
  for (j in 1:nper) {
    rankper[j, ] = sample(ranks)
  }
  for (j in 1:nper) {
    phit = 0
    pmiss = 0
    m = rankper[j, ]
    for (i in m) {
      if (test$gene[i] %in% hitgene == TRUE) {
        phit = phit + abs(test$qvalue[i]) ^ p / NR
      }
      if (test$gene[i] %in% hitgene == FALSE) {
        pmiss = pmiss + 1 / (length(genes) - length(hitgene))
      }
      per[j, rownames(per) == i] = phit - pmiss
    }
    esper = c(esper, per[j, max.col(abs(per[j, ]))])
  }
  if (ES >= 0) {
    esmean = mean(esper[esper >= 0])
  }
  if (ES < 0) {
    esmean = mean(esper[esper < 0])
  }
  report$esmean[report$geneset == k] = esmean
  #calculate p-value for ES
  if (ES >= 0) {
    esp = 1 - ecdf(esper[esper >= 0])(ES)
  }
  if (ES < 0) {
    esp = ecdf(esper[esper < 0])(ES)
  }
  report$pvalue[report$geneset == k] = esp
}
report$NES = report$ES / report$esmean
report = report[order(-report$ES), ]
print(report[1:20, ])
write.csv(report, file = "report.csv")
#get finish time
time.finish = Sys.time()
time.elapse = time.finish - time.start
sprintf("Start Time: %s", time.start)
sprintf("Finish Time: %s", time.finish)
sprintf("Time Elapse: %s", time.elapse)