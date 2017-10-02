source("functions_alpha.R")
foldernames=c("ROIs")
##clustermethod: 1 = ours, 2 = DBScan
##bestonly: 1 = only write the label file for the best scoring proposal
sapply(foldernames, function(foldername){
r = readLines(con=file.path(paste(foldername, "/editconfig.txt", sep="")))
get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
model=get("model")
{if (model=="Gaussian(prec)"){
  xlim = as.v(get("xlim"))
  ylim = as.v(get("ylim"))
  zlim = as.v(get("zlim"))

  histbins = as.v(get("histbins"))
  histvalues = as.v(get("histvalues"))
  if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
    useplabel=FALSE; pb=NULL; alpha=NULL
  }
  else {
    useplabel=TRUE; 
    pb=as.numeric(get("pbackground"))
    alpha=as.numeric(get("alpha"))
  }
  if (length(grep("bestonly",r))==0) bestonly=FALSE
  else bestonly=as.numeric(get("bestonly"))>0
  if (length(grep("rseq",r))==0) rseq=seq(30, 200, by=10)
  else {
      rparams=as.v(get("rseq"))
      rseq=seq(rparams[1], rparams[2], by=rparams[3])
  }
  if (length(grep("thseq",r))==0) thseq=seq(0, 500, by=25)
  else {
      thparams=as.v(get("thseq"))
      thseq=seq(thparams[1], thparams[2], by=thparams[3])
  }
  if (length(grep("clustermethod",r))==0) clustermethod="K"
  else {
      method=as.numeric(get("clustermethod"))
      if (method==1) clustermethod="K"
      else clustermethod="DBSCAN"
  }
}
else {stop("Haven't implemented anything else!")}}

o = order(histbins); histbins=histbins[o]; histvalues=histvalues[o]
f = approxfun(histbins, histvalues, yleft=histvalues[1], yright=histvalues[length(histvalues)])
cst=integrate(f, lower=histbins[o],upper=histbins[length(histbins)])$value
psd <- function(sd){
  log(f(sd))-log(cst) 
}
minsd = histbins[1]; maxsd = histbins[length(histbins)]

ld=list.dirs(foldername, recursive=FALSE)
ld=ld[ld!=foldername]

for (u in 1:length(ld)){

data= read.csv(file.path(paste(foldernames, "/",u,"/editdata.txt", sep="")))
pts = data[,1:3]; sdsx = data[,4]; sdsy=data[,5]; sdsz = data[,6];
#beta parameters estimate
N = dim(data)[1];
zvalues = data[,3];
#scaled z values
zfix = (zvalues - zlim[1])/(zlim[2]-zlim[1]);
samplemean = (1/N) * sum(zfix);
samplevar = (1/(N-1)) * sum((zfix-samplemean)^2)

paramc = samplemean * (((samplemean*(1-samplemean))/samplevar) - 1);
paramd = (1-samplemean) * (((samplemean*(1-samplemean))/samplevar) - 1);

res=Kclust(pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb, score=TRUE, rlabel=TRUE, report=TRUE, rseq=rseq, thseq=thseq)
writeRes(res, file.path(paste(foldernames, "/",u, "/r_vs_thresh.txt", sep="")), file.path(paste(foldernames, "/",u, "/labels", sep="")), bestonly=bestonly)

}


})
