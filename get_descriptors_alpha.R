source("functions_alpha.R")
library(scatterplot3d)
library(rgl)
library(plot3D)
foldernames=c("ROIs")

superplot=TRUE
makeplot=TRUE
skeleton=FALSE

sapply(foldernames, function(expname){
  nexpname=expname
  if (skeleton){
    nexpname=paste("R_", expname, sep="")
    dir.create(file.path(nexpname))
    file.copy(file.path(paste(expname, "/editconfig.txt", sep="")), paste(nexpname, "/editconfig.txt", sep=""))
    file.copy(file.path(paste(expname, "/editsim_params.txt", sep="")), paste(nexpname, "/editsim_params.txt", sep=""))
  }
  
  id=1
  id<<-id
  r = readLines(con=file.path(paste(nexpname, "/editconfig.txt", sep="")))
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
  }
  else {stop("Haven't implemented anything else!")}}
  
  all=list.files(expname)
  dirnames=all[file.info(file.path(paste(expname, "/", all, sep="")))$isdir]
  
  
  
  axes=FALSE;
  cex=1/(3*sqrt(length(dirnames)))
  if (makeplot & superplot) {
    ##These settings control image size and resolution
    png(file.path(paste(nexpname, "/together.png",sep="")), width=10, height=10, units="cm", res=1200)
    nrow=ceiling(sqrt(length(dirnames)))
    par(mfrow=c(nrow, nrow))
  }
  
  
  
  
  #res=sapply(dirnames, function(dirname){
    
  print('length(dirnames)')
  print(length(dirnames))  
  position=c()
  position<<-position
  #step of 15 up to 600 (40 bars) first line is the sum, second is the number of time sampled
  hist_position=c()
  #step of 30
  hist_position=matrix(0,2,20)
  hist_position<<-hist_position
  nmols=c()
  nmols<<-nmols
  radii=c()
  radii<<-radii
  nclusters=c()
  nclusters<<-nclusters
  percentage=c()
  percentage<<-percentage
  local_percentage=matrix(0,2,20)
  local_percentage<<-local_percentage
      
  for (dirname in 1:length(dirnames)){  
  #New loop as ordered should be here  
    print(zlim)
    coeff=zlim[2]/30
    coeff=trunc(coeff)+1
    if (coeff>20){
      coeff=20 
    }

    foldername=file.path(paste(expname, "/", dirname, sep=""))
    nfoldername=file.path(paste(nexpname, "/", dirname, sep=""))
    if (skeleton){
      dir.create(nfoldername)  
      file.copy(file.path(paste(foldername, "/editdata.txt", sep="")), file.path(paste(nfoldername, "/editdata.txt", sep="")))
    }
    data=read.csv(file.path(paste(nfoldername, "/editdata.txt", sep="")))
    
    pts = data[,1:3]; sdsx = data[,4]; sdsy = data[,5]; sdz = data[,6];
    if (skeleton){
      file.copy(file.path(paste(foldername, "/r_vs_thresh.txt", sep="")), file.path(paste(nfoldername, "/r_vs_thresh.txt", sep="")))
    }
    r = read.csv(file.path(paste(nfoldername, "/r_vs_thresh.txt",sep="")), header=FALSE, sep="\t")
    
    m = as.matrix(r)
    cs=(m[1,])[-1]
    thr=(m[,1])[-1]
    m = as.matrix(m[2:length(m[,1]),2:length(m[1,])])
    which.maxm <- function(mat){
      indcol <- rep(1:ncol(mat), each=nrow(mat))[which.max(mat)] 
      indrow <- rep(1:nrow(mat), ncol(mat))[which.max(mat)]
      c(indrow, indcol)
    }
    best=which.maxm(m)
    bestcs=cs[best[2]]
    bestthr=thr[best[1]]
    bfile=file.path(paste(foldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
    nbfile=bfile
    if (skeleton){
      dir.create(paste(nfoldername, "/labels", sep=""))    
      nbfile=file.path(paste(nfoldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
      file.copy(bfile, nbfile)
    }
    labelsbest = strsplit(readLines(nbfile),",")[[1]]
    labelbest = matrix(c(labelsbest), ncol=1)
    blabel = as.matrix(dist(labelbest))
    blabel[blabel > 0] <- -1
    blabel[blabel == 0] <- 1
    blabel[blabel < 0] <- 0
    
    write.csv(blabel, file=file.path(paste(foldername, "/blabel.csv", sep="")), row.names=FALSE, quote=FALSE)
    
    ##Some summaries
    wfile=file.path(paste(nfoldername, "/summary.txt", sep=""))
    cat("The best: clusterscale", bestcs, " thresh", bestthr, "labels.txt\nNumber of clusters:", nClusters(labelsbest), "\nPercentage in clusters: ", percentageInCluster(labelsbest), "%\nMean number of molecules per cluster: ", nMolsPerCluster(labelsbest), "\nMean radius: ", mean(clusterRadii(pts, labelsbest)), sep="", file=wfile)
   
    pdf(file.path(paste(foldername, "/bestplot.pdf", sep="")))
    scatterplot3d(pts, xlim=xlim, ylim=ylim, zlim=zlim, color = mkcols(labelsbest), sub="Best labels", xlab="x",ylab="y", zlab='z')
    dev.off()
    #Plots-interactive
    open3d(windowRect=c(250,250,1750,1750), zoom=0.65)
    plot3d(pts, xlab='x', ylab='y', zlab='z', col=mkcols(labelsbest), xlim=xlim, ylim=ylim, zlim=zlim, sub="Best", aspect=c(3,3,1), type = 's', size = 0.75)
    rgl.snapshot(file.path(paste(foldername, "/Bestplot.png", sep="")), fmt = "png", top = TRUE)
    rgl.close()
    
    #Descriptors
    s=length(labelsbest)
    n_cluster=nClustersi(labelsbest)
    
    nclusters=c(nclusters,nClusters(labelsbest))
    nclusters<<-nclusters
    percentage=c(percentage,percentageInCluster(labelsbest))
    percentage<<-percentage
    
    for (h in 1:coeff){
      hist_position[2,h]=hist_position[2,h]+1
      
      labelsbest_tempo=c()
      
      for (g in 1:length(labelsbest)){
        if ((pts[g,3]<=(h*30))& (pts[g,3]>=((h-1)*30))){
          labelsbest_tempo=c(labelsbest_tempo,labelbest[g])
        }
      }
      if (length(labelsbest_tempo)>1){
        local_percentage[2,h]=local_percentage[2,h]+1
        local_percentage[1,h]=local_percentage[1,h]+percentageInCluster(labelsbest_tempo)
        
      }
      
      local_percentage<<-local_percentage
      
    }
    hist_position<<-hist_position
    local_percentage<<-local_percentage
    eq_coeff=30*coeff

    n_cluster=c()
    n_cluster=nClustersi(labelsbest)
    for (k in 1:nClusters(labelsbest)){
      Coord=c()
      pp=0
      l=n_cluster[k]
      for (ll in 1:s){
        if (labelsbest[ll]==l){
          if (pp==0){
            Coord=pts[ll,]
            pp=1
          }else{
            
            Coord=rbind(Coord,pts[ll,])
          }
        }
      }
      av_x=mean(Coord[,1])
      av_y=mean(Coord[,2])
      av_z=mean(Coord[,3])
      
      H=0
      H=molsPerClusteri(labelsbest,l)
      J=0
      J=clusterRadiii(pts, labelsbest, l)
      
      
      A=c(av_x,av_y,av_z)
      
      if (av_z<=eq_coeff){
        which_bar=trunc(av_z/30)+1
        hist_position[1,which_bar]=hist_position[1,which_bar]+1/nClusters(labelsbest)
      }
      hist_position<<-hist_position
      position=rbind(position,A)
      position<<-position
      nmols=c(nmols,H)
      nmols<<-nmols
      radii=c(radii,J)
      radii<<-radii
    }
    

    
    
    #list(radii=clusterRadii(pts, labelsbest), nmols=molsPerCluster(labelsbest), nclusters=nClusters(labelsbest), percentageclustered=percentageInCluster(labelsbest))
   
    }

  #h=hist(nmols, breaks = seq(0,500,10), plot=FALSE)
  #pdf(file.path(paste(nexpname, "/nmols.pdf", sep="")))
  #plot(h, xlab="Number of molecules per cluster", ylab="Number of clusters", main="")
  #dev.off()
  f=file.path(paste(nexpname, "/nmols.txt", sep="")); cat(nmols, file=f, sep=","); cat("\n", file=f, append=TRUE)
  dev.off()
  write.csv(position,file = paste(nexpname, "/position.csv", sep=""), sep=",")
  write.csv(hist_position,file = paste(nexpname, "/hist_position.csv", sep=""), sep=",")
  write.csv(local_percentage,file = paste(nexpname, "/percentage_in_local.csv", sep=""), sep=",")
  #h=hist(radii, breaks = seq(0,250,10), plot=FALSE)
  #pdf(file.path(paste(nexpname, "/radii.pdf", sep="")))
  #plot(h, xlab="Cluster radius", ylab="Number of clusters", main="")
  #dev.off()
  f=file.path(paste(nexpname, "/radii.txt", sep="")); cat(radii, file=f, sep=","); cat("\n", file=f, append=TRUE)
  
  
  h=hist(nclusters, breaks = seq(0,50,2), plot=FALSE)
  pdf(file.path(paste(nexpname, "/nclusters.pdf", sep="")))
  plot(h, xlab="Number of clusters", ylab="Number of regions", main="")
  dev.off()
  f=file.path(paste(nexpname, "/nclusters.txt", sep="")); cat(nclusters, file=f, sep=","); cat("\n", file=f, append=TRUE)
  
  
  h=hist(percentage, breaks = seq(0,100,5), plot=FALSE)
  pdf(file.path(paste(nexpname, "/pclustered.pdf", sep="")))
  plot(h, xlab="Percentage clustered", ylab="Number of regions", main="")
  dev.off()
  f=file.path(paste(nexpname, "/pclustered.txt", sep="")); cat(percentage, file=f, sep=","); cat("\n", file=f, append=TRUE)   
  


  })