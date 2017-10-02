library(splancs)
library(igraph)

mcgaussprec3sd <- function(pts, sdsx, sdsy, sdsz, xlim = c(0,1), ylim=c(0,1), zlim=c(0,1), paramc, paramd, psd=function(sd){0},minsd=0.1,maxsd=100, grid=100){
  N = dim(pts)[1]
  fsd <- Vectorize(function(sd){
      wtsx = 1/(sd^2 + sdsx^2); tildeNx = sum(wtsx);
      wtsy = 1/(sd^2 + sdsy^2); tildeNy = sum(wtsy);
      wtsz = 1/(sd^2 + sdsz^2); tildeNz = sum(wtsz);
      mu = c(sum(wtsx * pts[,1])/tildeNx, sum(wtsy * pts[,2])/tildeNy, sum(wtsz * pts[,3])/tildeNz)
      totdist = sum(c(wtsx*(pts[,1] - mu[1])^2, wtsy*(pts[,2]-mu[2])^2, wtsz*(pts[,3] - mu[3])^2))    

    #numerical integration for z axis
    zfunc <- function(zmean) {(dbeta((zmean - zlim[1])/(zlim[2]-zlim[1]), paramc, paramd, log = FALSE)) * exp(-0.5 * tildeNz * ((zmean - mu[3])*(zmean - mu[3])))}
    zint = integrate(zfunc, zlim[1], zlim[2])
    
    ##x-axis
    log(pnorm(sqrt(tildeNx) * (xlim[2]-mu[1])) - pnorm(sqrt(tildeNx) * (xlim[1] - mu[1])))+
      ##y-axis
      log(pnorm(sqrt(tildeNy) * (ylim[2]-mu[2])) - pnorm(sqrt(tildeNy) * (ylim[1] - mu[2])))+
      ##z-axis
      log(zint$value)+ 
      ##marginalised (with factor 2pi taken from standardisation above)
      -(1.5*N - 1)*log(2 * pi)+.5*(sum(log(wtsx))+sum(log(wtsz))+sum(log(wtsy)))-totdist/2+
      ##size of area
      -log(diff(xlim)*diff(ylim))+
      ##cst -- left in standardisation above
      -1/2*log(tildeNx)-1/2*log(tildeNy)+
      ##prior on sd
      psd(sd)
  })
  ##discrete prior:
  x = seq(minsd, maxsd, length=grid)[-1];
  values = fsd(x); dx=x[2]-x[1]; m = max(values); int=sum(exp(values-m))*dx;
  log(int)+m
}


mkcols <- function(labels){
  t=table(labels)
  cnames=names(t[t>1])  
  colors=sample(rainbow(length(cnames)))
  s=sapply(labels, function(l){
    i=which(names(t)==l);
    if (t[i]==1){"grey"}
    else {colors[which(cnames==l)]}
  })
  s  
}

toroid <- function(pts, xlim, ylim, zlim, range){
  xd=xlim[2]-xlim[1]; yd=ylim[2]-ylim[1]; zd=zlim[2]-zlim[1];
  
  R=pts[pts[,1] >= (xlim[2]-range),,drop=FALSE]; Rshift=t(apply(R, 1, function(v) {v - c(xd, 0, 0)}))
  L=pts[pts[,1] <= range,,drop=FALSE]; Lshift=t(apply(L, 1, function(v) {v + c(xd, 0, 0)}))
  U=pts[pts[,2] >= ylim[2]-range,,drop=FALSE]; Ushift=t(apply(U, 1, function(v) {v - c(0, yd, 0)}))
  D=pts[pts[,2] <= range,,drop=FALSE]; Dshift=t(apply(D, 1, function(v) {v + c(0, yd, 0)}))
  Fw=pts[pts[,3] >= zlim[2]-range,,drop=FALSE]; Fwshift=t(apply(Fw, 1, function(v) {v - c(0, 0, zd)}))
  B=pts[pts[,3] <= range,,drop=FALSE]; Bshift=t(apply(B, 1, function(v) {v + c(0, 0, zd)}))
  
  UFw=pts[(pts[,2] >= ylim[2]-range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; UFwshift=t(apply(UFw, 1, function(v) {v + c(0, -yd, -zd)}))
  DFw=pts[(pts[,2] <= range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; DFwshift=t(apply(DFw, 1, function(v) {v + c(0, yd, -zd)}))
  LFw=pts[(pts[,1] <= range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; LFwshift=t(apply(LFw, 1, function(v) {v + c(xd, 0, -zd)}))
  RFw=pts[(pts[,1] >= xlim[2]-range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; RFwshift=t(apply(RFw, 1, function(v) {v + c(-xd, 0, -zd)}))
  UB=pts[(pts[,2] >= ylim[2]-range) & (pts[,3] <= range),,drop=FALSE]; UBshift=t(apply(UB, 1, function(v) {v + c(0, -yd, zd)}))
  DB=pts[(pts[,2] <= range) & (pts[,3] <= range),,drop=FALSE]; DBshift=t(apply(DB, 1, function(v) {v + c(0, yd, zd)}))
  LB=pts[(pts[,1] <= range) & (pts[,3] <= range),,drop=FALSE]; LBshift=t(apply(LB, 1, function(v) {v + c(xd, 0, zd)}))
  RB=pts[(pts[,1] >= xlim[2]-range) & (pts[,3] <= range),,drop=FALSE]; RBshift=t(apply(RB, 1, function(v) {v + c(-xd, 0, zd)}))
  LU=pts[(pts[,1] <= range) & (pts[,2] >= ylim[2]-range),,drop=FALSE]; LUshift=t(apply(LU, 1, function(v) {v + c(xd, -yd, 0)}))
  RU=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] >= ylim[2]-range),,drop=FALSE]; RUshift=t(apply(RU, 1, function(v) {v + c(-xd, -yd, 0)}))
  LD=pts[(pts[,1] <= range) & (pts[,2] <= range),,drop=FALSE]; LDshift=t(apply(LD, 1, function(v) {v + c(xd, yd, 0)}))
  RD=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] <= range),,drop=FALSE]; RDshift=t(apply(RD, 1, function(v) {v + c(-xd, yd, 0)}))
  
  LUB=pts[(pts[,1] <= range) & (pts[,2] >= ylim[2]-range) & (pts[,3] <= range),,drop=FALSE]; LUBshift=t(apply(LUB, 1, function(v) {v + c(xd, -yd, zd)}))
  RUB=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] >= ylim[2]-range) & (pts[,3] <= range),,drop=FALSE]; RUBshift=t(apply(RUB, 1, function(v) {v + c(-xd, -yd, zd)}))
  LUFw=pts[(pts[,1] <= range) & (pts[,2] >= ylim[2]-range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; LUFwshift=t(apply(LUFw, 1, function(v) {v + c(xd, -yd, -zd)}))
  RUFw=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] >= ylim[2]-range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; RUFwshift=t(apply(RUFw, 1, function(v) {v + c(-xd, -yd, -zd)}))
  LDB=pts[(pts[,1] <= range) & (pts[,2] <= range) & (pts[,3] <= range),,drop=FALSE]; LDBshift=t(apply(LDB, 1, function(v) {v + c(xd, yd, zd)}))
  RDB=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] <= range) & (pts[,3] <= range),,drop=FALSE]; RDBshift=t(apply(RDB, 1, function(v) {v + c(-xd, yd, zd)}))
  LDFw=pts[(pts[,1] <= range) & (pts[,2] <= range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; LDFwshift=t(apply(LDFw, 1, function(v) {v + c(xd, yd, -zd)}))
  RDFw=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] <= range) & (pts[,3] >= zlim[2]-range),,drop=FALSE]; RDFwshift=t(apply(RDFw, 1, function(v) {v + c(-xd, yd, -zd)}))
  
  
  if (length(Rshift)>0)
    pts=rbind(pts, Rshift)
  if (length(Lshift)>0)
    pts=rbind(pts, Lshift)
  if (length(Ushift)>0)
    pts=rbind(pts, Ushift)
  if (length(Dshift)>0)
    pts=rbind(pts, Dshift)
  if (length(Fwshift)>0)
    pts=rbind(pts, Fwshift)
  if (length(Bshift)>0)
    pts=rbind(pts, Bshift)
  
  if (length(LUBshift)>0)
    pts=rbind(pts, LUBshift)
  if (length(RUBshift)>0)
    pts=rbind(pts, RUBshift)
  if (length(LUFwshift)>0)
    pts=rbind(pts, LUFwshift)
  if (length(RUFwshift)>0)
    pts=rbind(pts, RUFwshift)
  if (length(LDBshift)>0)
    pts=rbind(pts, LDBshift)
  if (length(RDBshift)>0)
    pts=rbind(pts, RDBshift)
  if (length(RDFwshift)>0)
    pts=rbind(pts, RDFwshift)
  if (length(LDFwshift)>0)
    pts=rbind(pts, LDFwshift)
  
  if (length(UFwshift)>0)
    pts=rbind(pts, UFwshift)
  if (length(DFwshift)>0)
    pts=rbind(pts, DFwshift)
  if (length(LFwshift)>0)
    pts=rbind(pts, LFwshift)
  if (length(RFwshift)>0)
    pts=rbind(pts, RFwshift)
  if (length(UBshift)>0)
    pts=rbind(pts, UBshift)
  if (length(DBshift)>0)
    pts=rbind(pts, DBshift)
  if (length(LBshift)>0)
    pts=rbind(pts, LBshift)
  if (length(RBshift)>0)
    pts=rbind(pts, RBshift)
  if (length(LUshift)>0)
    pts=rbind(pts, LUshift)
  if (length(RUshift)>0)
    pts=rbind(pts, RUshift)
  if (length(LDshift)>0)
    pts=rbind(pts, LDshift)
  if (length(RDshift)>0)
    pts=rbind(pts, RDshift)
  
  pts
}





get_max<-function(r,mat,av){
  
 
  num_mat=length(mat)
  Max=matrix()

  if (num_mat>4){
    num_l=dim(mat)[1] 
    num_c=dim(mat)[2]
    mat_tempo=matrix(0,num_l,5)
    index=0
    mat_tempo[,1]=mat[,1]
    mat_tempo[,2]=mat[,2]
    mat_tempo[,3]=mat[,3]
    mat_tempo[,4]=mat[,4]
    
    Max=matrix(0,1,4)

    for (i in 1:num_l){
      for(y in 1:num_l){
        Z=(mat_tempo[i,1]-mat_tempo[y,1])^2+(mat_tempo[i,2]-mat_tempo[y,2])^2+(mat_tempo[i,3]-mat_tempo[y,3])^2
        mat_tempo[y,5]=sqrt(Z)
      }
      
      MAX=TRUE
      count=0
      for (j in 1:num_l){
        
        if (mat_tempo[j,5]<r){
          count=count+1
        }
        if (mat_tempo[j,5]<r && mat_tempo[j,4]>mat_tempo[i,4] && i!=j){
          MAX=FALSE
          break
        }
        if (mat_tempo[j,5]<r && mat_tempo[j,4]==mat_tempo[i,4] && i!=j && j<i){
          MAX=FALSE
          break
        }
      }
      
      if (count<(av+1) | count<3){
        MAX==FALSE
        }
      
      if (MAX==TRUE && mat_tempo[i,4]>0){
        if (index==0){
          Max=mat_tempo[i,1:4]
          index=index+1}else{
            Max=rbind(Max,mat_tempo[i,1:4])
          }
      }
    }
    if (index==0){
      Max=matrix()
    }
  }
  return(Max) 
}





get_tp_values<-function(r,Max_infot, mat){
  
  num_max=length(Max_infot)
  
  if (num_max==1){
    Max_with_TP=matrix()
    
  }
  
  if (num_max==4){
    Max_with_TP=c(0,0,0,0,0,0)
    Max_with_TP[1]=Max_infot[1] 
    Max_with_TP[2]=Max_infot[2] 
    Max_with_TP[3]=Max_infot[3]
    Max_with_TP[4]=Max_infot[4] 
    Max_with_TP[5]=0 
    Max_with_TP[6]=Max_infot[4] 
    
  }
  
  if (num_max>4){
    
    num_l=dim(mat)[1]
    num_c=dim(mat)[2]
    num_max=dim(Max_infot)[1]
    Max_info=matrix(0,num_max,5)
    Max_info[,1:4]=Max_infot
 
    for (h in 1:num_max){
      go=0
      for (g in 1:num_l){
        if (go==0){
          if (mat[g,1]==Max_info[h,1] && mat[g,2]==Max_info[h,2] && mat[g,3]==Max_info[h,3]){
            Max_info[h,5]=g
            go=1
          }
        }
      }
    }
    
    Max_with_TP=matrix()
    Max_with_TP=matrix(0,num_max,6)
    coord=matrix()
    coord=cbind(mat[,1],mat[,2],mat[,3])
    dist_mat=as.matrix(dist(coord))
    
    
    for (t in 1:num_max){
      
      
      mat_to_flood=matrix()
      labels_in_mat=1:num_l
      mat_to_flood=cbind(mat,labels_in_mat)
      num_max_sup=t-1
      saddle_point=0
      saddle_point_previous=0
      saddle_point_new=0
      flooding_coeff=Max_info[t,4]/2 
      flooding_previous=0
      flooding=0
      ##
      
      if (t==1){ # The Everest case
        Max_with_TP[t,1]=Max_info[t,1] 
        Max_with_TP[t,2]=Max_info[t,2] 
        Max_with_TP[t,3]=Max_info[t,3] 
        Max_with_TP[t,4]=Max_info[t,4] 
        Max_with_TP[t,5]=Max_info[t,5] 
        Max_with_TP[t,6]=Max_info[t,3] 
      }else{
        
        finish=FALSE
        first_thing_first=0
        once=0
        
        while (finish==FALSE){

          mat_to_flood[,4]=mat[,4]-flooding
          C=which(mat_to_flood[,4]>0)
          
          if (length(C)>0){
            G=graph.adjacency(dist_mat[C,C]<r) 
            lab=clusters(G,"weak")
            init=num_l+1
            fin=2*num_l
            labels=init:fin
            labels[C]=lab$membership
            
            CONNECT=FALSE
            
            for (i in 1:num_max_sup){ 
              if (labels[Max_info[t,5]]==labels[Max_info[i,5]]){
                first_thing_first=1
                CONNECT=TRUE
                break
              }else{
                CONNECT=FALSE
              }
            }
            
            if (CONNECT==FALSE && first_thing_first==0){ 
              saddle_point=0 
              first_thing_first=1
              finish=TRUE
              break
            }
            if (CONNECT==FALSE && first_thing_first!=0){ 
              flooding=flooding_previous-flooding_coeff 
            }
            if (CONNECT==TRUE){ 
              flooding=flooding_previous+flooding_coeff
            }
            
            
            flooding_previous=flooding
            flooding_coeff=flooding_coeff/2
            
            if (flooding_coeff<0.5){ 
                if (CONNECT==TRUE){
                  while (CONNECT==TRUE){
                    flooding=flooding+1
                    mat_to_flood[,4]=mat[,4]-flooding
                    C=which(mat_to_flood[,4]>0)
                    if (length(C)>1){
                      G=graph.adjacency(dist_mat[C,C]<r)
                      lab=clusters(G,"weak")
                      init=num_l+1
                      fin=2*num_l
                      labels=init:fin
                      labels[C]=lab$membership
                      clor=0
                      for (i in 1:num_max_sup){
                        if (clor==0){
                          if (labels[Max_info[t,5]]==labels[Max_info[i,5]]){
                            CONNECT=TRUE
                            clor=1
                            break
                          }else{
                            CONNECT=FALSE
                            break
                          }}
                        
                      }}else{
                        CONNECT=FALSE
                        break
                      }
                  }
                  saddle_point=flooding
                  finish=TRUE
                  break
                }
                
                if (CONNECT== FALSE){
                  saddle_point=flooding
                  finish=TRUE
                  break
                }
                
            } 
            
          }
          
        } 
        Max_with_TP[t,1]=Max_info[t,1] 
        Max_with_TP[t,2]=Max_info[t,2] 
        Max_with_TP[t,3]=Max_info[t,3] 
        Max_with_TP[t,4]=Max_info[t,4] 
        Max_with_TP[t,5]=Max_info[t,5] 
        Max_with_TP[t,6]=Max_info[t,4]-saddle_point 
      }
      } 
    } 

  return(Max_with_TP) 
}



above_th<-function(Max_TP_values,th){
  
  num_max=length(Max_TP_values)
  
  if (num_max==1){
    Max_above_th=matrix()
  }
  if (num_max==6){
    if (Max_TP_values[6]>th){
      Max_above_th=Max_TP_values
    }else{
      Max_above_th=matrix()
    }
  }
  if (num_max>6){
    num_l=dim(Max_TP_values)[1]
    index=0
    Max_above_th=matrix(0,1,6)
    if (num_l>1){
      for (i in 1:num_l){
        if (Max_TP_values[i,6]>th){
          index=index+1
          if (index==1){
            Max_above_th=Max_TP_values[i,]}else{
              Max_above_th=rbind(Max_above_th,Max_TP_values[i,])
            }
        }
      }
      if (index==0){
        Max_above_th=matrix()
      }}
  }
  return(Max_above_th)
}


get_clusters<-function(r, Max_above_th,  mat_croped, mat){
  
  num_l_croped=dim(mat_croped)[1]
  num_c_croped=dim(mat_croped)[2]
  num_l=dim(mat)[1]
  num_c=dim(mat)[2]
  num_max=length(Max_above_th)
  
  if (num_max==1){
    mat_with_clusters=matrix(0,num_l,4)
    mat_with_clusters[,1:3]=mat[,1:3]
    labels_in_mat=1:num_l
    mat_with_clusters=cbind(mat_with_clusters,labels_in_mat)
    
  }
  if (num_max==6){
    
    mat_with_clusters=matrix()
    mat_with_clusters=matrix(0,num_l,4)
    mat_with_clusters[,1:3]=mat[,1:3]
    mat_with_clusters_croped=matrix()
    mat_with_clusters_croped=matrix(0,num_l_croped,4)
    mat_with_clusters_croped[,1:3]=mat_croped[,1:3]
    coord=cbind(mat_croped[,1],mat_croped[,2],mat_croped[,3])
    dist_mat=as.matrix(dist(coord))
    
    
    for (g in 1:num_l_croped){
      if (mat_croped[g,1]==Max_above_th[1] && mat_croped[g,2]==Max_above_th[2] && mat_croped[g,3]==Max_above_th[3]){
        max_index=g
      }
    }
    
    C=which(mat_croped[,4]>0)
    G=graph.adjacency(dist_mat[C,C]<r)
    lab=clusters(G,"weak")
    init=num_l_croped+1
    fin=2*num_l_croped
    labels=init:fin
    labels[C]=lab$membership
    for (rr in 1:num_l_croped){
      if (labels[max_index]==labels[rr]){
        mat_with_clusters_croped[rr,4]=1
      }
    }
    
    
    for (k in 1:num_l_croped){
      yep=0
      for (l in 1:num_l){
        if (yep==0){
          if (mat_croped[k,1]==mat[l,1] && mat_croped[k,2]==mat[l,2] && mat_croped[k,3]==mat[l,3]){
            mat_with_clusters[l,4]=mat_with_clusters_croped[k,4]
            yep=1
          } 
        }
      }
    }
    
    index=1
    for (i in 1:num_l){
      if (mat_with_clusters[i,4]==0){
        index=index+1
        mat_with_clusters[i,4]=index
      }
    }

  }

  if (num_max>6){
    
    num_max=dim(Max_above_th)[1]
    Max_info=matrix(0,num_max,5)
    Max_info[,1:4]=Max_above_th[,1:4]

   
    for (h in 1:num_max){
      go=0
      for (g in 1:num_l_croped){
        if (go==0){
          if (mat_croped[g,1]==Max_info[h,1] && mat_croped[g,2]==Max_info[h,2] && mat_croped[g,3]==Max_info[h,3]){
            Max_info[h,5]=g
            go=1
          }
        }
      }
    }
    
    mat_with_clusters=matrix()
    mat_with_clusters=matrix(0,num_l,4)
    mat_with_clusters[,1:3]=mat[,1:3]
    mat_with_clusters_croped=matrix()
    mat_with_clusters_croped=matrix(0,num_l_croped,4)
    mat_with_clusters_croped[,1:3]=mat_croped[,1:3]
    
    coord=cbind(mat_croped[,1],mat_croped[,2],mat_croped[,3])
    dist_mat=matrix()
    dist_mat=as.matrix(dist(coord))
    
    
    for (t in 1:num_max){
      t=num_max-t+1
      mat_to_flood=matrix()
      mat_to_flood=mat_croped[,1:4]
      flooding_previous=0
      flooding=Max_above_th[t,4]-Max_above_th[t,6]
      
      mat_to_flood[,4]=mat_croped[,4]-flooding
      C=which(mat_to_flood[,4]>0)
        
        if (length(C)>0){
          G=graph.adjacency(dist_mat[C,C]<r)
          lab=clusters(G,"weak")
          init=num_l_croped+1
          fin=2*num_l_croped
          labels=init:fin
          labels[C]=lab$membership
          CONNECT=FALSE
          for (rr in 1:num_l_croped){
            if (labels[Max_info[t,5]]==labels[rr]){
              if (mat_with_clusters_croped[rr,4]>0){
                CONNECT=TRUE
                cc=mat_with_clusters_croped[rr,4]
              }
            }}
          if (CONNECT==TRUE){
            
            n=which(Max_info[,5]==cc)
            
            mat_to_flood=matrix()
            mat_to_flood=mat_croped[,1:4]
            flooding_previous=0
            flooding=Max_above_th[n,4]-Max_above_th[n,6]# the TP value
            
            mat_to_flood[,4]=mat_croped[,4]-flooding
            C=which(mat_to_flood[,4]>0)
            
            if (length(C)>0){
              G=graph.adjacency(dist_mat[C,C]<r)
              lab=clusters(G,"weak")
              init=num_l_croped+1
              fin=2*num_l_croped
              labels=init:fin
              labels[C]=lab$membership}
          }
          
          for (rr in 1:num_l_croped){
            if (labels[Max_info[t,5]]==labels[rr]){
              mat_with_clusters_croped[rr,4]=Max_info[t,5]
            }}}}
      
      
      
    for (w in 1:num_max){
      count=0
      for (k in 1:num_l_croped){
        if  (mat_with_clusters_croped[k,4]==Max_info[w,5]){
          count=count+1
        }
      }
      if (count<6){
        for (k in 1:num_l_croped){
          if  (mat_with_clusters_croped[k,4]==Max_info[w,5]){
            mat_with_clusters_croped[k,4]=0
          }
        }
      }
    }
      
    
    for (k in 1:num_l_croped){
      yep=0
      for (l in 1:num_l){
        if (yep==0){
          if (mat_croped[k,1]==mat[l,1] && mat_croped[k,2]==mat[l,2] && mat_croped[k,3]==mat[l,3]){
            mat_with_clusters[l,4]=mat_with_clusters_croped[k,4]
            yep=1
          } }}}

    
    index=max(Max_info[,5])
    for (i in 1:num_l){
      if (mat_with_clusters[i,4]==0){
        index=index+1
        mat_with_clusters[i,4]=index
      }
    }
    
  }
  return(mat_with_clusters) 
}




Kclust <- function(pts, sdsx=0, sdsy=0, sdsz=0, xlim, ylim, zlim, paramc, paramd, psd=NULL, minsd=NULL, maxsd=NULL, useplabel=TRUE, alpha=NULL, pb=.5, rseq=seq(10, 150, by=5), thseq=seq(0, 100, by=2.5), score=FALSE, rlabel=FALSE, report=TRUE){
  ##browser()
  getting_in=0
  getting_in<<-getting_in
  N = dim(pts)[1]    
  if (N==1){
    rs=c()
    ths=c()
    for (r in rseq){
      for (th in thseq){
        rs=c(rs, r); ths=c(ths,th)
      }
    }
    labels=rep(1,length(rs)); dim(labels)=c(length(rs),1)
    return(list(scores=rep(0,length(rs)), scale=rs, thresh=ths, labels=labels))
  }
  tor=toroid(pts, xlim, ylim, zlim, range=max(rseq))
  D=as.matrix(dist(tor))
  D=D[1:N, 1:N]
  scores=c()
  retlabels=c()
  rs=c()
  ths=c()
  
  coord=cbind(pts[,1],pts[,2],pts[,3])
  dist_mat=matrix()
  num_l=dim(pts)[1]
  dist_mat=as.matrix(dist(coord))

  for (r in rseq){
    K=apply(D, 1, function(v){sum(v <= r)-1})
    L=((diff(xlim)+2*max(rseq))*(diff(ylim)+2*max(rseq))*(diff(zlim)+2*max(rseq))* K /(4/3 * pi * (dim(tor)[1]-1)))^(1/3)
    
    NEW_R=TRUE
    
    num_l=dim(pts)[1]
    X=c()
    Y=c()
    z=c()
    X=runif(num_l, xlim[1], xlim[2])
    Y=runif(num_l, ylim[1], ylim[2])
    z=runif(num_l, zlim[1], zlim[2])
    pts_sigma=matrix()
    pts_sigma=cbind(X,Y,z)
    


    coord=cbind(pts_sigma[,1],pts_sigma[,2],pts_sigma[,3])
    dist_mat=matrix()
    num_l=dim(pts_sigma)[1]
    dist_mat=as.matrix(dist(coord))
    Min=0
    min_vect=c()
    for (i in 1:num_l){
      first=0
      Min_tempo=0
      for (j in 1:num_l){
        if (i!=j){
          if (first==0){
            Min_tempo=dist_mat[i,j]
            first=1
          }else{
            Min_tempo=min(Min_tempo,dist_mat[i,j])
          }
        }
      }
      min_vect=cbind(min_vect,Min_tempo) 
      Min=Min+Min_tempo
    }
    r_scanning=Min/num_l


    av_point_per_nm=num_l/(3000*3000*600)
    av_point_in_sphere=av_point_per_nm*((4/3)*pi*r^3)

    

    tor_sigma=toroid(pts_sigma, xlim, ylim, zlim, max(rseq))
    D_sigma=as.matrix(dist(tor_sigma))
    D_sigma=D_sigma[1:N, 1:N]
    K_sigma=apply(D_sigma, 1, function(v){sum(v <= r)-1})
    L_sigma=((diff(xlim)+2*max(rseq))*(diff(ylim)+2*max(rseq))*(diff(zlim)+2*max(rseq))* K_sigma /(4/3 * pi * (dim(tor)[1]-1)))^(1/3)
    L_sigma_tempo=0
    L_mean=0
    
    for (j in 1:num_l){
      L_mean=L_mean+L_sigma[j]
    }
    L_mean=L_mean/num_l
    
    for (j in 1:num_l){
      L_sigma_tempo=L_sigma_tempo+(L_sigma[j]-L_mean)^2
    }
    L_sigma_tempo=L_sigma_tempo/num_l
    sigma_TBU=sqrt(L_sigma_tempo)

            
      if (NEW_R==TRUE){ 
        
        Lr=t(L)
        num_l=dim(pts)[1]
        mat=matrix(0,num_l,4)
        mat[,1]=pts[,1]
        mat[,2]=pts[,2]
        mat[,3]=pts[,3]

        mat[,4]=Lr-(L_mean+sigma_TBU)

        num_l=dim(mat)[1]
        for (i in 1:num_l){
          mat[i,4]= max(0,mat[i,4])
        }
        index=0
        
        for (i in 1:num_l){
          if (mat[i,4]!= 0){
            index=index+1
            
          }
        }
        
        
        mat_croped=matrix(0,index,4)
        if (index==0){
          mat_croped=matrix()  
        }
        
        index=0
        num_l=dim(mat)[1]
        for (i in 1:num_l){
          if (mat[i,4]!= 0){
            index=index+1
            mat_croped[index,]=mat[i,] 
          }
        }

        
        max_list=matrix()
        max_list=get_max(r,mat_croped,av_point_in_sphere) 
        
        
        NN=length(max_list)
        if (NN>4){
          num_max=dim(max_list)[1]

          max_list_ordered_coeff=order(max_list[,4],decreasing=TRUE)

          max_list_tempo=matrix(0,num_max,4)
          
          for (i in 1:num_max){
            for (j in 1:num_max){
              if (i==max_list_ordered_coeff[j]){
                max_list_tempo[j,]=max_list[i,]
              }
            }
          }
          max_list=c()
          max_list=max_list_tempo 
          num_max=dim(max_list)[1]
          max_list<<-max_list

        }
        if (NN==1){
          
          max_list=matrix()
        }

        max_tp_list=matrix()

        max_tp_list=get_tp_values(r_scanning,max_list, mat_croped)

        NEW_R=FALSE
      }
      

    NN=length(max_tp_list)
    if (NN>6){
      num_tp_max=dim(max_tp_list)[1]

      max_tp_list_ordered_coeff=order(max_tp_list[,6],decreasing=TRUE)

      max_tp_list_tempo=matrix(0,num_max,6)
      
      for (i in 1:num_tp_max){
        for (j in 1:num_tp_max){
          if (i==max_tp_list_ordered_coeff[j]){
            max_tp_list_tempo[j,]=max_tp_list[i,]
          }
        }
      }
      max_tp_list=c()
      max_tp_list=max_tp_list_tempo 
      num_tp_max=dim(max_tp_list)[1]
      max_tp_list<<-max_tp_list

    }
    if (NN==1){
      num_tp_max=0
      max_list=matrix()
    }
    if (NN==6){
      num_tp_max=1
   
      
    }
    
      if (num_tp_max>=1){    
      th_prec=-1
      for (p in 1:num_tp_max){
        q=(num_tp_max-p+1)
      th=max_tp_list[q,6] 
      if (th!=th_prec && th>0){
        th_prec=th
      max_above_tpth_list=matrix()
      max_above_tpth_list=above_th(max_tp_list,th)
      mat_with_clusters_labels=matrix()
      mat_with_clusters_labels=get_clusters(r_scanning, max_above_tpth_list, mat_croped, mat)
      labels=matrix(0,num_l,1)
      labels=mat_with_clusters_labels[,4]
      s=0
      getting_in=getting_in+1
      getting_in<<-getting_in
      
      if (score){
        if (getting_in>1){
          getting_in=getting_in+1
          getting_in<<-getting_in
          s=scorewprec3sd(labels=labels, pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
        
           if (s>=scores){
            scores=s
            scores<<-scores
            ths_best=floor(th)
            rs_best=r
            ths_best<<-ths_best
            rs_best<<-rs_best
          }                 
        
        }else{
          s=scorewprec3sd(labels=labels, pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
          getting_in=getting_in+1
          getting_in<<-getting_in
          scores=s
          scores<<-scores
          ths_best=floor(th)
          rs_best=r
          ths_best<<-ths_best
          rs_best<<-rs_best
        }
      
      }
       
      if (report){cat("Scale:", r, "Thr:", th, "Score: ", s, "\n")}
      if (rlabel){
        if (getting_in>1){
          s=scorewprec3sd(labels=labels, pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
          if (s>=scores){
            retlabels_best=labels
            retlabels_best<<-retlabels_best
          }                 
          
        }else{
          retlabels_best=labels
          retlabels_best<<-retlabels_best
 
        }
        

      }
      }
    }
  }
  if (num_tp_max==0){
    labels=matrix(0,num_l,1)
    for (u in 1:num_l){
      labels[u]=u
    }
    s=0
    getting_in=getting_in+1
    getting_in<<-getting_in
    
    if (score){
      if (getting_in>1){
        getting_in=getting_in+1
        getting_in<<-getting_in
        s=scorewprec3sd(labels=labels, pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
      
        if (s>=scores){
          scores=s
          scores<<-scores
          ths_best=floor(th)
          rs_best=r
          ths_best<<-ths_best
          rs_best<<-rs_best
        }                 
        
      }else{
        s=scorewprec3sd(labels=labels, pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
        getting_in=getting_in+1
        getting_in<<-getting_in
        scores=s
        
        scores<<-scores
        ths_best=floor(th)
        rs_best=r
        ths_best<<-ths_best
        rs_best<<-rs_best
      }
      
    }
    
    if (report){cat("Scale:", r, "Thr:", th, "Score: ", s, "\n")}
    if (rlabel){
      if (getting_in>1){
        s=scorewprec3sd(labels=labels, pts=pts, sdsx=sdsx, sdsy=sdsy, sdsz=sdsz, xlim=xlim, ylim=ylim, zlim=zlim, paramc=paramc, paramd=paramd, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
        if (s>=scores){
          retlabels_best=labels
          retlabels_best<<-retlabels_best
        }                 
        
      }else{
        retlabels_best=labels
        retlabels_best<<-retlabels_best
        
      }
      
      
    }
    
  }
  }
  

  scores_best=scores
  scores_best<<-scores_best
  list(scores=scores_best, scale=rs_best, thresh=ths_best, labels=retlabels_best)
}
########################################################################################################

writeRes <- function(res, rfile, labdir, bestonly=FALSE){
  scale=unique(res[["scale"]]); scale=scale[order(as.numeric(scale))]
  thresh = unique(res[["thresh"]]); thresh=thresh[order(as.numeric(thresh))]
  cat("0", scale, sep="\t", file=rfile); cat("\n", file=rfile, append=TRUE)
  for (line in thresh){
    scales=res[["scale"]][res[["thresh"]]==line]; o=order(scales); scales=scales[o]
    scores=res[["scores"]][res[["thresh"]]==line]; scores=scores[o]
    cat(line, "\t", sep="", file=rfile, append=TRUE); cat(scores, sep="\t", append=TRUE, file=rfile); cat("\n", file=rfile, append=TRUE)
  }
  dir.create(labdir)
  if (bestonly) is=which.max(res[["scores"]])
  else is=(1:dim(res[["labels"]])[1])
  for (i in is){
    f=file.path(paste(labdir, "/clusterscale", res[["scale"]][i], "\ thresh", res[["thresh"]][i], "labels.txt", sep=""))
    cat(res[["labels"]][], file=f, sep=","); cat("\n", file=f, append=TRUE)
  }
}


#Descriptors

nClusters <- function(labels){
  sum(table(labels)>1)
}

nClustersi <- function(labels){
  s=length(labels)
  Nclus=c()
  M=0
  for (i in 1:s){
    if (i==1){
      M=as.numeric(labels[i])
    }else{
      if (M<as.numeric(labels[i])){
        M=as.numeric(labels[i]) 
      }
    }
  }
  
  
  for (j in 1:M){ 
    count=0
    for (i in 1:s){
      if (as.numeric(labels[i])==j) {
        count=count+1
      }
    }
    if (count>1){
      Nclus=c(Nclus, j) 
    }
  }
  return(Nclus)
  
}



percentageInCluster <- function(labels){
  Nb=sum(table(labels)==1)
  (length(labels)-Nb)/length(labels) * 100
}

molsPerCluster <- function(labels){
  ta=table(labels); ta[ta>1]
}

nMolsPerCluster <- function(labels){
  length(labels)*percentageInCluster(labels)/(100*nClusters(labels))
}

histnMols <- function(labels){
  ta=table(labels)[table(labels)>1]; h=hist(ta, plot=FALSE)
  plot(h, xlab="Number of molecules", ylab="Number of clusters", main="")
}

molsPerClusteri <- function(labels,i){
  s=length(labels); 
  count=0
  n=as.numeric(i)
  
  for (h in (1:s)){
    a=as.numeric(labels[h])
    if (a==n){
      count=count+1
    }
  }
  return(count)
}


clusterRadii <- function(pts, labels){
  radii=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
      mean(c(sd(pts[v,1]), sd(pts[v,2]), sd(pts[v,3])))
    }
  })
  radii[radii>=0]
}

clusterRadiii <- function(pts, labels, i){
  id=i
  radii=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
      mean(c(sd(pts[v,1]), sd(pts[v,2]), sd(pts[v,3])))
    }
  })
  radii[radii>=0]
  for (g in (1:length(radii))){
    l=as.numeric(row.names(radii)[g])
    if(l==id){
      RADII=radii[g]
    }    
  }
  return(RADII)
}


convexHullAreas <- function(pts, labels){
  areas=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
      i<-chull(pts[v,1],pts[v,2],pts[v,3])
      areapl(as.matrix(pts[v[i],]))
    }
  })
  areas[areas>=0]
}

plabel <- function(labels, alpha, pb){
  cnt <-tapply(1:length(labels), labels, length)
  cl =cnt[cnt!=1]
  B = length(labels)-sum(cl)
  Bcont = B*log(pb)+(1-B)*log(1-pb)
  ## Green 2001 p.357, Scand J Statist 28
  partcont=0
  if (length(cl) >0)
    partcont=length(cl)*log(alpha)+lgamma(alpha)+sum(lgamma(cl))-lgamma(alpha+sum(cl))
  Bcont+partcont
}


scorewprec3sd <- function(labels, pts, sdsx, sdsy, sdsz, xlim, ylim, zlim, paramc, paramd, psd, minsd, maxsd, useplabel=TRUE, alpha=NULL, pb=.5){  
  s=sum(tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)>1) mcgaussprec3sd(pts[v,], sdsx[v], sdsy[v], sdsz[v], xlim, ylim, zlim, paramc, paramd, psd=psd, minsd=minsd, maxsd=maxsd)
    else -log(diff(xlim)*diff(ylim)*diff(zlim))      
  }))  
  prlab=0
  if (useplabel){
    if (is.null(alpha)){
      cnt <-tapply(1:length(labels), labels, length)
      n =sum(cnt[cnt!=1])
      alpha=20
    }
    prlab=plabel(labels, alpha, pb)
  }
  s+prlab
}

