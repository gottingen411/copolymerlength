computedist = function(nA,nB,thetaa,thetab,nsamples,flag) {
  na=nA
  nb=nB
  n=na+nb
  a_set=as.integer(rep(0,na))
  b_set=as.integer(rep(1,nb))
  total_set=c(a_set,b_set)
  thetas=c(thetaa,thetab)*pi/180
  
  #bonds=c(1,1)
  n_samples=nsamples
  if(!flag){
  bonds=sin(thetas/2)*2
  lengths=rep(0,n_samples)
  #x=rep(0,n_samples)
  #y=rep(0,n_samples)
  for (i in 1:n_samples){
    loc_origin=matrix(c(0,0),nrow=2)
    directions=runif(n,min=-1,max=1)
    #directions=sample(c(-1,1),n,replace=TRUE)
    polymer=sample(total_set,n,replace=FALSE)
    
    
    for (j in 1:n){
      d=directions[j]*pi
      if (!polymer[j]){
        loc_cord=matrix(c(bonds[1]*cos(d),bonds[1]*sin(d)),nrow=2)
        cord=loc_cord+loc_origin
      } else{
        loc_cord=matrix(c(bonds[2]*cos(d),bonds[2]*sin(d)),nrow=2)
        cord=loc_cord+loc_origin
        
      }
      loc_origin=cord
      #x[j]=loc_origin[1]
      #y[j]=loc_origin[2]
    }
    lengths[i]=norm(cord,type='f')
    #x[i]=cord[1]
    #y[i]=cord[2]
  }
  
  }else{
    
    
    lengths=rep(0,n_samples)
    directions=as.integer(c(-1,1))
    for (i in 1:n_samples){
      loc_origin=matrix(c(0,0),nrow=2)
      rot_mat=diag(2)
      polymer=sample(total_set,n,replace=FALSE)
      dir_i=sample(directions,n,replace=TRUE)
      
      for (j in 1:n){
        d=dir_i[j]
        if (!polymer[j]){
          loc_cord=matrix(2*c(cos((pi-thetas[1])/2)**2,d*sin(thetas[1]/2)*sin((pi-thetas[1])/2)),nrow=2)
          cord=rot_mat %*% loc_cord+loc_origin
          rot_mat=rbind(c(cos(pi-thetas[1]),-sin(d*(pi-thetas[1]))),
                        c(sin(d*(pi-thetas[1])),cos(pi-thetas[1]))) %*% rot_mat
        } else{
          loc_cord=matrix(2*c(cos((pi-thetas[2])/2)**2,d*sin(thetas[2]/2)*sin((pi-thetas[2])/2)),nrow=2)
          cord=rot_mat %*% loc_cord+loc_origin
          rot_mat=rbind(c(cos(pi-thetas[2]),-sin(d*(pi-thetas[2]))),
                        c(sin(d*(pi-thetas[2])),cos(pi-thetas[2]))) %*% rot_mat
        }
        loc_origin=cord
      }
      lengths[i]=norm(loc_origin,type='f')
    }
  }
  return(lengths)
}
