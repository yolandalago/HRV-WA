
require("FMM")

HRV_WA<-function(data=data){
storeNSRDB<-list()
tableOut_m_pac<-c()
#for(i in 1:length(datosNSRDB)){#activate if more than one subject

	f<-data

	s<-seq(0,2*pi,length.out=(length(f)+1))#tiemposNSRDB[[i]]
	s<-s[-length(s)]
	m<-3
	i<-1#cancel if more than one subject
	for(j in 1:m){
		out<-computeA_heart(f,s,j)
		storeNSRDB[[m*(i-1)+j]]<-out
		tableOut_m<-cbind(rep(i,j),rep(j,j),1:j,rep(out[[1]]@M,j),
			coef(out[[1]])$wave,getFMMPeaks(out[[1]])$tpeakU,getFMMPeaks(out[[1]])$tpeakL,out[[1]]@R2)
		tableOut_m_pac<-rbind(tableOut_m_pac,tableOut_m)
	}
#}#activate if more than one subject
colnames(tableOut_m_pac)<-c("Pac","TotalComp","Comp","M","A","Alpha","Beta","Omega","t_U","t_L","R2_m")
fmm3NSRDB<-tableOut_m_pac[tableOut_m_pac[,"TotalComp"]==3,]



compNSRDB<-c()
#for(i in 1:length(datosNSRDB)){#activate if more than one subject
filasPac<-fmm3NSRDB[fmm3NSRDB[,"Pac"]==i,]

#Inizialitation
C1<-NA;C2<-NA

#Direct Wave
ordenFilas<-order(filasPac[,"R2_m"],decreasing=TRUE)
sigueC1<-TRUE
for(f in 1:nrow(filasPac)){
	fila<-ordenFilas[f]
	paramFila<-filasPac[fila,c("A","Beta","Omega","t_U","t_L","R2_m")]
	if(condC1_new(paramFila)[[1]] & sigueC1){
		C1<-fila
		sigueC1<-FALSE
		tuC1<-condC1(paramFila)[[2]]
	}
}

#Guided Wave

sigueC2<-TRUE
for(f in 1:nrow(filasPac)){
	fila<-ordenFilas[f]
	paramFila<-filasPac[fila,c("A","Beta","Omega","t_U","t_L","R2_m")]
	if(!is.na(C1)){
		if(fila!=C1 ){#& !is.na(C1)){#cambio
			if(condC2_new(paramFila,tuC1) & sigueC2){
				C2<-fila
				sigueC2<-FALSE
				tuC2<-paramFila[3]
			}
		}
	}else{
		if(condC2_new(paramFila,tuC1) & sigueC2){
			C2<-fila
			sigueC2<-FALSE
			tuC2<-paramFila[3]
		}
	}
}
compNSRDB<-rbind(compNSRDB,c(C1,C2))
#}




################################################################
#Outputs and plots

addFila<-c()
m<-3
nComp<-c()
r2_Fmm<-c()
r2_Cos<-c()
tu_Cos_NSRDB<-c()
tl_Cos_NSRDB<-c()
#for(i in 1:length(datosNSRDB)){

	

	#	Plots
	par(mfrow=c(1,2))#storeNSRDB[[m*(i-1)+j]]
	vPar<-cbind(storeNSRDB[[m*(i-1)+3]][[1]]@A,storeNSRDB[[m*(i-1)+3]][[1]]@alpha,storeNSRDB[[m*(i-1)+3]][[1]]@beta,storeNSRDB[[m*(i-1)+3]][[1]]@omega)
	
	if(!is.na(compNSRDB[i,1]) & !is.na(compNSRDB[i,2])){
		a<-genPar(vPar[compNSRDB[i,1],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
		b<-genPar(vPar[compNSRDB[i,2],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
		Madj<-optim(par=0.5,fn=ff,method="Brent",lower=0.1,upper=2.9)$par
		ajuste<-Madj+a+b#storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+(a-a[1])+(b-b[1])
	}else{
		if(is.na(compNSRDB[i,1])){
			a<-0
			b<-genPar(vPar[compNSRDB[i,2],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
			Madj<-optim(par=0.5,fn=ff,method="Brent",lower=0.1,upper=2.9)$par
			ajuste<-Madj+a+b#storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+(b-b[1])
		}else{
			b<-0
			a<-genPar(vPar[compNSRDB[i,1],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
			Madj<-optim(par=0.5,fn=ff,method="Brent",lower=0.1,upper=2.9)$par
			ajuste<-Madj+a+b#storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+(a-a[1])
		}
	}
	numFmm<-sum((ajuste-f)^2)
	denFmm<-sum((f-mean(f))^2)
	r2_Fmm[i]<-1-numFmm/denFmm
	nComp[i]<-sum(!is.na(compNSRDB[i,]),na.rm=TRUE)
	plot(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,storeNSRDB[[m*(i-1)+3]][[1]]@data,col=1,type="p",ylab="",xlab="",xaxt="n",
		main=paste("Pac_",i," Data and Prediction",sep=""))
	lines(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,ajuste,col=4,lwd=3)
	axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c("0","6","12","18","24"))
	legend("topright",c("RR Data","HRV Prediction"),lwd=c(NA,3),col=c(1,4),lty=c(NA,1),pch=c(1,NA))

	#COMPONENTES
	comp1<-NA;comp2<-NA;comp3<-NA
	if(!is.na(compNSRDB[i,1]) & !is.na(compNSRDB[i,2])){
		compA<-genPar(vPar[compNSRDB[i,1],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
		compA<-storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+compA-compA[1]
		compB<-genPar(vPar[compNSRDB[i,2],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
		compB<-storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+compB-compB[1]
		comp1<-compNSRDB[i,1]
		comp2<-compNSRDB[i,2]
	}else{
		if(is.na(compNSRDB[i,1])){
			compA<-genPar(vPar[compNSRDB[i,2],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
			compA<-storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+compA-compA[1]
			comp1<-compNSRDB[i,2]
			comp2<-NA #faltaba y no pintaba bien cosinor, por eso hago las carpetas bis
		}else{
			compA<-genPar(vPar[compNSRDB[i,1],],storeNSRDB[[m*(i-1)+3]][[1]]@timePoints)
			compA<-storeNSRDB[[m*(i-1)+3]][[1]]@data[1]+compA-compA[1]
			comp1<-compNSRDB[i,1]
			comp2<-NA #faltaba y no pintaba bien cosinor, por eso hago las carpetas bis
		}
	}

	plot(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,compA,col=3,lwd=3,type="l",ylab="",xlab="",xaxt="n",
		main=paste("Pac_",i," Components",sep=""),ylim=c(min(c(compA,compB),na.rm=TRUE),max(c(compA,compB),na.rm=TRUE)))
	if(!is.na(comp2))lines(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,compB,col="orange",lwd=3)

	axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c("0","6","12","18","24"))
	legend("topright",c("Direct","Guided"),lwd=c(3,2),col=c(3,"orange"),lty=c(1,1))


	#	Tabla
	
	subTabla<-fmm3NSRDB[fmm3NSRDB[,"Pac"]==i,]
	for(j in 1:ncol(compNSRDB)){
		if(!is.na(compNSRDB[i,j])){
			addFila<-rbind(addFila,c(unlist(subTabla[compNSRDB[i,j],3:ncol(subTabla)])))
		}else{
			addFila<-rbind(addFila,c(rep(NA,ncol(subTabla)-2)))
		}
	}

#}
	return(addFila)
}


#################	Functions ###########

condC1_new<-function(params){#"A","Beta","Omega","t_U","t_L","R2_m"
	a<-params[1]
	b<-params[2]
	o<-params[3]
	tu<-params[4]
	tl<-params[5]
	r2<-params[6]
	thr1<-pi/8
	if(o<0.2){
		if(b<pi/4 | b>7*pi/4)tu<-tl
	}else{
		if(o>=0.2 & o<0.4){
			if(b<pi/2 | b>3*pi/2)tu<-tl
		}else{
			if(b<3*pi/4 | b>5*pi/4)tu<-tl
		}
	}
	value<-r2>0.03 & a>0.075 & o>0.0275 & ((tu<pi & tu>thr1) | (tu>pi & tu<(2*pi-thr1)))
	return(list(value,tu))
}
#params=paramFila
#tu=tuC1
condC2_new<-function(params,tu){#"A","Beta","Omega","t_U","t_L","R2_m"
	a<-params[1]
	b<-params[2]
	o<-params[3]
	tu2<-params[4]
	tl2<-params[5]
	r2<-params[6]
	tu1<-tu
	thr1<-pi/8

	if(o<0.2){
		if(b<pi/4 | b>7*pi/4)tu2<-tl2
	}else{
		if(o>=0.2 & o<0.4){
			if(b<pi/2 | b>3*pi/2)tu2<-tl2
		}else{
			if(b<3*pi/4 | b>5*pi/4)tu2<-tl2
		}
	}
	if(!is.na(tu1)){
		if(tu1<pi){
			return( r2>0.03 & a>0.05 & o>0.0275 & (tu2>(tu1+thr1)%%(2*pi) &
										 ((tu2>pi & o>0.05 & (b>pi/2 & b<3*pi/2))| tu2>5*pi/4) &
										 tu2<(2*pi-thr1)) )
		}else{
			return( r2>0.03 & a>0.05 & o>0.0275 & (tu2<(tu1)%%(2*pi) & 
                                                             ((tu2<(3*pi/2) & o>0.05 & (b>pi/2 & b<3*pi/2))) & tu2>(thr1)) )
		}
	}else{
		return( r2>0.03 & a>0.05 & o>0.0275 & ((tu2<(3*pi/2) & tu2>thr1) | (((tu2>pi & o>0.05 & (b>pi/2 & b<3*pi/2))| tu2>5*pi/4)  & tu2<(2*pi-thr1)) ))
	}
}



ff<-function(M){
			x1<-M[1]
			return(sum(((x1+a+b)-storeNSRDB[[m*(i-1)+3]][[1]]@data)^2))
}


genPar<-function(p,t){
		a<-p[1];al<-p[2];b<-p[3];o<-p[4]
		return(a*cos(b+2*atan2(o*sin((t-al)/2),cos((t-al)/2))))
	}


computeA_heart<-function(f,s,m){
	nObs<-length(f)
	fit_m<-fitFMM(vData=f,timePoints=s,nback=m)
	return(list(fit_m))
}


 condC1<-function(params){#"A","Beta","Omega","t_U","t_L","R2_m"
a<-params[1]
b<-params[2]
o<-params[3]
tu<-params[4]
tl<-params[5]
r2<-params[6]
thr1<-pi/8
if(b<pi/4 | b>7*pi/4)tu<-tl
value<-r2>0.03 & a>0.1 & o>0.0275 & ((tu<pi & tu>thr1) | (tu>3*pi/2 & tu<(2*pi-thr1)))
return(list(value,tu))
}
