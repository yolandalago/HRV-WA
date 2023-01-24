

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

leoNsrdb<-fmm3NSRDB
add2Nsrdb<-c()
#uniPacNsrdb<-unique(leoNsrdb[,"Pac"])#activate if more than one subject
#for(i in 1:length(uniPacNsrdb)){#activate if more than one subject
	#filas<-which(uniPacNsrdb[i]==leoNsrdb[,"Pac"])#activate if more than one subject
	#mSmall<-(leoNsrdb[filas,])[order(leoNsrdb[filas,"R2_m"],decreasing=TRUE),]#activate if more than one subject
	mSmall<-fmm3NSRDB
	asig<-0
	noche<-FALSE
	k<-1
	while(asig<2 & k<=3){
		if(mSmall[k,"Omega"]>0.035 & mSmall[k,"R2_m"]>0.03 & mSmall[k,"A"]>0.05 ){
			if(mSmall[k,"Alpha"]>pi & noche==FALSE ){#noche no esignada o dia
				asig<-asig+1
				noche<-TRUE
				add2Nsrdb<-rbind(add2Nsrdb,mSmall[k,])
				k<-k+1
			}else{#ya se ha asignado una noche, o bien es dia
				if(mSmall[k,"Alpha"]<(pi-pi/8)){
					asig<-asig+1
					add2Nsrdb<-rbind(add2Nsrdb,mSmall[k,])
					k<-k+1
				}else{#si no es dia sigue sumando
					k<-k+1
				}
			}
		}else{
			k<-k+1
		}#si es ruido sigue
	}
#}



################################################################
#Outputs and plots

m<-3
addR2FmmNSRDB<-c()
addR2CosNSRDB<-c()

tu_Cos_NSRDB<-c()
tl_Cos_NSRDB<-c()
r2_Fmm<-c()
nComp<-c()
i=1
compNSRDB<-matrix(c(add2Nsrdb[,3]),1,2)
x11()
	

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
	numFmm<-sum((ajuste-data)^2)
	denFmm<-sum((data-mean(data))^2)
	r2_Fmm[i]<-1-numFmm/denFmm
	nComp[i]<-sum(!is.na(compNSRDB[i,]),na.rm=TRUE)
	plot(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,storeNSRDB[[m*(i-1)+3]][[1]]@data,col=1,type="p",ylab="",xlab="",xaxt="n",
		main=paste("Pac_",i,"_Fitted R2: ",round(r2_Fmm[i],3),sep=""))
	lines(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,ajuste,col=3,lwd=3)
	axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c("0","6","12","18","24"))

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

	plot(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,compA,col="grey",lwd=3,type="l",ylab="",xlab="",xaxt="n",
		main=paste("Pac_",i,"_Components",sep=""),ylim=c(min(c(compA,compB),na.rm=TRUE),max(c(compA,compB),na.rm=TRUE)))
	if(!is.na(comp2))lines(storeNSRDB[[m*(i-1)+3]][[1]]@timePoints,compB,col="lightblue",lwd=3)

	axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c("0","6","12","18","24"))
	legend("topright",c("Pricipal","Secondary"),lwd=rep(3,2),col=c("grey","lightblue"),lty=rep(1,2))

return(add2Nsrdb[,-c(1:3)])
}
#################	Functions ###########


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


