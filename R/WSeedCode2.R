#' List of functions

require('magrittr')
eps<-10^(-5);
is.odd <- function(x) x %% 2 != 0
globalVariables(c("fn.rep","num.fold"))


#' freq.less function
#'
#' This function finds the count the x-sample observes is greater than the one from y-sample based on ranks. (description)
#' @param fx,fy This is where any argument goes
#' @details This is where any details go (delete if unnecessary)
#' @return c(less.count, more.count) (return value here)
#' @note This is where any additional notes go (delete if unnecessary)
#' @references This is where any references go (delete if unnecessary)
#' @export
#' @examples
#' freq.less(fx=c(1,2,4,9,0,0,NA),fy=c(1,4,9,NA))

freq.less<-function(fx,fy)
{
  if(all(is.na(fx))|all(is.na(fy))) return(c(NA,NA));

  n.fx<-length(fx);n.fy<-length(fy); n.all<-n.fx+n.fy;

  fx[abs(fx)<eps]<-0;fy[abs(fy)<eps]<-0;

  less.count<-sum(outer(fx,fy,'<')+outer(fx,fy,'==')*.5,na.rm=TRUE)
  more.count<-sum(outer(fx,fy,'>')+outer(fx,fy,'==')*.5,na.rm=TRUE)

  return(c(less.count, more.count))
}


#' multi.prob function
#'
#' This function creates a vector and then adds data to the vector (description)
#' @param fsam (Arguments here)
#' @details (details here; delete if unnecessary)
#' @return prob.vec (return value)
#' @note (notes here; delete if unnecessary)
#' @references (references here; delete if unnecessary)
#' @export
#' @examples
#' multi.prob(fsam=NULL)

multi.prob<-function(fsam)
{
  n.col<-length(fsam);
  prob.vec<-NULL;
  for(fi in 2:n.col)
  {temp.f<-freq.less(fsam[[fi-1]],fsam[[fi]])
  prob.vec<-c(prob.vec,temp.f);}
  return(prob.vec)
}


#' simu.cases function
#'
#' This function... (description)
#' @param N.boot.pow,n.sub,n.sam,mean.sub,sd.sub (argument here)
#' @details (details here; DIU)
#' @return list.sam (return value)
#' @note (notes here; DIU)
#' @references (references here; DIU)
#' @importFrom stats "rnorm"
#' @export
#' @examples
#' simu.cases(N.boot.pow=10,n.sub=10,n.sam=5,mean.sub=c(1,2,1,2,1),sd.sub=rep(0.2,n.sam=5))

simu.cases<-function(N.boot.pow,n.sub=10,n.sam=5,mean.sub=c(1,2,1,2,1),sd.sub=rep(0.2,n.sam)) ### pure simulations
{
  sam<-NULL;
  for (i in 1:n.sam) sam<-cbind(sam,rnorm(n.sub*N.boot.pow,mean.sub[i],sd.sub[i]));
  list.sam <- lapply(split(sam,rep(1:N.boot.pow,rep(n.sub,N.boot.pow))),matrix,nrow=n.sub)
  return(list.sam);
}


#' simu.ustat.pattern function
#'
#' This function will simulate the use of the pattern (description here)
#' @param mean.prob.vec,effn.subs,n.rep (argument here)
#' @details (details here; DIU)
#' @return simu.tab (return value)
#' @note (notes here; DIU)
#' @references (references here; DIU)
#' @importFrom stats "qnorm" "rnorm"
#' @export
#' @examples
#' simu.ustat.pattern(mean.prob=c(.82,.1,.6))

simu.ustat.pattern<-function(mean.prob.vec, effn.subs=rep(5,length(mean.prob.vec)),n.rep=10^2)
{
  mean.prob.vec<-replace(mean.prob.vec, mean.prob.vec<eps, eps);mean.prob.vec<-replace(mean.prob.vec, mean.prob.vec>1-eps, 1-eps);
  n.col=length(mean.prob.vec);
  h<-cumsum(2^.5*qnorm(mean.prob.vec[1]));

  simu.norm.sam1.all<-matrix(rnorm(effn.subs[1]*n.rep, mean=0),nrow=n.rep);
  simu.norm.sam2.all<-matrix(rnorm(effn.subs[2]*n.rep, mean=h),nrow=n.rep);
  simu.tab<-lapply(c(1:n.rep), function(x) multi.prob(list(simu.norm.sam1.all[x,],simu.norm.sam2.all[x,])))
  return(simu.tab)

}


#' chi.stat function
#'
#' This function... (description here)
#'
#' This function adds 50 to the user input.
#' @param ftab (argument here)
#' @details (details here; DIU)
#' @return chi.val (return value)
#' @note (notes here; DIU)
#' @references (references here; DIU)
#' @export
#' @examples
#' chi.stat(ftab=rbind(c(20,10,20),c(15,15,20)))

chi.stat<-function(ftab)
{
  tot<-sum(ftab);
  expv<-outer(rowSums(ftab)/tot,  colSums(ftab)/tot, '*')*tot;
  signal<-(colMeans(expv)<eps)*(1:dim(ftab)[2]);
  indx<-setdiff(signal,0);
  ftemp<-((ftab-expv)^2/expv)
  chi.val<-ifelse(length(indx)==0, sum(ftemp), sum(ftemp[,-indx]));

  return(chi.val);
}


#' gen.decision function
#'
#' This function will test to see if two functions are the same or not the same (description)
#' @param est.prob,effn.subsam1,effn.subsam2,fn.rep,alpha (argument here)
#' @details (details here; DIU)
#' @return result (return value)
#' @note (notes here; DIU)
#' @references (references here; DIU)
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' gen.decision(rbind(c(20,5,10,15,20,5),c(15,10,15,10,20,5)),effn.subsam1=rep(5,4),effn.subsam2=rep(5,4),fn.rep=10^3,alpha=.05)

gen.decision<-function(est.prob, effn.subsam1=rep(5,8),effn.subsam2=rep(5,8),fn.rep=10^3,alpha=.05)

{
  n.col=length(effn.subsam1);
  base.n=rbind(rep(effn.subsam1[-n.col]*effn.subsam1[-1],rep(2,n.col-1)),
               rep(effn.subsam2[-n.col]*effn.subsam2[-1],rep(2,n.col-1)));

  mean.prob<-colSums(est.prob)/colSums(base.n);
  id.na<-((mean.prob%>%is.na)*c(1:length(mean.prob)))%>%setdiff(0);
  if(length(id.na)==0) id.na=NA;
  id.keep<-setdiff(c(1:length(mean.prob)),id.na)


  mean.prob.simu<-matrix(mean.prob[id.keep],ncol=2,byrow=TRUE);
  id.effb<-unique(round(id.keep/2+eps));



  effn.sub.simu1<-cbind(effn.subsam1[-n.col],effn.subsam1[-1])[id.effb,];
  effn.sub.simu2<-cbind(effn.subsam2[-n.col],effn.subsam2[-1])[id.effb,];

  num.fold<-dim(mean.prob.simu)[1];
  if(num.fold==1) { simu.sam1= simu.ustat.pattern(mean.prob.simu[1,], effn.subs=effn.sub.simu1,n.rep=fn.rep)
  simu.sam2= simu.ustat.pattern(mean.prob.simu[1,], effn.subs=effn.sub.simu2,n.rep=fn.rep)
  }   else {simu.sam1=list(rep(NULL,fn.rep));
  for(i in 1:num.fold) simu.sam1=mapply(c, simu.sam1, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=effn.sub.simu1[i,],n.rep=fn.rep), SIMPLIFY=FALSE)
  simu.sam2=list(rep(NULL,fn.rep));
  for(i in 1:num.fold) simu.sam2=mapply(c, simu.sam2, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=effn.sub.simu1[i,],n.rep=fn.rep), SIMPLIFY=FALSE)
  }

  simu.tab.list<-mapply(rbind,simu.sam1,simu.sam2,SIMPLIFY=FALSE)

  aa.temp2<-unlist(lapply(simu.tab.list, chi.stat));

  f.cri<-sort(aa.temp2)[(1-alpha)*fn.rep];
  f.stat<-chi.stat(est.prob[,id.keep]);
  pval<-sum(aa.temp2>f.stat)/fn.rep;
  result<-c(f.cri,f.stat,pval); names(result)<-c('critical.value(5%)', 'chi-stat','pvalue')
  return(result)

}


#' pow.ana.gen.decision function
#'
#' This function determines how powerful the test is (description)
#' @param mean.prob1,mean.prob2,effn.subsam1,effn.subsam2,N.rep,boot.rep (arguments here)
#' @details (details here; DIU)
#' @return c(mean(out[,1]),sd(out[,1]),sum(out[,3]<.05)/N.rep) (return value here)
#' @note (notes here; DIU)
#' @references (references here; DIU)
#' @importFrom stats "sd"
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' fn.rep<-10^3
#' num.fold<-1
#' pow.ana.gen.decision(mean.prob1=c(.8,.2,.5),mean.prob2=c(.4,.3,.2),N.rep=10^2, boot.rep=10^3)

pow.ana.gen.decision<-function(mean.prob1=c(.8,.2,.5),mean.prob2=c(.8,.2,.5),effn.subsam1=c(5,5,5,5),effn.subsam2=c(5,5,5,5),N.rep=10^3, boot.rep=10^3)
{
  n.col=length(mean.prob1);

  simu.sam1=list(rep(NULL,N.rep));mean.prob.simu=cbind(mean.prob1,1-mean.prob1)
  for(i in 1:n.col) simu.sam1=mapply(c, simu.sam1, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=c(effn.subsam1[i],effn.subsam1[i+1]),n.rep=fn.rep), SIMPLIFY=FALSE)
  simu.sam2=list(rep(NULL,N.rep));mean.prob.simu=cbind(mean.prob2,1-mean.prob2)
  for(i in 1:num.fold) simu.sam2=mapply(c, simu.sam2, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=c(effn.subsam2[i],effn.subsam2[i+1]),n.rep=fn.rep), SIMPLIFY=FALSE)

  simu.tab.list<-mapply(rbind,simu.sam1,simu.sam2,SIMPLIFY=FALSE)
  dec<-lapply(c(1:N.rep),  function(x) {#print(x); flush.console();
    gen.decision(matrix(unlist(simu.tab.list[x]),nrow=2), effn.subsam1=effn.subsam1,effn.subsam2=effn.subsam1,fn.rep=boot.rep,alpha=.05);})




  out<-matrix(unlist(dec),nrow=N.rep, byrow=TRUE);
  return(c(mean(out[,1]),sd(out[,1]),sum(out[,3]<.05)/N.rep));
}


#' sub.test function
#'
#' This function compares two samples together (description)
#' @param sam1,sam2 (arguments here)
#' @details (details here; DIU)
#' @return result (return value here)
#' @note (notes here; DIU)
#' @references (references here: DIU)
#' @export
#' @examples
#' fn<-file.choose(); ### C:\Users\wangy\Documents\Research Since 2013\Research 2016\Sub group trend comparison
#' dat<-read.csv(fn, header=TRUE);
#' attach(dat)
#'
#' Lev.TN<-levels(TreatmentName);
#' Lev.Line<-levels(Line);
#' n<-dim(dat)[1]
#'
#' idx<-(TreatmentName==Lev.TN[1])*(Line==Lev.Line[9])*(1:n)
#' boxplot(seedwt[idx]~Env[idx], main=paste(Lev.TN[1],'and',Lev.Line[9]),xlab="ENV levles",ylab='seedwt');

#' idx2<-(TreatmentName==Lev.TN[2])*(Line==Lev.Line[9])*(1:n)
#' boxplot(seedwt[idx2]~Env[idx2], main=paste(Lev.TN[2],'and',Lev.Line[9]),xlab="ENV levles",ylab='seedwt');
#'
#' temp<-seedwt[idx];lab<-Env[idx]; uni.lab<-unique(lab)
#' sam.Ga.Mo352<-lapply(1:length(uni.lab), function(x) temp[lab==uni.lab[x]])
#' sam.Ga.Mo352 ### some are missing
#' avg.sam.Ga.Mo352<-unlist(lapply(sam.Ga.Mo352,mean,na.rm=TRUE));avg.sam.Ga.Mo352
#' sd.sam.Ga.Mo352<-unlist(lapply(sam.Ga.Mo352,sd,na.rm=TRUE));sd.sam.Ga.Mo352
#' size.sam.Ga.Mo352<-unlist(lapply(sam.Ga.Mo352,length));size.sam.Ga.Mo352
#' ####
#' #### find the sd of second group;
#' temp2<-seedwt[idx2];lab2<-Env[idx2]; uni.lab2<-unique(lab2)
#' sam.nh.Mo352<-lapply(1:length(uni.lab2), function(x) temp2[lab2==uni.lab2[x]])
#' sam.nh.Mo352 ### some are missing and some are zero
#' avg.sam.nh.Mo352<-unlist(lapply(sam.nh.Mo352,mean,na.rm=TRUE));avg.sam.nh.Mo352
#' sd.sam.nh.Mo352<-unlist(lapply(sam.nh.Mo352,sd,na.rm=TRUE));sd.sam.nh.Mo352;
#' size.sam.nh.Mo352<-unlist(lapply(sam.nh.Mo352,length));size.sam.nh.Mo352
#' sub.test(sam.nh.Mo352,sam.Ga.Mo352)


sub.test<-function(sam1,sam2)
{
  est.prob.1<-multi.prob(sam1); est.prob.2<-multi.prob(sam2);
  est.prob<-rbind(est.prob.1,est.prob.2)

  effn.subsam1=unlist(lapply(c(1:length(sam1)),function(x) length(setdiff(sam1[[x]],NA))))
  effn.subsam2=unlist(lapply(c(1:length(sam2)),function(x) length(setdiff(sam2[[x]],NA))))

  result<- gen.decision(est.prob,effn.subsam1=effn.subsam1,effn.subsam2=effn.subsam2,fn.rep=10^4,alpha=.05)
  return(result)

}


#' experimental data
#' This assigns the data to variables so that they can be called and read

fn<-file.choose(); ### C:\Users\wangy\Documents\Research Since 2013\Research 2016\Sub group trend comparison
seeddata<-read.csv(fn, header=TRUE);
attach(seeddata)

Lev.TN<-levels(TreatmentName);
Lev.Line<-levels(Line);
n<-dim(seeddata)[1]


#' visual presentation of the data
#'

############ work on Line of Mo352 with GA vs nohormone ##########

#' Visual presentation of data
#' Shows how we generated boxplots
#' Example code after data is introduced
par(mfrow=c(1,2))
idx<-(TreatmentName==Lev.TN[1])*(Line==Lev.Line[9])*(1:n)
boxplot(seedwt[idx]~Env[idx], main=paste(Lev.TN[1],'and',Lev.Line[9]),xlab="ENV levles",ylab='seedwt');

idx2<-(TreatmentName==Lev.TN[2])*(Line==Lev.Line[9])*(1:n)
boxplot(seedwt[idx2]~Env[idx2], main=paste(Lev.TN[2],'and',Lev.Line[9]),xlab="ENV levles",ylab='seedwt');



