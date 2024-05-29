.check.args<-function(param.args,chain.args,inpar.args){
  tmp<-c(param.args,chain.args,inpar.args)
  dups<-duplicated(tmp)
  if(any(dups)){
    dups<-tmp[dups]
    stop('Multiple instances of argument ',paste(dups,collapse=', '),' found')
  }
}

.proc.x<-function(x,fit){
  if(!is.list(x)){
    x<-list(x)
  }
  .expand.par(do.call(par.c,c(x,fit=list(fit))))
}

.Ops.ls<-function(){
  c("+", "-", "*", "/", "^", "%%", "%/%",
    "&", "|", "!","==", "!=", "<", "<=", ">=", ">")
}

.get.param.names<-function(nms,overwrites){
  param.names<-gsub(paste0('%(\\',
                           paste(.Ops.ls(),collapse='|\\'),
                           ')%'),
                    ' \\1 ',
                    nms)
  if(!is.null(overwrites)){
    if(length(overwrites)==length(param.names)){
      param.names<-overwrites
    }else{
      warning('Mistmatch between number of parameters to be plotted and provided parameter name overwrites: parameter name overwrites were ignored')
    }
  }
  param.names
}

#bg only sets plot background in calls to par(); otherwise, it is used to set coloring for points!
.recyclables<-function(){
  c('bg','cex','col','lty','lwd','pch',
    'border','density','angle',
    'alpha','border.alpha',
    'lines.col','lines.lty','lines.lwd',
    'lines.alpha',
    'ribbon.col','ribbon.lty','ribbon.lwd',
    'ribbon.border','ribbon.density','ribbon.angle',
    'ribbon.alpha','ribbon.border.alpha')
}

#to prevent errors, may be good to hand-curate a list of allowable graphical parameters to be recycled...
#wow, surprisingly short list!
.make.args.fxn<-function(...,
                         param.args,chain.args,inpar.args,
                         n.params,n.chains,n.inpar,
                         inheritance,
                         defaults,
                         trace=FALSE){
  args.ls<-list(...)
  if(trace){
    inds<-names(args.ls)%in%c("border","density","angle","border.alpha")
    names(args.ls)[inds]<-paste0("ribbon.",names(args.ls)[inds])
  }
  args.ls<-args.ls[!(names(args.ls)%in%c('x','y',
                                         'x0','y0','x1','y1',
                                         'plot','warn.unused',
                                         'add'))]
  #defaults--breaks if back-reference is used earlier in sequence than it's value is set!
  #(for example, problems would result if border='.col'   is specified before col's default)
  #changed to more elegant expression based solution
  #sorting would be the way to solve this issue, but would be overly-complicated I think
  #probably just something to warn about when using function
  #the gist here is to set defaults to alternative values for any base defaults that work poorly for the parameter plot
  for(i in names(defaults)){
    if(is.null(args.ls[[i]])){
      if(is.expression(defaults[[i]])){
        args.ls[[i]]<-tryCatch(with(args.ls,eval(defaults[[i]])),error=function(e) NULL)
      }else{
        args.ls[[i]]<-defaults[[i]]
      }
    }
  }
  #compiling names list and inheriting arguments
  args.nms<-list(plot=c(names(formals(plot.default)),
                        names(formals(axis)),
                        names(formals(box)),
                        names(formals(plot.window)),
                        names(formals(title))),
                 gen=names(par()),
                 dens=names(formals(density.default)),
                 hist=names(formals(hist.default)),
                 poly=names(formals(polygon)),
                 legend=names(formals(legend2)))
  #add some custom parameters
  args.nms[['gen']]<-c(args.nms[['gen']],'alpha')
  args.nms[['poly']]<-c(args.nms[['poly']],'border.alpha')
  args.nms<-lapply(args.nms,function(ii) ii[ii!='...'])
  for(i in names(inheritance)){
    if(!is.null(inheritance[[i]])){
      #if you want it to inherit things...
      #first figure out which args have been PRESPECified using appropriate prefixes...
      #make sure to ONLY inherit arguments for which nothing has been prespecified
      incl.args<-unlist(args.nms[inheritance[[i]]])
      args.nms[[i]]<-paste0(i,'.',incl.args)
      prespec<-names(args.ls)%in%args.nms[[i]]
      excl.args<-gsub(paste0('^',i,'\\.'),'',names(args.ls)[prespec])
      incl.args<-incl.args[!(incl.args%in%excl.args)]
      inds<-names(args.ls)%in%incl.args
      #not make these arguments get inherited by adding the appropriate prefix
      if(any(inds)){
        add.ls<-args.ls[inds]
        names(add.ls)<-paste0(i,'.',names(add.ls))
        args.ls<-c(args.ls,add.ls)
      }
    }else{
      #otherwise, don't do anything! but make sure the only arguments are those with the appropriate prefix...
      args.nms[[i]]<-paste0(i,'.',args.nms[[i]])
    }
  }
  #recycling
  recycle.args<-list(param.args=param.args,
                     chain.args=chain.args,
                     inpar.args=inpar.args)
  recyclables<-.recyclables()
  recycle<-function(x,ind){
    if(!is.null(x)){
      out<-switch(ind,
                  rep(rep(x,length.out=n.params),n.chains),
                  rep(rep(x,length.out=n.chains),each=n.params),
                  rep(rep(x,length.out=n.inpar),n.params*n.chains))
      if(ind<3&n.inpar){
        out<-rep(out,each=n.inpar)
      }
      out
    }
  }
  for(i in seq_along(recycle.args)){
    inds<-recycle.args[[i]]%in%recyclables
    if(length(inds)){
      recycle.args[[i]]<-recycle.args[[i]][inds]
      args.ls[recycle.args[[i]]]<-lapply(args.ls[recycle.args[[i]]],recycle,ind=i)
    }
  }
  #take care of alpha stuff
  vec1<-c('',paste0(names(inheritance),'.'))
  vec2<-paste0(c('alpha','border.alpha','bg.alpha'))
  vec3<-paste0(c('col','border','bg'))
  nn<-seq_along(vec2)
  for(i in nn){
    for(j in nn){
      alpha.nm<-paste0(vec1[i],vec2[j])
      tmp.alpha<-args.ls[[alpha.nm]]
      if(!is.null(tmp.alpha)){
        col.nm<-paste0(vec1[i],vec3[j])
        tmp.col<-args.ls[[col.nm]]
        args.ls[[col.nm]]<-alter.cols(tmp.col,tmp.alpha)
        args.ls[[alpha.nm]]<-NULL
      }
    }
  }
  args.ls<-args.ls[lengths(args.ls)>0] #drop NULL arguments--might cause problems down the line
  args.nms<-c(args.nms,recyc=list(unlist(recycle.args)))
  args.inds<-args.nms
  for(i in seq_along(args.inds)){
    args.inds[[i]]<-names(args.ls)%in%args.nms[[i]]
  }
  args.inds[['noncol']]<-!(names(args.ls)%in%c('col','ribbon.col'))
  args.inds[['nonborder']]<-!(names(args.ls)%in%c('border','ribbon.border'))
  for(i in names(inheritance)){
    names(args.ls)<-gsub(paste0('^',i,'\\.'),'',names(args.ls))
  }
  out<-function(form,ind=NULL){
    res<-args.ls[with(args.inds,eval(form[[2]]))]
    if(!is.null(ind)){
      res<-lapply(res,'[',ind)
    }
    res
  }
  out
}

.get.lines<-function(x,add.lines,quant.lines,n.tot){
  lines<-rep(list(add.lines),n.tot)
  if(!is.null(add.lines)){
    means<-is.na(add.lines)
    if(any(means)){
      tmp.lines<-asplit(.call.op('means',
                                 list(chains=x),
                                 list('.',add.lines[!means])),
                        -1)
      for(i in seq_len(n.tot)){
        lines[[i]][means]<-tmp.lines[[i]]
      }
    }
    if(quant.lines){
      tmp.lines<-asplit(.call.op('quantiles',
                                 list(chains=x),
                                 list('.',add.lines[!means])),
                        -1)
      for(i in seq_len(n.tot)){
        lines[[i]][!means]<-tmp.lines[[i]]
      }
    }
  }
  lines
}

#this function could definitely use some cleaning up...
.get.cuts<-function(x,
                    p=0.05,
                    lower.quant=NULL,upper.quant=NULL,
                    lower.cut=NULL,upper.cut=NULL,
                    n.params,n.chains,
                    trace=FALSE,ribbon=TRUE,add.lines=matrix(nrow=1,ncol=0),
                    n.windows=NULL,window.size=NULL){
  n<-dim(x)[1]
  p.seq<-seq_len(n.params)
  c.seq<-seq_len(n.chains)
  if(is.null(p)){
    lower.prob<-upper.prob<-rep(NA,length.out=n.params)
  }else{
    lower.prob<-rep(p/2,length.out=n.params)
    upper.prob<-1-lower.prob
  }
  cuts.list<-list(lower.quant,upper.quant,lower.cut,upper.cut)
  if(!trace){
    cuts.list<-c(cuts.list,list(lower.cut),list(upper.cut))
  }
  for(i in seq_along(cuts.list)){
    if(i==3){
      #not using %quantiles% since probs can differ by parameter!
      if(!trace){
        prob2cut<-unlist(lapply(seq_len(n.params),
                                function(ii) apply(x[,ii,,drop=FALSE],
                                                   -1,
                                                   .int.quant,
                                                   n=n,
                                                   p=c(lower.prob[ii],
                                                       upper.prob[ii]),
                                                   sorted=FALSE)))
        inds<-2*(rep(c.seq,each=n.params)+(p.seq-1)*n.chains)
        lower.prob<-pmin(prob2cut[inds-1],prob2cut[inds])
        upper.prob<-pmax(prob2cut[inds-1],prob2cut[inds])
      }else{
        tmp.prob<-matrix(nrow=n.params,ncol=2)
        tmp.prob[,1]<-pmin(lower.prob,upper.prob)
        tmp.prob[,2]<-pmax(lower.prob,upper.prob)
        tmp.prob<-tmp.prob[rep(p.seq,each=n.chains),,drop=FALSE]
        prob2cut<-.roll.quant(x,nw=n.windows,w=window.size,p=if(ribbon) cbind(tmp.prob,add.lines) else add.lines)
      }
    }
    tmp.cut<-cuts.list[[i]]
    if(!is.null(tmp.cut)){
      tmp.cut<-rep(tmp.cut,length.out=n.params)
      nonNAs<-!is.na(tmp.cut)
      if(!i%%2){
        upper.prob[nonNAs]<-tmp.cut[nonNAs]
      }else{
        lower.prob[nonNAs]<-tmp.cut[nonNAs]
      }
    }
  }
  if(!trace){
    lower.prob[is.na(lower.prob)]<-par()$usr[1]-1
    upper.prob[is.na(upper.prob)]<-par()$usr[2]+1
    cuts<-asplit(matrix(c(pmin(lower.prob,upper.prob),pmax(lower.prob,upper.prob),add.lines),n.params*n.chains,2+ncol(add.lines)),1)
  }else{
    cuts<-list(prob2cut[1],prob2cut[-1])
  }
  cuts
}

.make.param.plot<-function(x,
                           param.names,n.params,
                           def.xlab=NULL,smooth=NULL,
                           ...){
  args.ls<-list(...)
  args.ls[['pch']]<-NULL
  if(!hasArg(xlim)){
    if(is.null(smooth)){
      args.ls[['xlim']]<-c(1,dim(x)[1])
    }else{
      args.ls[['xlim']]<-range(unlist(lapply(x,'[[','x')),na.rm=TRUE)
    }
  }
  if(!hasArg(ylim)){
    if(is.null(smooth)){
      args.ls[['ylim']]<-range(as.vector(x),na.rm=TRUE)
    }else{
      args.ls[['ylim']]<-c(0,max(unlist(lapply(x,'[[','y')),na.rm=TRUE))
    }
  }
  if(n.params>5|any(nchar(param.names)>50)){
    tmp.lab<-def.xlab
  }else{
    tmp.nms<-param.names
    modes<-unlist(lapply(tmp.nms,mode))
    not.expr<-modes=='numeric'|modes=='character'
    tmp.nms[not.expr]<-lapply(tmp.nms[not.expr],function(ii) paste0('"',ii,'"'))
    tmp.lab<-str2expression(paste(as.character(tmp.nms),collapse="*'; '*"))
  }
  if(!hasArg(xlab)){
    if(is.null(smooth)){
      args.ls[['xlab']]<-'Iteration'
    }else{
      args.ls[['xlab']]<-tmp.lab
    }
  }
  if(!hasArg(ylab)){
    if(is.null(smooth)){
      args.ls[['ylab']]<-tmp.lab
    }else{
      args.ls[['ylab']]<-if(smooth) 'Density' else 'Frequency'
    }
  }
  do.call(plot,
          c(x=list(0),
            pch=NA,
            args.ls))
}

.get.y<-function(ind,x,y,xout,smooth){
  ind[ind==1|!ind]<-NA
  ind2<-ind-1
  if(smooth){
    out<-(y[ind]-y[ind2])/(x[ind]-x[ind2])*(xout-x[ind2])+y[ind2]
  }else{
    out<-y[ind2]
  }
  out[is.na(out)]<-0
  out
}

.fix.pts<-function(x,y,smooth){
  n.pts<-length(x)
  if(smooth){
    xs<-c(x[1],x,x[n.pts])
    ys<-c(0,y,0)
  }else{
    xs<-rep(x,each=2)
    ys<-c(0,rep(y,each=2),0)
  }
  list(xs,ys)
}

.trim.prof<-function(x,y,smooth,xran,yran){
  n.pts<-length(y)
  if(smooth){
    yran<-
      zz<-y<yran/1000&
      (x<xran[1]|x>xran[2])
    zz<-rle(zz) #based on plotting window--should change with zoom level and such!
    #a more complicated scheme might also condition on whether x values are observed past a certain range
    #but I think this works for most purposes!
    #I think I did this correctly...
  }else{
    zz<-rle(!y)
  }
  if(zz$values[1]){
    l.pts<-seq_len(zz$lengths[1])
    if(smooth) l.pts<-l.pts[-zz$lengths[1]]
  }else{
    l.pts<-integer(0)
  }
  len<-length(zz$values)
  if(zz$values[len]){
    r.pts<-seq.int(n.pts-zz$lengths[len]+1,n.pts)
    if(smooth) r.pts<-r.pts[-1]
  }else{
    r.pts<-integer(0)
  }
  inds<-c(l.pts,r.pts)
  if(length(inds)){
    y<-y[-inds]
    if(smooth){
      x<-x[-inds]
    }else{
      if(length(l.pts)){
        inds[-l.pts]<-inds[-l.pts]+1
      }else{
        inds<-inds+1
      }
      x<-x[-inds]
    }
  }
  list(list(x,y),.fix.pts(x,y,smooth))
}

.cut.prof<-function(xy,cuts,smooth){
  x<-xy[[1]]
  y<-xy[[2]]
  n.pts<-length(x)
  cut.pts<-findInterval(cuts,x)
  part.x<-x
  part.y<-y
  l<-cuts[1]
  r<-cuts[2]
  l.pt<-cut.pts[1]
  r.pt<-cut.pts[2]
  if(cut.pts[1]>0){
    inds<-seq.int(l.pt+1,n.pts)
    r.pt<-r.pt-l.pt+1
    part.y<-c(.get.y(inds[1],part.x,part.y,l,smooth),part.y[inds])
    if(!smooth) part.y<-part.y[-length(part.y)]
    part.x<-c(l,part.x[inds])
  }
  if(cut.pts[2]<n.pts){
    inds<-seq.int(1,r.pt)
    part.y<-c(part.y[inds],.get.y(length(inds)+1,part.x,part.y,r,smooth))
    if(!smooth) part.y<-part.y[-length(part.y)]
    part.x<-c(part.x[inds],r)
  }
  .fix.pts(part.x,part.y,smooth)
}

.make.param.plot.legend<-function(param.names,chain.names,
                                  n.params,n.chains,n.inpar,
                                  get.args,
                                  include.chain.legend,
                                  trace=FALSE,ribbon=TRUE,ribbon.first=NULL){
  legend.args<-get.args(~legend)
  legend.args$ribbon.first <- ribbon.first
  if(is.null(legend.args$x)){
    legend.args$x<-'topright'
  }
  if(is.null(legend.args$legend)){
    add.legend<-TRUE
  }
  if(n.params>1){
    inds<-if(n.inpar) n.inpar*(seq_len(n.params)-1)+1 else seq_len(n.params)
    if(add.legend){
      legend.args$legend<-param.names
    }
  }else{
    inds<-integer(0)
  }
  if(n.chains>1&include.chain.legend){
    inds<-c(inds,
            if(n.inpar) seq.int(1,n.params*n.inpar*n.chains,n.params*n.inpar) else seq.int(1,n.params*n.chains,n.params))
    if(add.legend){
      legend.args$legend<-c(legend.args$legend,chain.names)
    }
  }
  if(!is.null(legend.args$legend)){
    #works for prof.plot--will have to revise for trace.plot
    if(trace){
      if(ribbon){
        tmp.args<-c(get.args(~ribbon&!recyc),get.args(~ribbon&recyc,inds))
      }else{
        tmp.args<-list()
      }
    }else{
      tmp.args<-c(get.args(~(gen|poly)&!recyc),get.args(~(gen|poly)&recyc,inds))
    }
    #keep only arguments that CAN be recycled--otherwise you get unexpected inheritance!
    #this should cover all symbol plotting stuff for now, but it may have to be modified in the future!
    tmp.args<-tmp.args[names(tmp.args)%in%.recyclables()]
    names(tmp.args)[names(tmp.args)=='col']<-'fill'
    names(tmp.args)[names(tmp.args)=='lty']<-'fill.lty'
    names(tmp.args)[names(tmp.args)=='lwd']<-'fill.lwd'
    if(trace){
      tmp.args<-c(tmp.args,get.args(~gen&!recyc),get.args(~gen&recyc,inds))
      if(is.null(legend.args$merge)) legend.args$merge<-TRUE
      if(is.null(legend.args$seg.len)) legend.args$seg.len<-0.8
    }
    #overwrites
    tmp.args<-tmp.args[!(names(tmp.args)%in%names(legend.args))]
    legend.args<-c(legend.args,tmp.args)
    #found it necessary to add this to prevent certain extraneous parameters from getting passed and breaking the function
    do.call(legend2,legend.args[names(legend.args)%in%names(formals(legend2))])
  }
}

#should now work with unsorted vectors, NA-containing vectors, NA ps, ps greater than 1 or less than 0
#should be able to be integrated with quantile functions
.int.quant<-function(x,n,p,sorted=TRUE){
  if(!sorted){
    x<-sort(x,method='quick',na.last=TRUE)
  }
  if(is.na(x[n])){
    nas<-is.na(x)
    x<-x[!nas]
    n<-n-sum(nas)
  }
  if(!n){
    rep(NA,length(p))
  }else{
    x<-c(x,0)
    p[p>1]<-1
    p[p<0]<-0
    tmp<-p*(n-1)+1
    j<-floor(tmp)
    g<-tmp-j
    (1-g)*x[j]+g*x[j+1]
  }
}

.get.windows<-function(w,n,nn){
  halfw<-round(w/2)
  inds<-seq.int(-halfw,halfw)
  foo<-function(i){
    out<-i+inds
    out[out>=1&out<=n]
  }
  lapply(nn,foo)
}

.roll.quant<-function(x,p,nw,w){
  n<-dim(x)[1]
  nn<-round(seq(1,n,length.out=nw))
  winds<-.get.windows(w,n,nn)
  ws<-lengths(winds)
  x<-asplit(x,-1)
  ranks<-lapply(x,rank)
  np<-ncol(p)
  holder<-rep(NA,np)
  p.seq<-seq_len(np)
  foo<-function(i){
    x<-x[[i]]
    ranks<-ranks[[i]]
    nas<-is.na(p[i,])
    means<-nas&p.seq>2
    int.foo<-function(j){
      inds<-winds[[j]]
      x<-x[inds]
      ranks<-ranks[inds]
      sorts<-sort.int(ranks,method='quick',index.return=TRUE)$ix
      holder[!nas]<-.int.quant(x[sorts],ws[j],p[i,!nas])
      holder[means]<-mean(x)
      holder
    }
    do.call(rbind,lapply(seq_len(nw),int.foo))
  }
  out<-c(list(nn),lapply(seq_along(x),foo))
  out
}

#slower
# .roll.quant<-function(x,w,p){
#   n<-dim(x)[1]
#   nn<-seq_len(n)
#   winds<-.get.windows(w)
#   ws<-lengths(winds)
#   x<-asplit(x,-1)
#   foo<-function(i){
#     x<-x[[i]]
#     int.foo<-function(j){
#       inds<-winds[[j]]
#       x<-x[inds]
#       .int.quant(sort(x,method='radix'),ws[j],p)
#     }
#     do.call(rbind,lapply(nn,int.foo))
#   }
#   out<-lapply(seq_along(x),foo)
# }

#slowest
# .roll.quant<-function(x,w,p){
#   n<-dim(x)[1]
#   nn<-seq_len(n)
#   winds<-.get.windows(w)
#   ws<-lengths(winds)
#   x<-asplit(x,-1)
#   foo<-function(i){
#     x<-x[[i]]
#     ranks<-ranks[[i]]
#     int.foo<-function(j){
#       inds<-winds[[j]]
#       x<-x[inds]
#       ranks<-ranks[inds]
#       x<-data.table(x=x,r=ranks,key='r')$x
#       .int.quant(x,ws[j],p)
#     }
#     do.call(rbind,lapply(nn,int.foo))
#   }
#   out<-lapply(seq_along(x),foo)
# }

# microbenchmark::microbenchmark(.roll.quant(fit%c%1:5,100,100,c(0.025,0.975)),times=10) #takes about a second--not too bad
# #a nice thing is that increasing window sizes seem to have relatively minor effects on its speed--going from 100 to 500 only doubled the time
# #might be a nice option to offer skips?
# test<-.roll.quant(fit%c%1:5,c(0.025,0.975))
# plot((fit%c%1)[,1],type='l')
# matplot(test[[1]],x=seq.int(1,1500,length.out=100),type='l',add=TRUE)
# #yeah, offering skips is simple, doesn't sacrifice visual quality much, and offers tremendous speed benefits