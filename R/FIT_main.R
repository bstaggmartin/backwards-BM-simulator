#new plan (only thinking univariate-ly for now):
##trait.data and trait.se are separate objects
##Still must be labelled
##Can be named vectors, rownamed column matrices, or data.frames with a "tip.label" column
##NAs in trait.data are taken to be missing data (coded internally as 0)
##NAs in SE are taken to indicate unfixed SEs to be estimated during model fitting (coded internally as -1)
##Infs in SE are taken to indicate missing data (coded internally as -2)
##After coercing trait.data and trait.se to rownamed matrices...
####First make sure tree and trait.data is alright--lump together tips separated by 0 distance, remove trait.data with no
######corresponding tip
####Priority of conflicts for SE handling:
########tips with intra-tip observations (>1 observations) and non-unfixed SEs return warnings and set SEs to unfixed
########tips with any observations (>0 observations) and infinite SEs return warnings and set SEs to unfixed
########tips with no observations (0 observations) and a non-infinite SE return warnings and set SEs to infinite
####After making these corrections, check for any double SE specification in which there are two valid SEs specified for
######the same tip, and return error if found
##Additionally, might be good to make T/F "switch" that instead estimates fixed SE for each tip using intra-tip observations and
###using these as fixed SEs, rather than estimating them during model fitting
####Will need to have an additional argument specifying how to handle tips with single observations in this case
########Probably just something you can set to NA (for unfixing these) or some number (for default SE)
########Maybe NULL will use max SD found when calculating SEs for other tips?
##Handle priors with ...?
##Technically, standard errors in this context are standard variance --> should probably square given SEs to make this less
###confusing.

#for some reason, can't use '~/../' in output directory???

#' Prepare input data for an Evolving Rates model
#'
#'
#' This function processes tree and trait data in preparation for fitting an evorates model. Also takes additional model inputs
#' (e.g., priors), though these may be overwritten later by \code{run.evorates()}.
#'
#'
#' @param tree An object of class "\code{phylo}"
#' @param trait.data Three options:
#' \itemize{
#' \item{A named vector of trait values.}
#' \item{A rownamed matrix of trait values. For the moment, must be 1 column since multivariate models are not
#'       yet supported.}
#' \item{A data.frame with 2 columns: 1 with numeric data which is interpreted as trait values and 1 with
#'       string/factor data interpreted as names. If more than 1 of either kind of column is found, returns error.
#'       Will also use rownames if no string/factor column is found, but this limits data to consist of only 0-1
#'       observations per tip.}
#' }
#' Both multiple observations for a single tip and missing observations are allowed.
#' In all cases, the associated names must generally match the tip labels found in \code{tree} (\code{tree$tip.label})
#' exactly. The only exception are numbers corresponding to internal node indices (these can be viewed by running:
#' \code{plot(tree); nodelabels()}), which are used to assign trait values to internal nodes. As of now, the function
#' actually completely ignores \code{tree$node.label}--I should probably fix that and plan to in the future. Be aware
#' that this feature slightly alters \code{tree} because it works by grafting extra tips of length 0 to the appropriate
#' nodes.
#' @param trait.se A vector, matrix, or data.frame of trait value standard errors which must be unambiguously 
#' labeled (see \code{trait.data}). Alternatively, a single, unlabeled number that will be applied to all tips in
#' \code{tree} (e.g., you could set it to 0 to specify all trait values are known without error). If \code{NULL}
#' (the default), standard error is estimated for all tips in the tree while fitting the model. There are several
#' things to note here: 
#' \itemize{
#' \item{Unlike \code{trait.data}, there can only be 1 trait standard error per tip.}
#' \item{You can use \code{NA} to specify tips that you want to estimate standard error for and \code{Inf} to
#'       specify tips that have missing trait values.}
#' \item{Generally, tips with multiple observations and no observations are automatically assigned standard errors
#'       of \code{NA} and \code{Inf}, respectively.}
#' \item{Any tips with unspecified standard error default to \code{NA}.}
#' }
#' Any conflicting/impossible standard error specifications are corrected and return warnings.
#' @param constrain.Rsig2 \code{TRUE} or \code{FALSE}: should the rate varinace (\code{R_sig2}) parameter be constrained to 0,
#' resulting in a simple Brownian Motion or early/late burst model? Defaults to \code{FALSE}. See Details for a
#' definition of model parameters.
#' @param trend \code{TRUE} or \code{FALSE}: should a trend in rates over time (\code{R_mu}) be estimated,
#' resulting in an early/late burst or trended evorates model? Defaults to \code{FALSE}.
#' @param lik.power A single number between 0 and 1 specifying what power to raise the likelihood function to.
#' This is useful for sampling from "power posteriors" that shift the posterior to look more like the prior (indeed,
#' you can set this to 0 to sample from the prior itself). Useful for model diagnostics and calculating things
#' like Bayes Factors. Technically, you can set \code{lik.power} above 1 (a technique called "data cloning") but
#' this is not beneficial for evorates models and we don't recommend it.
#' @param sampling.scale \code{TRUE} or \code{FALSE}: should provided prior parameters (see \code{...}) be
#' interpreted on the raw scale of the data or on the transformed scale? All data passed to the Stan-based HMC
#' sampler are transformed such that \code{tree}'s total height is 1 and the variance of the trait data
#' is 1. Defaults to \code{FALSE}, such that prior parameters are interpreted on the untransformed scale.
#' @param ... Prior arguments (see details for further information on what parameters mean):
#' \itemize{
#' \item{In general, prior standard deviations and, in some cases, means, can be tweaked by passing arguments
#'       named "<parameter_name>_sd" and "<parameter_name>_mean", respectively. For example:
#'       \code{R_mu_mean = 1} or \code{R_sig2_sd = 1}.}
#' \item{Prior on rate variance (\code{R_sig2}): Follows a half-Cauchy distribution (basically a half-normal distribution 
#'       with extremely fat tails) with adjustable standard deviation. By default, standard deviation is set to 5 divided
#'       by the maximum height of \code{tree} (i.e., 5 on the sampling scale).}
#' \item{Prior on trend (\code{R_mu}): Follows a normal distribution with adjustable mean and standard deviation. By
#'       default, mean is set to 0 (also 0 on sampling scale), and standard deviation is set to 10 divided by the
#'       maximum height of \code{tree} (i.e., 10 on the sampling scale).}
#' \item{Prior on rate at the root (\code{R_0}): Follows a normal distribution with adjustable mean and standard deviation.
#'       By default, mean is set to variance of \code{trait.data} divided by the maximum height of \code{tree} on the
#'       natural log scale (i.e., 0 on sampling scale), and standard deviation is set to 10 (also 10 on sampling scale).}
#' \item{Prior on tip error variance (\code{Y_sig2}): Follows a truncated t distribution (basically a normal
#'       distribution with potentially fatter tails) with adjustable mean, standard deviation, and degrees of freedom.
#'       The distribution is truncated at 0 to the left such that Y_sig2 estimates are always positive (note that setting
#'       the mean to 0 results in a half-t distribution). By default, mean is set to 0, standard deviation is set to 2
#'       times the variance of \code{trait.data} (i.e., 2 on the sampling scale), and degrees of freedom is set to 1. The
#'       degrees of freedom can be tweaked with the argument "Ysig2_df" and may be any positive number, including \code{Inf},
#'       with lower numbers leading to a less informative priors with fatter tails. Notably, setting "Ysig2_df" to \code{Inf}
#'       makes this prior a truncated normal and much more informative!}
#' }
#' 
#' 
#' @return A list of \code{call}, mainly containing information on the inputted tree and trait data, \code{trans.const},
#' containing transformation constants used to scale data (see \code{sampling.scale}), and \code{dat}, containing the
#' processed data that will be passed to directly to the Stan-based Hamiltonian Monte Carlo sampler.
#' 
#' 
#' #' @details Parameter definitions:
#' \itemize{
#' \item{Rate variance (\code{R_sig2}, also denoted \eqn{\sigma^2_{\sigma^2}}{\sigma^2[\sigma^2]}): Determines how much random variation
#'       accumulates in rates over time. Specifically, independent lineages evolving for a length of time \code{t} would exhibit log-normally
#'       distributed rates with standard deviation \code{sqrt(t * R_sig2)}. Another way to think about this is that the 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time is given by \code{exp(c(-1,1) * 1.96 * sqrt(R_sig2))}, assuming \code{R_mu = 0.}}
#' \item{Trend (\code{R_mu}, also denoted \eqn{\mu_{\sigma^2}}{\mu[\sigma^2]}): Determines whether median rates tend to decrease (if
#'       negative) or increase (if positive) over time. Specifically, the median fold-change in rate for a given lineage is 
#'       \code{exp(t * R_mu)} over a length of time \code{t}. This is distinct from changes in average rates, which depends on both
#'       \code{R_mu} and \code{R_sig2}.}
#' \item{Note on combining rate variance and trend parameters: The 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time, given a trend, is simply
#'       \code{exp(R_mu) * exp(c(-1,1) * 1.96 * sqrt(R_sig2))} or equivalently \code{exp(R_mu + c(-1,1) * 1.96 * sqrt(R_sig2))}.
#'       This distribution of rate change is right-skewed due to the \code{exp()} function. Because of this, median changes in rates
#'       over time will always be lower than average changes. The "average change" parameter, denoted \code{R_del} or
#'       \eqn{\delta_{\sigma^2}}{\delta[\sigma^2]}, is given by \code{R_mu + R_sig2 / 2}. Accordingly, the average fold-change in rate
#'       for a given lineage is \code{exp(t * R_del)} over a length of time \code{t}.}
#' \item{Tip error variance (\code{Y_sig2}, also denoted \eqn{\sigma^2_y}{\sigma^2[y]}): Determines the among of "error" around observed
#'       trait values for tips without fixed standard errors. Specifically, raw observations for a given tip are sampled from a normal
#'       distribution centered at that tip's "true trait value" with standard deviation \code{sqrt(Y_sig2)}.}
#' \item{Rate at the root (\code{R_0}, also denoted \eqn{\sigma^2_0}{\sigma^2[0]}): The natural log of the rate at the root of the entire phylogeny.}
#' \item{Branchwise rates (\code{R_i}, also denoted \eqn{ln \bar{\sigma^2_i}}{ln \sigma^2[i]}, where i is the index of an edge in \code{tree}): 
#'       The natural log of the average rate along branches of the phylogeny. Note that rates are always shifting over time under this model,
#'       and these quantities are thus averages. The true rate value at any particular time point along a branch is a related, but separate,
#'       quantity. These will be NA for branches of length 0.}
#' \item{Background rate (\code{bg_rate}): The natural log of the average trait evolution rate for the entire phylogeny, given by the average of
#'       \code{exp(R_i)}, with each entry weighted by its respective branch length. NAs are ignored in this calculation since they correspond
#'       to branches of length 0.}
#' \item{Branchwise rate deviations (\code{Rdev_i}, also denoted ln \eqn{\bar{\sigma^2_{dev,i}}}{ln \sigma^2[dev,i]}): Determines whether a
#'       branches exhibit relatively "fast" (if positive) or "slow" rates (if negative). Generally, these are the differences
#'       between the branchwise rates and \emph{geometric} background rate on the natural log scale,
#'       though this depends on \code{remove.trend)}. Here, the geometric background rate is defined as the weighted average of \code{R_i},
#'       with weights corresponding to branch lengths. This helps prevent some issues with comparing highly right-skewed distributions.
#'       If \code{remove.trend = TRUE}, then branchwise rates are first "detrended" prior to calculating background rates and deviations
#'       (\code{R_i - (-log(abs(R_mu)) - log(l_i) + log(abs(exp(R_mu * t1_i) - exp(R_mu * t0_i))))}, where \code{l} is a vector of
#'       branch lengths and \code{t0}/\code{t1} are vectors of start and end times of each branch). This basically make branchwise rate deviations
#'       determine whether branches exhibit slow/fast rates \emph{given} the overall trend in rates through time. Otherwise, branchwise rate
#'       deviations will simmply indicate slow/fast branches occur at the root/tips of a tree in the presence of a strong trend.}
#' }
#' 
#' 
#' @seealso \link{fit.evorates} is a convenient wrapper for \link{input.evorates}, \link{run.evorates}, and \link{output.evorates}
#' 
#' 
#' @examples
#' #get whale/dolphin tree/trait data
#' data("cet_fit")
#' tree <- cet_fit$call$tree
#' X <- cet_fit$call$trait.data
#' 
#' #prepare data for evorates model
#' input <- input.evorates(tree = tree, trait.data = trait.data)
#' #run the evorates model sampler (takes a couple minutes)
#' run <- run.evorates(input.evorates.obj = input, out.file = "test_fit", chains = 1)
#' #process the output from the object
#' fit <- output.evorates(run.evorates.obj = run)
#' #process the output from the files
#' fit <- output.evorates(run.evorates.obj = "test_fit")
#' 
#' #see fit.evorates() documentation for further examples
#' 
#' 
#' @export
input.evorates<-function(tree,trait.data,trait.se=NULL,constrain.Rsig2=FALSE,trend=FALSE,lik.power=1,sampling.scale=FALSE,
                         ...){
  ##INITIAL DATA/PRIOR CLEAN-UP##
  #make sure any prior scale parameters are positive
  prior.list<-do.call(.make.prior.list,list(...))
  #coerce trait.data and trait.se to rownamed matrices
  trait.data<-.coerce.trait.mat(trait.data)
  if(!is.null(trait.se)){
    if(length(trait.se==1)&is.null(names(trait.se))&is.numeric(trait.se)){
      trait.se<-matrix(trait.se,length(tree$tip.label),1)
      rownames(trait.se)<-tree$tip.label
    }else{
      trait.se<-.coerce.trait.mat(trait.se)
    }
  }
  #coerce tree to be compatible (collapse internal 0-length edges to polytomies, identical tips to single tips)
  tmp<-.coerce.tree(tree,trait.data,trait.se)
  tree<-tmp$tree
  trait.data<-tmp$trait.data
  trait.se<-tmp$trait.se
  #check for "nits" (labels in trait info found to be NNNot IIIn TTTree)
  tmp<-sapply(rownames(trait.data),function(ii) ii%in%tree$tip.label)
  conflicts_nits<-rownames(trait.data)[!tmp]
  if(length(conflicts_nits)>0){
    trait.data<-trait.data[tmp,,drop=FALSE]
    warning('No tip in tree matched with label(s) ',
            paste(conflicts_nits,collapse=', '),
            ' in trait.data: value(s) with aforementioned label(s) removed.',
            immediate.=TRUE)
  }
  tmp<-sapply(rownames(trait.se),function(ii) ii%in%tree$tip.label)
  conflicts_nits<-rownames(trait.se)[!tmp]
  if(length(conflicts_nits)>0){
    trait.se<-trait.se[tmp,,drop=FALSE]
    warning('No tip in tree matched with label(s) ',
            paste(conflicts_nits,collapse=', '),
            ' in trait.se: standard errors with aforementioned label(s) removed.',
            immediate.=TRUE)
  }
  
  
  ##TRAIT INFO##
  n<-length(tree$tip.label)
  e<-nrow(tree$edge)
  parsed.trait.data<-split(trait.data,rownames(trait.data))
  #return NA is every trait val is NA, otherwise ignore NA entries
  foo<-function(x){
    if(all(is.na(x))){
      NA
    }else{
      mean(x,na.rm=TRUE)
    }
  }
  X<-as.matrix(sapply(parsed.trait.data,foo)[tree$tip.label])
  n_obs<-lengths(parsed.trait.data)[tree$tip.label]
  X[is.na(X)]<-n_obs[is.na(X)]<-0
  perc_mis<-sum(n_obs==0)/n
  if(perc_mis==1){
    stop('No trait information provided: please double-check input and labelling.')
  }
  #will have to change this a bit for the multivariate case
  if(perc_mis>0.25){
    warning('You are about to run a model with greater than ',
            floor(perc_mis*100),
            '% of the tips missing trait information: this can be done, but it could indicate some problems with input/labelling.',
            immediate.=TRUE)
  }
  has_intra<-n_obs>1
  mis<-n_obs==0
  p_SE<-rep(-1,n)
  p_SE[mis]<- -2
  names(X)<-names(n_obs)<-names(p_SE)<-tree$tip.label
  if(length(trait.se)>0){
    #make sure SE info doesn't have any conflicts with trait info
    p_SE<-.coerce.trait.se(trait.se,p_SE,has_intra,mis)
  }
  p_SE[p_SE>0]<-p_SE[p_SE>0]^2
  inv_n_obs<-1/n_obs
  n_contra<-sum(lengths(parsed.trait.data)-1)
  foo<-function(x){
    l<-length(x)
    if(l==1){
      NULL
    }else{
      tmp.seq<-seq_along(x)
      rbind(x[-1]-cumsum(x[-l])/tmp.seq[-l],tmp.seq[-1]/tmp.seq[-l])
    }
  }
  tmp<-do.call(cbind,lapply(parsed.trait.data,foo))
  if(length(tmp)==0){
    contra<-contra_var<-vector('numeric',0)
  }else{
    contra<-array(tmp[1,])
    contra_var<-array(tmp[2,])
  }
  which_mis_SE<-array(which(p_SE==-1))
  n_mis_SE<-length(which_mis_SE)
  
  
  ##FORM CALL##
  out.trait.se<-as.matrix(p_SE)
  out.trait.se[out.trait.se==-1]<-NA
  out.trait.se[out.trait.se==-2]<-Inf
  out.trait.se<-sqrt(out.trait.se)
  call<-list(trait.data=trait.data,trait.se=out.trait.se,tree=tree)
  
  
  ##TRANSFORMATION CONSTANTS##
  trans.const<-.get.trans.const(tree,X,mis)
  tree$edge.length<-tree$edge.length/trans.const$hgt
  for(i in 1:ncol(X)){
    X[,i]<-X[,i]/sqrt(trans.const$X_sig2[i])
  }
  #not prepped for multivariate data
  p_SE[p_SE>0]<-p_SE[p_SE>0]/trans.const$X_sig2[1]
  if(length(contra)>1){
    contra<-contra/sqrt(trans.const$X_sig2[1])
  }
  if(!sampling.scale){
    for(i in c("Ysig2_prior_mu","Ysig2_prior_sig")){
      if(!is.null(prior.list[[i]])){
        prior.list[[i]]<-prior.list[[i]]/trans.const$X_sig2
      }
    }
    if(!is.null(prior.list$R0_prior_mu)){
      prior.list$R0_prior_mu<-prior.list$R0_prior_mu+log(trans.const$hgt)-log(mean(trans.const$X_sig2))
    }
    for(i in c('Rsig2_prior_sig','Rmu_prior_mu','Rmu_prior_sig')){
      if(!is.null(prior.list[[i]])){
        prior.list[[i]]<-prior.list[[i]]*trans.const$hgt
      }
    }
  }
  
  
  ##TREE INFO##
  edge.key<-which(tree$edge.length!=0)
  call$edge.key<-edge.key
  poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
  if(length(poly.nodes)>0){
    d_poly<-lapply(poly.nodes,function(ii) which(ii==tree$edge[,1]))
    tmp<-sort(unlist(lapply(d_poly,function(ii) ii[seq(2,length(ii)-1)])))
    tmp<-tmp+0:(length(tmp)-1)
    real_e<-(1:(2*n-2))[-tmp]+1
    tree<-multi2di(tree,random=FALSE)
  }else{
    real_e<-1:e+1
  }
  des_e<-sapply(1:nrow(tree$edge),function(ii) which(tree$edge[,1]==tree$edge[ii,2])+1)
  tip_e<-which(lengths(des_e)==0)
  tip_e<-tip_e[order(tree$edge[tip_e,2])]
  des_e[tip_e]<-list(rep(-1,2))
  des_e<-matrix(unlist(des_e),ncol=2,byrow=TRUE)
  root.edges<-which(tree$edge[,1]==n+1)+1
  des_e<-rbind(root.edges,des_e)
  prune_T<-c(0,tree$edge.length)
  which_0tip<-prune_T[real_e]==0
  eV<-edge.vcv(tree)
  if(sum(which_0tip)>0){
    tree<-drop.tip(tree,tree$tip.label[tree$edge[(real_e-1)[which_0tip],2]],trim.internal=FALSE,collapse.singles=FALSE)
    real_e<-real_e[!which_0tip]
    e<-length(real_e)
  }
  eV<-eV[real_e-1,real_e-1]
  #test--below now only works for polytomies, but not necessarily 0 tips
  # edge.quants<-sample(10,e,replace=TRUE)
  # cols=rainbow(10)
  # plot(call$tree,edge.color=cols[edge.quants],show.tip.label=FALSE)
  # vec<-rep('gray',nrow(tree$edge))
  # vec[real_e[edge.key]-1]<-cols[edge.quants]
  # plot(tree,edge.color=vec,show.tip.label=FALSE)
  #seems to work!
  postorder<-((2*n-2):1)[-((2*n-2)-tip_e)]
  tip_e<-tip_e+1
  XX<-vector('logical',2*n-1)
  XX[tip_e]<-ifelse(p_SE==-2,FALSE,TRUE)
  mis_code<-vector('numeric',n-1)
  which_non_mis<-vector('list',n-1)
  counter<-0
  for(i in postorder){
    des_X<-XX[des_e[i,]]
    counter<-counter+1
    mis_code[counter]<-sum(des_X)
    if(mis_code[counter]==0){
      XX[i]<-FALSE
      which_non_mis[[counter]]<-numeric(0)
    }else if(mis_code[counter]==1){
      XX[i]<-TRUE
      which_non_mis[[counter]]<-which(des_X)
    }else if(mis_code[counter]==2){
      XX[i]<-TRUE
      which_non_mis[[counter]]<-c(1,2)
    }
  }
  n_mis_2<-sum(mis_code==2)
  which_non_mis<-unlist(which_non_mis[mis_code==1])
  if(is.null(which_non_mis)){
    which_non_mis<-vector('integer',0)
  }
  n_mis_1<-length(which_non_mis)
  
  
  ##PUTTING IT ALL TOGETHER##
  dat<-list('n'=n,'e'=e,'eV'=eV,'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'postorder'=postorder,
            'mis_code'=mis_code,'n_mis_2'=n_mis_2,'n_mis_1'=n_mis_1,'which_non_mis'=array(which_non_mis),
            'X'=as.vector(X),'p_SE'=p_SE,'inv_n_obs'=inv_n_obs,'n_contra'=n_contra,'contra'=contra,'contra_var'=contra_var,
            'n_mis_SE'=n_mis_SE,'which_mis_SE'=which_mis_SE)
  dat<-c(dat,prior.list)
  dat$constr_Rsig2<-as.numeric(constrain.Rsig2)
  dat$constr_Rmu<-as.numeric(!trend)
  if(lik.power<0){
    warning('lik.power must be 0 or positive: set to default of 1.',
            immediate.=TRUE)
    lik.power<-1
  }
  dat$lik_power<-lik.power
  list(call=call,trans.const=trans.const,dat=dat)
}

#' Run sampler for an Evolving Rates model
#'
#'
#' This function runs a Stan-based Hamiltonian Monte Carlo (HMC) sampler given input data and priors, but only
#' stores, rather than processes, the output.
#'
#'
#' @param input.evorates.obj The output of \code{input.evorates()}
#' @param return.as.obj \code{TRUE}or \code{FALSE}: should results be passed back to the R
#' environment as a \code{evorates_fit} object? Defaults to \code{TRUE}. If \code{FALSE}, results are saved to files
#' with automatically generated names instead (see below).
#' @param out.file A directory to save results to. If unspecified, an automatic directory is
#' generated based on the current date and time and R's working directory. The function will generate csv files
#' for each chain of the HMC sampler using Stan's built-in functionality (note that these will thus be on the
#' transformed scale), as well as a separate RDS file giving additional information about the data. Defaults to
#' \code{NULL}, which means no results are saved to file, though this will be changed automatically if
#' \code{return.as.obj = FALSE}.
#' @param check.overwrite \code{TRUE} or \code{FALSE}: should files in the directory specified by \code{out.file}
#' be checked to
#' prevent accidentally overwriting existing files? This part of the code is not thoroughly tested and might take
#' a long time for folders with many files, so some users may wish to just switch it off. Defaults to \code{TRUE}.
#' @param constrain.Rsig2 \code{TRUE} or \code{FALSE}: should the rate varinace (\code{R_sig2}) parameter be constrained to 0,
#' resulting in a simple Brownian Motion or early/late burst model? Defaults to \code{FALSE}. See Details for a
#' definition of model parameters.
#' @param trend \code{TRUE} or \code{FALSE}: should a trend in rates over time (\code{R_mu}) be estimated,
#' resulting in an early/late burst or trended evorates model? Defaults to \code{FALSE}.
#' @param lik.power A single number between 0 and 1 specifying what power to raise the likelihood function to.
#' This is useful for sampling from "power posteriors" that shift the posterior to look more like the prior (indeed,
#' you can set this to 0 to sample from the prior itself). Useful for model diagnostics and calculating things
#' like Bayes Factors. Technically, you can set \code{lik.power} above 1 (a technique called "data cloning") but
#' this is not beneficial for evorates models and we don't recommend it.
#' @param ... Other optional arguments:
#' \itemize{
#' \item{Prior arguments (see details for further information on what parameters mean):
#' \itemize{
#' \item{In general, prior standard deviations and, in some cases, means, can be tweaked by passing arguments
#'       named "<parameter_name>_sd" and "<parameter_name>_mean", respectively. For example:
#'       \code{R_mu_mean = 1} or \code{R_sig2_sd = 1}. Note that all prior arguments are necessarily on the
#'       sampling scale when passed to \code{run.evorates()}}
#' \item{Prior on rate variance (\code{R_sig2}): Follows a half-Cauchy distribution (basically a half-normal distribution 
#'       with extremely fat tails) with adjustable standard deviation. By default, standard deviation is set to 5 divided
#'       by the maximum height of \code{tree} (i.e., 5 on the sampling scale).}
#' \item{Prior on trend (\code{R_mu}): Follows a normal distribution with adjustable mean and standard deviation. By
#'       default, mean is set to 0 (also 0 on sampling scale), and standard deviation is set to 10 divided by the
#'       maximum height of \code{tree} (i.e., 10 on the sampling scale).}
#' \item{Prior on rate at the root (\code{R_0}): Follows a normal distribution with adjustable mean and standard deviation.
#'       By default, mean is set to variance of \code{trait.data} divided by the maximum height of \code{tree} on the
#'       natural log scale (i.e., 0 on sampling scale), and standard deviation is set to 10 (also 10 on sampling scale).}
#' \item{Prior on tip error variance (\code{Y_sig2}): Follows a truncated t distribution (basically a normal
#'       distribution with potentially fatter tails) with adjustable mean, standard deviation, and degrees of freedom.
#'       The distribution is truncated at 0 to the left such that Y_sig2 estimates are always positive (note that setting
#'       the mean to 0 results in a half-t distribution). By default, mean is set to 0, standard deviation is set to 2
#'       times the variance of \code{trait.data} (i.e., 2 on the sampling scale), and degrees of freedom is set to 1. The
#'       degrees of freedom can be tweaked with the argument "Ysig2_df" and may be any positive number, including \code{Inf},
#'       with lower numbers leading to a less informative priors with fatter tails. Notably, setting "Ysig2_df" to \code{Inf}
#'       makes this prior a truncated normal and much more informative!}
#' }}
#' \item{Additional arguments to pass to \code{rstan::sampling()}, most commonly:
#' \itemize{
#' \item{\code{chains} to specify the number of HMC chains (defaults to 4)}
#' \item{\code{iter} to  specify the number of iterations in HMC chains (defaults to 2000)}
#' \item{\code{warmup} to specify the number of warmup iterations (defaults to \code{floor(iter/2)})}
#' \item{\code{thin} to specify which iterations to keep in results (defaults to 1 or no thinning)}
#' \item{\code{cores} to specify the number of computer cores to use (defaults to \code{getOption("mc.cores", 1L)})}
#' \item{\code{refresh} to control when progress is reported (defaults to \code{max(iter/10, 1)}, and can be
#' suppressed by setting to 0 or less)}
#' \item{There are other things users might want to mess with, like \code{seed}, \code{init}, and \code{control}.
#' See \code{?rstan::sampling} and \code{?rstan::stan} for further details.}
#' }}
#' }
#' 
#' 
#' @return A list if \code{return.as.obj = TRUE}, containing the unprocessed model fit in \code{stanfit} format and
#' updated \code{input.evorates.obj}. Otherwise, the directory results were saved to (see \code{out.file} for details).
#' 
#' 
#' #' @details Parameter definitions:
#' \itemize{
#' \item{Rate variance (\code{R_sig2}, also denoted \eqn{\sigma^2_{\sigma^2}}{\sigma^2[\sigma^2]}): Determines how much random variation
#'       accumulates in rates over time. Specifically, independent lineages evolving for a length of time \code{t} would exhibit log-normally
#'       distributed rates with standard deviation \code{sqrt(t * R_sig2)}. Another way to think about this is that the 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time is given by \code{exp(c(-1,1) * 1.96 * sqrt(R_sig2))}, assuming \code{R_mu = 0.}}
#' \item{Trend (\code{R_mu}, also denoted \eqn{\mu_{\sigma^2}}{\mu[\sigma^2]}): Determines whether median rates tend to decrease (if
#'       negative) or increase (if positive) over time. Specifically, the median fold-change in rate for a given lineage is 
#'       \code{exp(t * R_mu)} over a length of time \code{t}. This is distinct from changes in average rates, which depends on both
#'       \code{R_mu} and \code{R_sig2}.}
#' \item{Note on combining rate variance and trend parameters: The 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time, given a trend, is simply
#'       \code{exp(R_mu) * exp(c(-1,1) * 1.96 * sqrt(R_sig2))} or equivalently \code{exp(R_mu + c(-1,1) * 1.96 * sqrt(R_sig2))}.
#'       This distribution of rate change is right-skewed due to the \code{exp()} function. Because of this, median changes in rates
#'       over time will always be lower than average changes. The "average change" parameter, denoted \code{R_del} or
#'       \eqn{\delta_{\sigma^2}}{\delta[\sigma^2]}, is given by \code{R_mu + R_sig2 / 2}. Accordingly, the average fold-change in rate
#'       for a given lineage is \code{exp(t * R_del)} over a length of time \code{t}.}
#' \item{Tip error variance (\code{Y_sig2}, also denoted \eqn{\sigma^2_y}{\sigma^2[y]}): Determines the among of "error" around observed
#'       trait values for tips without fixed standard errors. Specifically, raw observations for a given tip are sampled from a normal
#'       distribution centered at that tip's "true trait value" with standard deviation \code{sqrt(Y_sig2)}.}
#' \item{Rate at the root (\code{R_0}, also denoted \eqn{\sigma^2_0}{\sigma^2[0]}): The natural log of the rate at the root of the entire phylogeny.}
#' \item{Branchwise rates (\code{R_i}, also denoted \eqn{ln \bar{\sigma^2_i}}{ln \sigma^2[i]}, where i is the index of an edge in \code{tree}): 
#'       The natural log of the average rate along branches of the phylogeny. Note that rates are always shifting over time under this model,
#'       and these quantities are thus averages. The true rate value at any particular time point along a branch is a related, but separate,
#'       quantity. These will be NA for branches of length 0.}
#' \item{Background rate (\code{bg_rate}): The natural log of the average trait evolution rate for the entire phylogeny, given by the average of
#'       \code{exp(R_i)}, with each entry weighted by its respective branch length. NAs are ignored in this calculation since they correspond
#'       to branches of length 0.}
#' \item{Branchwise rate deviations (\code{Rdev_i}, also denoted ln \eqn{\bar{\sigma^2_{dev,i}}}{ln \sigma^2[dev,i]}): Determines whether a
#'       branches exhibit relatively "fast" (if positive) or "slow" rates (if negative). Generally, these are the differences
#'       between the branchwise rates and \emph{geometric} background rate on the natural log scale,
#'       though this depends on \code{remove.trend)}. Here, the geometric background rate is defined as the weighted average of \code{R_i},
#'       with weights corresponding to branch lengths. This helps prevent some issues with comparing highly right-skewed distributions.
#'       If \code{remove.trend = TRUE}, then branchwise rates are first "detrended" prior to calculating background rates and deviations
#'       (\code{R_i - (-log(abs(R_mu)) - log(l_i) + log(abs(exp(R_mu * t1_i) - exp(R_mu * t0_i))))}, where \code{l} is a vector of
#'       branch lengths and \code{t0}/\code{t1} are vectors of start and end times of each branch). This basically make branchwise rate deviations
#'       determine whether branches exhibit slow/fast rates \emph{given} the overall trend in rates through time. Otherwise, branchwise rate
#'       deviations will simmply indicate slow/fast branches occur at the root/tips of a tree in the presence of a strong trend.}
#' }
#' 
#' 
#' @seealso \link{fit.evorates} is a convenient wrapper for \link{input.evorates}, \link{run.evorates}, and \link{output.evorates}
#' 
#' 
#' @examples
#' #get whale/dolphin tree/trait data
#' data("cet_fit")
#' tree <- cet_fit$call$tree
#' trait.data <- cet_fit$call$trait.data
#' 
#' #prepare data for evorates model
#' input <- input.evorates(tree = tree, trait.data = trait.data)
#' #run the evorates model sampler (takes a couple minutes)
#' run <- run.evorates(input.evorates.obj = input, out.file = "test_fit", chains = 1)
#' #process the output from the object
#' fit <- output.evorates(run.evorates.obj = run)
#' #process the output from the files
#' fit <- output.evorates(run.evorates.obj = "test_fit")
#' 
#' #see fit.evorates() documentation for further examples
#' 
#' 
#' @export
run.evorates<-function(input.evorates.obj,return.as.obj=TRUE,out.file=NULL,check.overwrite=TRUE,
                       constrain.Rsig2=FALSE,trend=FALSE,lik.power=1,...){
  stan.args.list<-list(...)
  if(length(stan.args.list)>0){
    prior.args.names<-c('^(R[\\._]?0).*((mu)|(mean))','(R[\\._]?0).*((sig)|(sd))',
                        '^((intra)|(Y))[\\._]?((sig2)|(var))',
                        '^R[\\._]?((sig2)|(var))',
                        '^((R[\\._]?mu)|(trend)).*((mu)|(mean))','^((R[\\._]?mu)|(trend)).*((sig)|(sd))')
    which.prior.args<-apply(do.call(cbind,lapply(prior.args.names,grepl,x=names(stan.args.list))),1,any)
    stan.args.list<-stan.args.list[!which.prior.args]
    stan.args.list<-stan.args.list[!(names(stan.args.list)%in%
                                       c('object','data','pars','include','sample_file','chain_id','save_warmup'))]
  }
  if(!is.null(stan.args.list$chains)){
    nchain<-stan.args.list$chains
  }else{
    nchain<-4
  }
  prior.list<-do.call(.make.prior.list,list(...))
  prep<-.prep.run(input.evorates.obj,return.as.obj,out.file,check.overwrite,
                  constrain.Rsig2,trend,lik.power,prior.list,nchain)
  if(return.as.obj){
    ret<-do.call(sampling,
                 c(object=list(stanmodels[[prep$stanobj]]),
                   data=list(prep$input.evorates.obj$dat),
                   pars=list(prep$exclude.pars),
                   include=list(FALSE),
                   sample_file=list(prep$out.file),
                   stan.args.list))
  }else{
    do.call(sampling,
            c(object=list(stanmodels[[prep$stanobj]]),
              data=list(prep$input.evorates.obj$dat),
              pars=list(prep$exclude.pars),
              include=list(FALSE),
              sample_file=list(prep$out.file),
              stan.args.list))
  }
  if(!is.null(prep$out.file)){
    saveRDS(prep$input.evorates.obj,prep$file.strings[length(prep$file.strings)])
    message('HMC samples were saved as:\n\t',
            paste(basename(prep$file.strings[-length(prep$file.strings)]),'\n',collapse='\t'),
            'and input data/transformation constants were saved as:\n\t',
            basename(prep$file.strings[length(prep$file.strings)]),'\n',
            'in ',normalizePath(prep$directory,winslash='/'))
    #return output file basename
    out<-gsub('_\\d+\\.csv$','',prep$file.strings[1])
  }
  if(return.as.obj){
    out<-c(stanfit=ret,prep$input.evorates.obj)
  }
  out
}

#' Process the output of an Evolving Rates model
#'
#'
#' This function processes raw Stan-based Hamiltonian Monte Carlo (HMC) samples and associated information and returns
#' all this information in a (relatively) user-friendly format.
#'
#'
#' @param run.evorates.obj Either the output of \code{run.evorates()} or a string specifying files to load. These
#' files always consist of some base name followed by a variable suffix: either "_<X>.csv", where <X> is a chain
#' number, or "_info". These files contain Stan-based HMC samples/information and auxiliary information related to
#' data input/transformation constants, respectively. In general, the string should consist of \emph{just} the
#' base name, not the variable suffix, though the function does attempt to remove suffixes if detected.
#' @param stanfit,call,trans.const,dat Generally, you should \emph{not mess with these} unless you know what you're doing!
#' These are the basic components the function uses to form its output. \code{stanfit} is the Stan-based HMC
#' samples/information, while the rest make up the auxiliary information regarding input and transformation constants.
#' You can use these arguments to overwrite certain parts of \code{run.evorates.obj} for debugging purposes.
#' @param include.warmup \code{TRUE} or \code{FALSE}: should warmup be included in posterior samples? Warmup is
#' always included for parameters used to tune the HMC chain, but warmup may be included or excluded for actual
#' estimated parameters. Defaults to \code{FALSE}.
#' @param report.quantiles A vector of posterior distribution quantiles to return (should be between 0 and 1). Set to
#' \code{NULL} to not return any quantiles. Defaults to 2.5\%, 50\%, and 97.5\% quantiles.
#' @param report.means \code{TRUE} or \code{FALSE}: should posterior distribution means be returned? Defaults to
#' \code{TRUE}.
#' @param report.devs \code{TRUE} or \code{FALSE}: should the difference between branchwise rates 
#' and the overall "background rate" on the natural log scale ("rate deviations") be returned?
#' The background rate is simply the mean branchwise rate, weighted by each branch's respective length.
#' If \code{TRUE}, also calculates the posterior probability rate deviations are greater than 0. These
#' additional parameters help give a sense of which branches in \code{tree} exhibit anomalous trait evolution rates.
#' Defaults to \code{TRUE}, but is automatically switched to \code{FALSE} when fitting a Brownian
#' Motion model (\code{constrain.Rsig2 = FALSE & trend = FALSE}).
#' @param remove.trend \code{TRUE} or \code{FALSE}: should the rate deviation calculations remove (i.e., "account for") the
#' effect of trends in rates over time? This may be helpful since strong trends in rates can mask otherwise anomalously
#' high or low trait evolution rates in certain parts of the tree. Defaults to \code{TRUE}, but has no effect if no trend was
#' fitted or if \code{report.devs} is \code{FALSE}. If \code{TRUE}, note that \code{report.devs} will be switched to \code{FALSE}
#' if fitting early/late burst model (\code{constrain.Rsig2 = FALSE & trend = TRUE}).
#' 
#' 
#' @return An object of class "\code{evorates_fit}" if \code{return.as.obj = TRUE}. Otherwise, the directory
#' results were saved to (see \code{out.file} for details). An \code{evorates_fit} object is a list of at least
#' 5 components:
#' \itemize{
#' \item{\code{call}, which contains information on the final tree, trait values, trait standard errors, and prior
#' parameters passed to Stan's HMC sampling algorithm (on untransformed scale for better interpretability,
#' see \code{sampling.scale}).}
#' \item{\code{sampler.control}, which contains various information on the HMC run, including the number of chains,
#' iterations, warmup, thinning rate, etc.}
#' \item{\code{sampler.params}, an array of class "\code{param_block}", containing parameters/diagnostics that were used to tune
#' the behavior of the HMC while Stan ran it, as well as the (log) prior (\code{prior__}), likelihood (\code{lik__}),
#' and posterior probability (\code{post__}) of each iteration in the HMC. See Stan manual for more information on
#' what the parameters mean. The likelihood is not raised to \code{lik.power} here, but the log posterior probability
#' is calculated while accounting for \code{lik.power}. This always includes warmup iterations, though these can be
#' discarded using \code{exclude.warmup()} or \code{combine.chains()}.}
#' \item{\code{param.diags}, a \code{param_block} array, containing diagnostics for each parameter
#' estimated during the fit, including the
#' initial value of the HMC chain (\code{init}), the bulk effective sample size (\code{bulk_ess}), the tail
#' effective sample size (\code{tail_ess}), and the Rhat (\code{Rhat}). See \code{?rstan::Rhat} for more details
#' on what these diagnostics mean. Generally, you want effective sample sizes to be greater than 400ish and Rhat to be
#' less than 1.01ish. The functions \code{check.ess()} and \code{check.mix()} will check these thresholds for you
#' automatically.}
#' \item{\code{chains}, another \code{param_block} array containing sampled parameter values for each parameter estimated during the fit. See
#' details for further information on what each parameter means.}
#' \item{The object optionally contains more \code{param_block} arrays of posterior distribution quantiles (\code{quantiles})
#' and means (\code{means}), and posterior probabilities (\code{post.prob}, see \code{report.devs} and \code{remove.trend}).}
#' }
#' All \code{param_block} arrays' dimensions go in the order of iterations/diagnostics/quantiles, then parameters, then chains.
#' 
#' 
#' @details Parameter definitions:
#' \itemize{
#' \item{Rate variance (\code{R_sig2}, also denoted \eqn{\sigma^2_{\sigma^2}}{\sigma^2[\sigma^2]}): Determines how much random variation
#'       accumulates in rates over time. Specifically, independent lineages evolving for a length of time \code{t} would exhibit log-normally
#'       distributed rates with standard deviation \code{sqrt(t * R_sig2)}. Another way to think about this is that the 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time is given by \code{exp(c(-1,1) * 1.96 * sqrt(R_sig2))}, assuming \code{R_mu = 0.}}
#' \item{Trend (\code{R_mu}, also denoted \eqn{\mu_{\sigma^2}}{\mu[\sigma^2]}): Determines whether median rates tend to decrease (if
#'       negative) or increase (if positive) over time. Specifically, the median fold-change in rate for a given lineage is 
#'       \code{exp(t * R_mu)} over a length of time \code{t}. This is distinct from changes in average rates, which depends on both
#'       \code{R_mu} and \code{R_sig2}.}
#' \item{Note on combining rate variance and trend parameters: The 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time, given a trend, is simply
#'       \code{exp(R_mu) * exp(c(-1,1) * 1.96 * sqrt(R_sig2))} or equivalently \code{exp(R_mu + c(-1,1) * 1.96 * sqrt(R_sig2))}.
#'       This distribution of rate change is right-skewed due to the \code{exp()} function. Because of this, median changes in rates
#'       over time will always be lower than average changes. The "average change" parameter, denoted \code{R_del} or
#'       \eqn{\delta_{\sigma^2}}{\delta[\sigma^2]}, is given by \code{R_mu + R_sig2 / 2}. Accordingly, the average fold-change in rate
#'       for a given lineage is \code{exp(t * R_del)} over a length of time \code{t}.}
#' \item{Tip error variance (\code{Y_sig2}, also denoted \eqn{\sigma^2_y}{\sigma^2[y]}): Determines the among of "error" around observed
#'       trait values for tips without fixed standard errors. Specifically, raw observations for a given tip are sampled from a normal
#'       distribution centered at that tip's "true trait value" with standard deviation \code{sqrt(Y_sig2)}.}
#' \item{Rate at the root (\code{R_0}, also denoted \eqn{\sigma^2_0}{\sigma^2[0]}): The natural log of the rate at the root of the entire phylogeny.}
#' \item{Branchwise rates (\code{R_i}, also denoted \eqn{ln \bar{\sigma^2_i}}{ln \sigma^2[i]}, where i is the index of an edge in \code{tree}): 
#'       The natural log of the average rate along branches of the phylogeny. Note that rates are always shifting over time under this model,
#'       and these quantities are thus averages. The true rate value at any particular time point along a branch is a related, but separate,
#'       quantity. These will be NA for branches of length 0.}
#' \item{Background rate (\code{bg_rate}): The natural log of the average trait evolution rate for the entire phylogeny, given by the average of
#'       \code{exp(R_i)}, with each entry weighted by its respective branch length. NAs are ignored in this calculation since they correspond
#'       to branches of length 0.}
#' \item{Branchwise rate deviations (\code{Rdev_i}, also denoted ln \eqn{\bar{\sigma^2_{dev,i}}}{ln \sigma^2[dev,i]}): Determines whether a
#'       branches exhibit relatively "fast" (if positive) or "slow" rates (if negative). Generally, these are the differences
#'       between the branchwise rates and \emph{geometric} background rate on the natural log scale,
#'       though this depends on \code{remove.trend)}. Here, the geometric background rate is defined as the weighted average of \code{R_i},
#'       with weights corresponding to branch lengths. This helps prevent some issues with comparing highly right-skewed distributions.
#'       If \code{remove.trend = TRUE}, then branchwise rates are first "detrended" prior to calculating background rates and deviations
#'       (\code{R_i - (-log(abs(R_mu)) - log(l_i) + log(abs(exp(R_mu * t1_i) - exp(R_mu * t0_i))))}, where \code{l} is a vector of
#'       branch lengths and \code{t0}/\code{t1} are vectors of start and end times of each branch). This basically make branchwise rate deviations
#'       determine whether branches exhibit slow/fast rates \emph{given} the overall trend in rates through time. Otherwise, branchwise rate
#'       deviations will simmply indicate slow/fast branches occur at the root/tips of a tree in the presence of a strong trend.}
#' }
#' 
#' 
#' @seealso \link{fit.evorates} is a convenient wrapper for \link{input.evorates}, \link{run.evorates}, and \link{output.evorates}
#' 
#' 
#' @examples
#' #get whale/dolphin tree/trait data
#' data("cet_fit")
#' tree <- cet_fit$call$tree
#' trait.data <- cet_fit$call$trait.data
#' 
#' #prepare data for evorates model
#' input <- input.evorates(tree = tree, trait.data = trait.data)
#' #run the evorates model sampler (takes a couple minutes)
#' run <- run.evorates(input.evorates.obj = input, out.file = "test_fit", chains = 1)
#' #process the output from the object
#' fit <- output.evorates(run.evorates.obj = run)
#' #process the output from the files
#' fit <- output.evorates(run.evorates.obj = "test_fit")
#' 
#' #see fit.evorates() documentation for further examples
#' 
#' 
#' @export
output.evorates<-function(run.evorates.obj,stanfit=NULL,call=NULL,trans.const=NULL,dat=NULL,include.warmup=FALSE,
                          report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.devs=TRUE,remove.trend=TRUE){
  ##GATHER UP COMPONENTS##
  list.from.input<-list(stanfit=stanfit,call=call,trans.const=trans.const,dat=dat)
  list.from.run<-list(stanfit=NULL,call=NULL,trans.const=NULL,dat=NULL)
  #get any components specified by run.evorates object
  if(is.character(run.evorates.obj)){
    simple.name<-gsub('_\\d+\\.csv$|\\.csv$|_info$','',run.evorates.obj)
    directory<-dirname(simple.name)
    file.list<-list.files(directory)
    file.strings<-grep(paste0('^',basename(simple.name),'_\\d+\\.csv'),file.list,value=TRUE)
    if(length(file.strings)>0){
      file.strings<-paste0(simple.name,
                           substr(file.strings,regexpr('_\\d+\\.csv',file.strings),nchar(file.strings)))
    }else{
      file.strings<-NULL
    }
    file.strings<-c(file.strings,paste0(simple.name,'_info'))
    if(length(file.strings>1)){
      list.from.run$stanfit<-rstan::read_stan_csv(file.strings[-length(file.strings)])
    }
    info<-try(readRDS(file.strings[length(file.strings)]))
    if(!inherits(info,'try-error')){
      for(i in names(info)){
        if(!is.null(info[[i]])){
          list.from.run[[i]]<-info[[i]]
        }
      }
    }
  }else if(is.list(run.evorates.obj)){
    for(i in names(run.evorates.obj)){
      if(!is.null(run.evorates.obj[[i]])){
        list.from.run[[i]]<-run.evorates.obj[[i]]
      }
    }
  }
  #merge list from function input with list from evorates run object
  #function input overwrites evorates run object if it isn't NULL
  for(i in names(list.from.input)){
    if(!is.null(list.from.input[[i]])){
      list.from.run[[i]]<-list.from.input[[i]]
    }
  }
  stanfit<-list.from.run$stanfit
  call<-list.from.run$call
  trans.const<-list.from.run$trans.const
  dat<-list.from.run$dat
  if(is.null(stanfit)){
    stop('no MCMC samples were provided: please recheck provided filenames and/or R objects')
  }
  if(is.null(call)){
    stop('no information on the call (i.e., trait data/standard errors and tree) were provided: please recheck provided filenames and/or R objects')
  }
  if(is.null(trans.const)){
    trans.const<-.get.trans.const(call$tree,call$trait.data,call$trait.se)$trans.const
  }
  if(is.null(dat)){
    stop('no information on the inputted stan data were provided: please recheck provided filenames and/or R objects. This can be re-obtained by running input.evorates() with the SAME trait data, tree, and prior/lik.power settings!')
  }
  
  
  ##FORM CALL##
  checks<-setNames(vector('logical',3),c('Ysig2','Rsig2','Rmu'))
  for(i in names(checks)){
    tmp.check<-try(extract(stanfit,i),silent=TRUE)
    if(inherits(tmp.check,'try-error')){
      checks[i]<-FALSE
    }else{
      checks[i]<-TRUE
    }
  }
  has.intra<-checks[1]
  constrain.Rsig2<-!checks[2]
  trend<-checks[3]
  lik.power<-dat$lik_power
  if(constrain.Rsig2&report.devs){
    if(!trend){
      report.devs<-FALSE
      warning('report.devs was set to FALSE: rate deviations are meaningless with no rate heterogeneity',
              immediate.=TRUE)
    }else if(trend&remove.trend){
      report.devs<-FALSE
      warning('report.devs was set to FALSE: rate deviations are meaningless in the presence of only a trend if remove.trend is TRUE',
              immediate.=TRUE)
    }
  }
  out<-list()
  class(out)<-'evorates_fit'
  out$call<-c(call,
              'constrain.Rsig2'=list(constrain.Rsig2),
              'trend'=list(trend),
              'lik.power'=lik.power)
  
  
  ##SAMPLER ARGUMENTS/PARAMETERS##
  sampler.args<-stanfit@stan_args[[1]]
  sampler.args<-sampler.args[-which(names(sampler.args)=='chain_id')]
  sampler.args$chains<-length(stanfit@sim[[1]])
  nchain<-sampler.args$call.chains<-sampler.args$chains
  sampler.args$call.iter<-sampler.args$iter
  sampler.args$call.warmup<-sampler.args$warmup
  sampler.args$iter<-floor(sampler.args$iter/sampler.args$thin)
  sampler.args$warmup<-floor(sampler.args$warmup/sampler.args$thin)
  sampler.args$call.thin<-sampler.args$thin
  if(include.warmup){
    niter<-sampler.args$iter
    excl.iter<-niter+1
  }else{
    niter<-sampler.args$iter-sampler.args$warmup+1
    excl.iter<-seq_len(sampler.args$warmup)[-1]
  }
  sampler.params<-array(NA,c(sampler.args$iter,9,nchain),
                        list(iterations=NULL,
                             parameters=c(names(attr(stanfit@sim$samples[[1]],'sampler_params')),
                                          paste0(c('prior','lik','post'),'__')),
                             chains=paste('chain',1:nchain)))
  for(i in 1:nchain){
    sampler.params[,1:6,i]<-unlist(attr(stanfit@sim$samples[[i]],'sampler_params'))
  }
  sampler.params[,7,]<-extract(stanfit,"prior",permute=FALSE,inc_warmup=TRUE)
  sampler.params[,8,]<-extract(stanfit,"lik",permute=FALSE,inc_warmup=TRUE)
  sampler.params[,9,]<-sampler.params[,7,]+dat$lik_power*sampler.params[,8,]
  out$sampler.control<-sampler.args #I think this is excluding a bunch of things somehow...
  out$sampler.params<-sampler.params
  out$sampler.params<-.add.par.class(out$sampler.params)
  attr(out$sampler.params,'param_type')<-'chains'
  
  
  ##FORM CHAINS##
  e<-nrow(call$tree$edge)
  #will have to change back to some previous stuff for multivar extension
  out$chains<-array(dim=c(niter,5+e,nchain))
  #might consider changing param names for better intepretability...
  param.names<-c('R_0','Y_sig2','R_sig2','R_mu','bg_rate',paste('R_',1:e,sep=''))
  dimnames(out$chains)<-list(iterations=NULL,
                             parameters=param.names,
                             chains=paste('chain',1:nchain))
  R0<-extract(stanfit,"R0",permute=FALSE,inc_warmup=TRUE)[-excl.iter,,,drop=FALSE]-
    log(trans.const$hgt)+log(mean(trans.const$X_sig2))
  out$chains[,'R_0',]<-R0
  out$call$R0_prior_mu<-dat$R0_prior_mu-log(trans.const$hgt)+log(mean(trans.const$X_sig2))
  out$call$R0_prior_sig<-dat$R0_prior_sig
  if(has.intra){
    #will need to change for multivar
    Ysig2<-extract(stanfit,"Ysig2",permute=FALSE,inc_warmup=TRUE)[-excl.iter,,,drop=FALSE]*
      trans.const$X_sig2
    out$chains[,'Y_sig2',]<-Ysig2
    out$call$Ysig2_prior_sig<-dat$Ysig2_prior_sig*trans.const$X_sig2
    if(is.null(dat$Ysig2_prior_df)){
      #if this is a pre-t update evorates model, Ysig2 prior was cauchy and df should therefore be 1
      out$call$Ysig2_prior_df<-1
    }else{
      #A negative/0 df indicates Ysig2 prior was normal (infinite df)
      out$call$Ysig2_prior_df<-if(dat$Ysig2_prior_df>0) dat$Ysig2_prior_df else Inf
    }
    if(is.null(dat$Ysig2_prior_mu)){
      #if this is a pre-truncated update evorates model, Ysig2 prior was half-t and mu should therefore be 0
      out$call$Ysig2_prior_mu<-0
    }else{
      out$call$Ysig2_prior_mu<-dat$Ysig2_prior_mu*trans.const$X_sig2
    }
  }
  if(!constrain.Rsig2){
    Rsig2<-extract(stanfit,"Rsig2",permute=FALSE,inc_warmup=TRUE)[-excl.iter,,,drop=FALSE]/
      trans.const$hgt
    out$chains[,'R_sig2',]<-Rsig2
    out$call$Rsig2_prior_sig<-dat$Rsig2_prior_sig/trans.const$hgt
  }
  if(trend){
    Rmu<-extract(stanfit,"Rmu",permute=FALSE,inc_warmup=TRUE)[-excl.iter,,,drop=FALSE]/
      trans.const$hgt
    out$chains[,'R_mu',]<-Rmu
    out$call$Rmu_prior_mu<-dat$Rmu_prior_mu/trans.const$hgt
    out$call$Rmu_prior_sig<-dat$Rmu_prior_sig/trans.const$hgt
  }
  if(!constrain.Rsig2|trend){
    R<-extract(stanfit,"R",permute=FALSE,inc_warmup=TRUE)[-excl.iter,,,drop=FALSE]-
      log(trans.const$hgt)+log(mean(trans.const$X_sig2))
    wgts<-call$tree$edge.length[call$edge.key]/sum(call$tree$edge.length)
    bg.rate<-log(apply(R,c(1,2),function(ii) sum(exp(ii)*wgts)))
    out$chains[,'bg_rate',]<-bg.rate
    out$chains[,paste0('R_',call$edge.key),]<-aperm(R,c(1,3,2))
  }
  excl.inds<-apply(out$chains,2,function(ii) all(is.na(ii)))
  #ensure it includes all edge params if rate heterogeneity exists!
  if(!constrain.Rsig2|trend){
    excl.inds<-excl.inds&c(rep(TRUE,5),rep(FALSE,e))
  }
  out$chains<-out$chains[,!excl.inds,,drop=FALSE]
  out$chains<-.add.par.class(out$chains)
  attr(out$chains,'param_type')<-'chains'
  
  
  #add rate deviation chains
  if(report.devs){
    e.seq<-seq_len(Nedge(out))
    Rs<-.call.select(out$chains,paste0('R_',e.seq))
    el<-out$call$tree$edge.length
    if(trend&remove.trend){
      er<-edge.ranges(out)
      Rmu<-.call.select(out$chains,'R_mu')
      Rs<-Rs-(-log(abs(Rmu))-log(el)+log(abs(exp(Rmu*er[,2])-exp(Rmu*er[,1]))))
    }
    tmp<-Rs
    tmp[is.na(tmp)]<-0
    bg.rate<-sum(el/sum(el)*tmp)
    Rs<-Rs-bg.rate
    names(Rs)<-paste0('Rdev_',e.seq)
    out$chains<-.combine.par(list(out$chains,Rs))
    Rs[Rs==0]<-NA
    Rs<-Rs>0
    #add rate deviation posterior probabilities
    #should never get instances where deviations are perfectly 0...but just in case
    out$post.probs<-.call.op('means',list(chains=Rs,sampler.params=1),'.',FALSE)
    out$post.probs[is.na(out$post.probs)]<-0.5
  }
  
  
  ##EXTRA POSTERIOR DISTRIBUTION INFO##
  #add parameter diagnostics
  inits<-out$chains[1,,,drop=FALSE]
  if(!include.warmup){
    out$chains<-.call.select(out$chains,list(NULL,-1))
  }
  out$diagnostics<-.call.op('diagnostics',out,'.',FALSE)
  out$diagnostics['inits',,]<-inits

  #add quantiles
  if(!is.null(report.quantiles)){
    out$quantiles<-.call.op('quantiles',out,list('.',report.quantiles),FALSE)
  }
  #add means
  if(report.means){
    out$means<-.call.op('means',out,'.',FALSE)
  }
  
  
  out
}

#' Fit an Evolving Rates model
#'
#'
#' This function processes tree and trait data, runs a Stan-based Hamiltonian Monte Carlo (HMC) sampler to fit
#' these data to an evorates model, and returns the output of this sampler in a (relatively) user-friendly format.
#'
#'
#' @param tree An object of class "\code{phylo}"
#' @param trait.data Three options:
#' \itemize{
#' \item{A named vector of trait values.}
#' \item{A rownamed matrix of trait values. For the moment, must be 1 column since multivariate models are not
#' yet supported.}
#' \item{A data.frame with 2 columns: 1 with numeric data which is interpreted as trait values and 1 with
#' string/factor data interpreted as names. If more than 1 of either kind of column is found, returns error.
#' Will also use rownames if no string/factor column is found, but this limits data to consist of only 0-1
#' observations per tip.}
#' }
#' Both multiple observations for a single tip and missing observations are allowed.
#' In all cases, the associated names must generally match the tip labels found in \code{tree} (\code{tree$tip.label})
#' exactly. The only exception are numbers corresponding to internal node indices (these can be viewed by running:
#' \code{plot(tree); nodelabels()}), which are used to assign trait values to internal nodes. As of now, the function
#' actually completely ignores \code{tree$node.label}--I should probably fix that and plan to in the future. Be aware
#' that this feature slightly alters \code{tree} because it works by grafting extra tips of length 0 to the appropriate
#' nodes.
#' @param trait.se A vector, matrix, or data.frame of trait value standard errors which must be unambiguously 
#' labeled (see \code{trait.data}). Alternatively, a single, unlabeled number that will be applied to all tips in
#' \code{tree} (e.g., you could set it to 0 to specify all trait values are known without error). If \code{NULL}
#' (the default), standard error is estimated for all tips in the tree while fitting the model. There are several
#' things to note here: 
#' \itemize{
#' \item{Unlike \code{trait.data}, there can only be 1 trait standard error per tip.}
#' \item{You can use \code{NA} to specify tips that you want to estimate standard error for and \code{Inf} to
#' specify tips that have missing trait values.}
#' \item{Generally, tips with multiple observations and no observations are automatically assigned standard errors
#' of \code{NA} and \code{Inf}, respectively.}
#' \item{Any tips with unspecified standard error default to \code{NA}.}
#' }
#' Any conflicting/impossible standard error specifications are corrected and return warnings.
#' @param constrain.Rsig2 \code{TRUE} or \code{FALSE}: should the rate varinace (\code{R_sig2}) parameter be constrained to 0,
#' resulting in a simple Brownian Motion or early/late burst model? Defaults to \code{FALSE}. See Details for a
#' definition of model parameters.
#' @param trend \code{TRUE} or \code{FALSE}: should a trend in rates over time (\code{R_mu}) be estimated,
#' resulting in an early/late burst or trended evorates model? Defaults to \code{FALSE}.
#' @param lik.power A single number between 0 and 1 specifying what power to raise the likelihood function to.
#' This is useful for sampling from "power posteriors" that shift the posterior to look more like the prior (indeed,
#' you can set this to 0 to sample from the prior itself). Useful for model diagnostics and calculating things
#' like Bayes Factors. Technically, you can set \code{lik.power} above 1 (a technique called "data cloning") but
#' this is not beneficial for evorates models and we don't recommend it.
#' @param sampling.scale \code{TRUE} or \code{FALSE}: should provided prior parameters (see \code{...}) be
#' interpreted on the raw scale of the data or on the transformed scale? All data passed to the Stan-based HMC
#' sampler are transformed such that \code{tree}'s total height is 1 and the variance of the trait data
#' is 1. Defaults to \code{FALSE}, such that prior parameters are interpreted on the untransformed scale.
#' @param return.as.obj \code{TRUE}or \code{FALSE}: should results be passed back to the R
#' environment as a \code{evorates_fit} object? Defaults to \code{TRUE}. If \code{FALSE}, results are saved to files
#' with automatically generated names instead (see below).
#' @param out.file A directory to save results to. If unspecified, an automatic directory is
#' generated based on the current date and time and R's working directory. The function will generate csv files
#' for each chain of the HMC sampler using Stan's built-in functionality (note that these will thus be on the
#' transformed scale), as well as a separate RDS file giving additional information about the data. Defaults to
#' \code{NULL}, which means no results are saved to file, though this will be changed automatically if
#' \code{return.as.obj = FALSE}.
#' @param check.overwrite \code{TRUE} or \code{FALSE}: should files in the directory specified by \code{out.file}
#' be checked to
#' prevent accidentally overwriting existing files? This part of the code is not thoroughly tested and might take
#' a long time for folders with many files, so some users may wish to just switch it off. Defaults to \code{TRUE}.
#' @param include.warmup \code{TRUE} or \code{FALSE}: should warmup be included in posterior samples? Warmup is
#' always included for parameters used to tune the HMC chain, but warmup may be included or excluded for actual
#' estimated parameters. Defaults to \code{FALSE}.
#' @param report.quantiles A vector of posterior distribution quantiles to return (should be between 0 and 1). Set to
#' \code{NULL} to not return any quantiles. Defaults to 2.5\%, 50\%, and 97.5\% quantiles.
#' @param report.means \code{TRUE} or \code{FALSE}: should posterior distribution means be returned? Defaults to
#' \code{TRUE}.
#' @param report.devs \code{TRUE} or \code{FALSE}: should the difference between branchwise rates 
#' and the overall "background rate" on the natural log scale ("rate deviations") be returned?
#' The background rate is simply the mean branchwise rate, weighted by each branch's respective length.
#' If \code{TRUE}, also calculates the posterior probability rate deviations are greater than 0. These
#' additional parameters help give a sense of which branches in \code{tree} exhibit anomalous trait evolution rates.
#' Defaults to \code{TRUE}, but is automatically switched to \code{FALSE} when fitting a Brownian
#' Motion model (\code{constrain.Rsig2 = FALSE & trend = FALSE}).
#' @param remove.trend \code{TRUE} or \code{FALSE}: should the rate deviation calculations remove (i.e., "account for") the
#' effect of trends in rates over time? This may be helpful since strong trends in rates can mask otherwise anomalously
#' high or low trait evolution rates in certain parts of the tree. Defaults to \code{TRUE}, but has no effect if no trend was
#' fitted or if \code{report.devs} is \code{FALSE}. If \code{TRUE}, note that \code{report.devs} will be switched to \code{FALSE}
#' if fitting early/late burst model (\code{constrain.Rsig2 = FALSE & trend = TRUE}).
#' @param ... Other optional arguments:
#' \itemize{
#' \item{Prior arguments (see details for further information on what parameters mean):
#' \itemize{
#' \item{In general, prior standard deviations and, in some cases, means, can be tweaked by passing arguments
#'       named "<parameter_name>_sd" and "<parameter_name>_mean", respectively. For example:
#'       \code{R_mu_mean = 1} or \code{R_sig2_sd = 1}.}
#' \item{Prior on rate variance (\code{R_sig2}): Follows a half-Cauchy distribution (basically a half-normal distribution 
#'       with extremely fat tails) with adjustable standard deviation. By default, standard deviation is set to 5 divided
#'       by the maximum height of \code{tree} (i.e., 5 on the sampling scale).}
#' \item{Prior on trend (\code{R_mu}): Follows a normal distribution with adjustable mean and standard deviation. By
#'       default, mean is set to 0 (also 0 on sampling scale), and standard deviation is set to 10 divided by the
#'       maximum height of \code{tree} (i.e., 10 on the sampling scale).}
#' \item{Prior on rate at the root (\code{R_0}): Follows a normal distribution with adjustable mean and standard deviation.
#'       By default, mean is set to variance of \code{trait.data} divided by the maximum height of \code{tree} on the
#'       natural log scale (i.e., 0 on sampling scale), and standard deviation is set to 10 (also 10 on sampling scale).}
#' \item{Prior on tip error variance (\code{Y_sig2}): Follows a truncated t distribution (basically a normal
#'       distribution with potentially fatter tails) with adjustable mean, standard deviation, and degrees of freedom.
#'       The distribution is truncated at 0 to the left such that Y_sig2 estimates are always positive (note that setting
#'       the mean to 0 results in a half-t distribution). By default, mean is set to 0, standard deviation is set to 2
#'       times the variance of \code{trait.data} (i.e., 2 on the sampling scale), and degrees of freedom is set to 1. The
#'       degrees of freedom can be tweaked with the argument "Ysig2_df" and may be any positive number, including \code{Inf},
#'       with lower numbers leading to a less informative priors with fatter tails. Notably, setting "Ysig2_df" to \code{Inf}
#'       makes this prior a truncated normal and much more informative!}
#' }}
#' \item{Additional arguments to pass to \code{rstan::sampling()}, most commonly:
#' \itemize{
#' \item{\code{chains} to specify the number of HMC chains (defaults to 4)}
#' \item{\code{iter} to  specify the number of iterations in HMC chains (defaults to 2000)}
#' \item{\code{warmup} to specify the number of warmup iterations (defaults to \code{floor(iter/2)})}
#' \item{\code{thin} to specify which iterations to keep in results (defaults to 1 or no thinning)}
#' \item{\code{cores} to specify the number of computer cores to use (defaults to \code{getOption("mc.cores", 1L)})}
#' \item{\code{refresh} to control when progress is reported (defaults to \code{max(iter/10, 1)}, and can be
#' suppressed by setting to 0 or less)}
#' \item{There are other things users might want to mess with, like \code{seed}, \code{init}, and \code{control}.
#' See \code{?rstan::sampling} and \code{?rstan::stan} for further details.}
#' }}
#' }
#' 
#' 
#' @return An object of class "\code{evorates_fit}" if \code{return.as.obj = TRUE}. Otherwise, the directory
#' results were saved to (see \code{out.file} for details). An \code{evorates_fit} object is a list of at least
#' 5 components:
#' \itemize{
#' \item{\code{call}, which contains information on the final tree, trait values, trait standard errors, and prior
#' parameters passed to Stan's HMC sampling algorithm (on untransformed scale for better interpretability,
#' see \code{sampling.scale}).}
#' \item{\code{sampler.control}, which contains various information on the HMC run, including the number of chains,
#' iterations, warmup, thinning rate, etc.}
#' \item{\code{sampler.params}, an array of class "\code{param_block}", containing parameters/diagnostics that were used to tune
#' the behavior of the HMC while Stan ran it, as well as the (log) prior (\code{prior__}), likelihood (\code{lik__}),
#' and posterior probability (\code{post__}) of each iteration in the HMC. See Stan manual for more information on
#' what the parameters mean. The likelihood is not raised to \code{lik.power} here, but the log posterior probability
#' is calculated while accounting for \code{lik.power}. This always includes warmup iterations, though these can be
#' discarded using \code{exclude.warmup()} or \code{combine.chains()}.}
#' \item{\code{param.diags}, a \code{param_block} array, containing diagnostics for each parameter
#' estimated during the fit, including the
#' initial value of the HMC chain (\code{init}), the bulk effective sample size (\code{bulk_ess}), the tail
#' effective sample size (\code{tail_ess}), and the Rhat (\code{Rhat}). See \link[rstan]{Rhat} for more details
#' on what these diagnostics mean. Generally, you want effective sample sizes to be greater than 100 times the number of
#' chains and Rhat to be less than 1.01. The functions \code{check.ess()} and \code{check.mix()} will check these thresholds for you
#' automatically.}
#' \item{\code{chains}, another \code{param_block} array containing sampled parameter values for each parameter estimated during the fit. See
#' details for further information on what each parameter means.}
#' \item{The object optionally contains more \code{param_block} arrays of posterior distribution quantiles (\code{quantiles})
#' and means (\code{means}), and posterior probabilities (\code{post.prob}, see \code{report.devs} and \code{remove.trend}).}
#' }
#' All \code{param_block} arrays' dimensions go in the order of iterations/diagnostics/quantiles, then parameters, then chains.
#' 
#' 
#' @details Parameter definitions:
#' \itemize{
#' \item{Rate variance (\code{R_sig2}, also denoted \eqn{\sigma^2_{\sigma^2}}{\sigma^2[\sigma^2]}): Determines how much random variation
#'       accumulates in rates over time. Specifically, independent lineages evolving for a length of time \code{t} would exhibit log-normally
#'       distributed rates with standard deviation \code{sqrt(t * R_sig2)}. Another way to think about this is that the 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time is given by \code{exp(c(-1,1) * 1.96 * sqrt(R_sig2))}, assuming \code{R_mu = 0.}}
#' \item{Trend (\code{R_mu}, also denoted \eqn{\mu_{\sigma^2}}{\mu[\sigma^2]}): Determines whether median rates tend to decrease (if
#'       negative) or increase (if positive) over time. Specifically, the median fold-change in rate for a given lineage is 
#'       \code{exp(t * R_mu)} over a length of time \code{t}. This is distinct from changes in average rates, which depends on both
#'       \code{R_mu} and \code{R_sig2}.}
#' \item{Note on combining rate variance and trend parameters: The 95\% credible interval 
#'       for fold-changes in rates over 1 unit of time, given a trend, is simply
#'       \code{exp(R_mu) * exp(c(-1,1) * 1.96 * sqrt(R_sig2))} or equivalently \code{exp(R_mu + c(-1,1) * 1.96 * sqrt(R_sig2))}.
#'       This distribution of rate change is right-skewed due to the \code{exp()} function. Because of this, median changes in rates
#'       over time will always be lower than average changes. The "average change" parameter, denoted \code{R_del} or
#'       \eqn{\delta_{\sigma^2}}{\delta[\sigma^2]}, is given by \code{R_mu + R_sig2 / 2}. Accordingly, the average fold-change in rate
#'       for a given lineage is \code{exp(t * R_del)} over a length of time \code{t}.}
#' \item{Tip error variance (\code{Y_sig2}, also denoted \eqn{\sigma^2_y}{\sigma^2[y]}): Determines the among of "error" around observed
#'       trait values for tips without fixed standard errors. Specifically, raw observations for a given tip are sampled from a normal
#'       distribution centered at that tip's "true trait value" with standard deviation \code{sqrt(Y_sig2)}.}
#' \item{Rate at the root (\code{R_0}, also denoted \eqn{\sigma^2_0}{\sigma^2[0]}): The natural log of the rate at the root of the entire phylogeny.}
#' \item{Branchwise rates (\code{R_i}, also denoted \eqn{ln \bar{\sigma^2_i}}{ln \sigma^2[i]}, where i is the index of an edge in \code{tree}): 
#'       The natural log of the average rate along branches of the phylogeny. Note that rates are always shifting over time under this model,
#'       and these quantities are thus averages. The true rate value at any particular time point along a branch is a related, but separate,
#'       quantity. These will be NA for branches of length 0.}
#' \item{Background rate (\code{bg_rate}): The natural log of the average trait evolution rate for the entire phylogeny, given by the average of
#'       \code{exp(R_i)}, with each entry weighted by its respective branch length. NAs are ignored in this calculation since they correspond
#'       to branches of length 0.}
#' \item{Branchwise rate deviations (\code{Rdev_i}, also denoted ln \eqn{\bar{\sigma^2_{dev,i}}}{ln \sigma^2[dev,i]}): Determines whether a
#'       branches exhibit relatively "fast" (if positive) or "slow" rates (if negative). Generally, these are the differences
#'       between the branchwise rates and \emph{geometric} background rate on the natural log scale,
#'       though this depends on \code{remove.trend)}. Here, the geometric background rate is defined as the weighted average of \code{R_i},
#'       with weights corresponding to branch lengths. This helps prevent some issues with comparing highly right-skewed distributions.
#'       If \code{remove.trend = TRUE}, then branchwise rates are first "detrended" prior to calculating background rates and deviations
#'       (\code{R_i - (-log(abs(R_mu)) - log(l_i) + log(abs(exp(R_mu * t1_i) - exp(R_mu * t0_i))))}, where \code{l} is a vector of
#'       branch lengths and \code{t0}/\code{t1} are vectors of start and end times of each branch). This basically make branchwise rate deviations
#'       determine whether branches exhibit slow/fast rates \emph{given} the overall trend in rates through time. Otherwise, branchwise rate
#'       deviations will simmply indicate slow/fast branches occur at the root/tips of a tree in the presence of a strong trend.}
#' }
#' 
#' 
#' @seealso \link{fit.evorates} is a convenient wrapper for \link{input.evorates}, \link{run.evorates}, and \link{output.evorates}
#' 
#' 
#' @examples
#' #get whale/dolphin tree/trait data
#' data("cet_fit")
#' tree <- cet_fit$call$tree
#' trait.data <- cet_fit$call$trait.data
#' 
#' #fit data to evorates model (takes a couple minutes)
#' fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1)
#' 
#' #specifying parameter constraints
#' #add trend parameter
#' trend.fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1,
#'                           trend = TRUE)
#' #only trend parameter (early burst model)
#' EB.fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1,
#'                        trend = TRUE, constrain.Rsig2 = TRUE)
#' #no Rsig2 or trend (Brownian Motion model)
#' BM.fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1,
#'                        constrain.Rsig2 = TRUE)
#' 
#' #specifying trait standard error
#' #by default, estimates standard error, but you can set all tips to have no standard error
#' noSE.fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1,
#'                          trait.se = 0)
#' #or even set it for specific tips
#' SE.fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1
#'                        trait.se = setNames(c(1, 2), c("Orcinus_orca", "Balaenoptera_musculus")))
#'                        
#' #specifying priors
#' #make Rsig2 prior tighter (normally it's 5)
#' Risg2.fit <- fit.evorates(tree  = tree, trait.data = trait.data, chains = 1,
#'                           R_sig2_prior_sd = 1, sampling.scale = TRUE)
#' #maybe you REALLY believe rates should decrease by ~1% every million years?
#' Rmu.fit <- fit.evorates(tree = tree, trait.data = trait.data, chains = 1,
#'                         trend = TRUE,
#'                         R_mu_prior_mean = log(1 - 0.01), R_mu_prior_sd = log(1 + 0.005))
#' 
#' 
#' @export
fit.evorates<-function(tree,trait.data,trait.se=NULL,constrain.Rsig2=FALSE,trend=FALSE,lik.power=1,sampling.scale=FALSE,
                       return.as.obj=TRUE,out.file=NULL,check.overwrite=TRUE,
                       include.warmup=FALSE,
                       report.quantiles=c(0.025,0.5,0.975),report.means=TRUE,report.devs=TRUE,remove.trend=TRUE,
                       ...){
  input<-input.evorates(tree,trait.data,trait.se,constrain.Rsig2,trend,lik.power,sampling.scale,...)
  stan.args.list<-list(...)
  if(length(stan.args.list)>0){
    prior.args.names<-c('^(R[\\._]?0).*((mu)|(mean))','(R[\\._]?0).*((sig)|(sd))',
                        '^((intra)|(Y))[\\._]?((sig2)|(var))','^((intra)|(Y))[\\._]?((sig2)|(var)).*((nu)|(df))',
                        '^R[\\._]?((sig2)|(var))',
                        '^((R[\\._]?mu)|(trend)).*((mu)|(mean))','^((R[\\._]?mu)|(trend)).*((sig)|(sd))')
    which.prior.args<-apply(do.call(cbind,lapply(prior.args.names,grepl,x=names(stan.args.list))),1,any)
    stan.args.list<-stan.args.list[!which.prior.args]
    stan.args.list<-stan.args.list[!(names(stan.args.list)%in%
                                       c('object','data','pars','include','sample_file','chain_id','save_warmup'))]
  }
  if(!is.null(stan.args.list$chains)){
    nchain<-stan.args.list$chains
  }else{
    nchain<-4
  }
  prep<-.prep.run(input,return.as.obj,out.file,check.overwrite,
                  constrain.Rsig2,trend,lik.power,NULL,nchain)
  if(return.as.obj){
    ret<-do.call(sampling,
                 c(object=list(stanmodels[[prep$stanobj]]),
                   data=list(prep$input.evorates.obj$dat),
                   pars=list(prep$exclude.pars),
                   include=list(FALSE),
                   sample_file=list(prep$out.file),
                   stan.args.list))
  }else{
    do.call(sampling,
            c(object=list(stanmodels[[prep$stanobj]]),
              data=list(prep$input.evorates.obj$dat),
              pars=list(prep$exclude.pars),
              include=list(FALSE),
              sample_file=list(prep$out.file),
              stan.args.list))
  }
  if(!is.null(prep$out.file)){
    saveRDS(prep$input.evorates.obj,prep$file.strings[length(prep$file.strings)])
    message('MCMC samples were saved as:\n\t',
            paste(basename(prep$file.strings[-length(prep$file.strings)]),'\n',collapse='\t'),
            'and input data/transformation constants were saved as:\n\t',
            basename(prep$file.strings[length(prep$file.strings)]),'\n',
            'in ',normalizePath(prep$directory,winslash='/'))
    out<-gsub('_\\d+\\.csv$','',prep$file.strings[1])
  }
  if(return.as.obj){
    run<-c(stanfit=ret,prep$input.evorates.obj)
    out<-output.evorates(run,NULL,NULL,NULL,NULL,include.warmup,
                         report.quantiles,report.means,report.devs,remove.trend)
  }
  out
}