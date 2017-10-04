#function for message on load
.onAttach <- function(libname, pkgname){
  packageStartupMessage("HHG Ver. 2.2 - package for non parametric tests of independence and equality of distributions.")
  packageStartupMessage("type vignette(\'HHG\') or ?HHG for documentation, examples and a quickstart guide.")
  packageStartupMessage("use suppressPackageStartupMessages(library(HHG)) to suppress this message.")
}

.ml_by_variant = function(variant,mmin,mmax){
  m_param = NULL
  l_param = NULL
  if(variant %in% c('ADP','ADP-EQP','DDP')){
    m_param =  (mmin:mmax)
    l_param =  (mmin:mmax)
  }else if(variant %in% c('ADP-ML','ADP-EQP-ML')){
    for(i in (mmin:mmax)){
      for(j in (mmin:mmax)){
        m_param = c(m_param,i)
        l_param = c(l_param,j)
      }
    }
  }
  titles = paste0('m.',as.character(m_param),'X',as.character(l_param))
  ret=list()
  ret$m_param = m_param
  ret$l_param = l_param
  ret$titles = titles
  return(ret)
}


.check_compress_args = function(compress.p,compress.p0,compress.p1){
  if(is.null(compress.p) | is.null(compress.p0) | is.null(compress.p1)){
    stop('compress.p, compress.p0 and compress.p1 should be between 0 and 1, with p>=0.95, p0<=0.001 and p1<=0.0001')
  }
  if(is.na(compress.p) | is.na(compress.p0) | is.na(compress.p1)){
    stop('compress.p, compress.p0 and compress.p1 should be between 0 and 1, with p>=0.95, p0<=0.001 and p1<=0.0001')
  }
  if( compress.p <0.95 | compress.p > 1 | compress.p0 <0 | compress.p0>0.001 | compress.p1 <0 | compress.p1 > 0.0001 )
  {
    stop('compress.p, compress.p0 and compress.p1 should be between 0 and 1, with p>=0.95, p0<=0.001 and p1<=0.0001')
  }
}

# function for computing the statistic of the univariate df test of independence
hhg.univariate.ind.stat = function(x, y, variant = 'ADP',aggregation.type='sum',score.type='LikelihoodRatio', mmax = max(floor(sqrt(length(x))/2),2), mmin =2, w.sum = 0, w.max = 2,nr.atoms = nr_bins_equipartition(length(x))){
  correct.mi.bias = F # currently mi bias correction is not supported in interface
  type='Independence'
  flag_perform_ADP_multiple_partitions = F #flag used to check whether all partition sizes can be computed at single call
  
  #input checks:
  if(is.null(x) |is.null(y)){
    stop("x & y should be vectors of doubles." )
  }
  if(!(correct.mi.bias==TRUE || correct.mi.bias==FALSE) || w.sum!=as.integer(w.sum) || w.max!=as.integer(w.max) ){
    stop("Correct bias should be boolean, w.sum & w.max should be integers.")
  }
  if(length(x)!= length(y)){
    stop("X & Y not the same length")
  }
  if(!is.numeric(mmin) | !is.numeric(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin !=round(mmin) | mmax !=round(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin<2 | mmax<mmin){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  
  .hhg.univariate.check.inputs(type,variant,length(x),score.type,aggregation.type,nr.atoms = nr.atoms)
  max_not_available=TRUE
  if (((variant == 'DDP') && ((mmax <= 4)))) {
    max_not_available=FALSE
  }
  if (((variant %in% c('ADP','ADP-ML', 'ADP-EQP', 'ADP-EQP-ML')) && ((mmax <= 3)))) {
    max_not_available=FALSE
  }

  if(max_not_available == TRUE  & is.element(aggregation.type,c('max'))){
    stop(" Maximum based scores only available for ADP with m=2 and DDP with m=2,3,4")
  }
  if(max_not_available == TRUE  & is.element(aggregation.type,c('both'))){
    warning(" Maximum based scores only available for ADP with m=2 and DDP with m=2,3,4")
  }
  #if(is.element(variant,c('ADP-EQP','ADP-ML','ADP-EQP-ML')) & is.element(aggregation.type,c('max','both'))){
  #  stop(" Maximum based scores not available for 'ADP-EQP' , 'ADP-ML' or 'ADP-EQP-ML'. ")
  #}
  if(is.element(variant,c('ADP','ADP-EQP','ADP-ML','ADP-EQP-ML')) & correct.mi.bias == F & is.element(aggregation.type, c('sum','both'))){
    flag_perform_ADP_multiple_partitions = T
  }
  
  if(correct.mi.bias == T & is.element(variant,c('ADP-EQP','ADP-ML','ADP-EQP-ML'))){
    stop('The following variants do not support MI bias correction: ADP-EQP, ADP-ML and ADP-EQP-ML')
  }
  
  if(is.element(variant,c('ADP-EQP','ADP-EQP-ML')) & (nr.atoms<2 | nr.atoms>length(x))){
    stop('for ADP-EQP and ADP-EQP-ML variants, nr.atoms should be at least two and no more than sample size.')
  }
  
  if(is.element(variant,c('ADP-EQP','ADP-EQP-ML')) & mmax>nr.atoms){
    stop('mmax should be smaller or equal to nr.atoms')
  }
  
  if(nr.atoms!=round(nr.atoms)){
    stop('nr.atoms should be an integer')
  }
  do.nonregular = F
  do.nonregular.type = NA

  #check if we need to run maximum variants
  do.HHGRcpp.maximum.variants = F
  if(!max_not_available){
    if(is.element(variant,c('ADP','ADP-ML','ADP-EQP','ADP-EQP-ML')) & is.element(aggregation.type,c('max','both'))){
      do.HHGRcpp.maximum.variants  = T
    }  
  }
  skip.ADP.only.MAX = F
  if(is.element(variant,c('ADP','ADP-ML','ADP-EQP','ADP-EQP-ML')) & is.element(aggregation.type,c('max'))){
    skip.ADP.only.MAX = T
  }
  
  #statistics computation:
  stat.sl = rep(NA,(mmax)-mmin+1)
  stat.sc = rep(NA,(mmax)-mmin+1)
  stat.ml = rep(NA,(mmax)-mmin+1)
  stat.mc = rep(NA,(mmax)-mmin+1)
  
  titles = NULL
  if(flag_perform_ADP_multiple_partitions){
    m_param = NULL
    l_param = NULL
    eqp_param=NULL
    create_params = .ml_by_variant(variant,mmin,mmax)
    m_param = create_params$m_param
    l_param = create_params$l_param
    titles = create_params$titles
    if(variant %in% c('ADP','ADP-ML')){
      eqp_param = 0
    }else if(variant %in% c('ADP-EQP','ADP-EQP-ML')){
      eqp_param = 1
    }
    

    res = .hhg.test.adp_mk(x,y,M=m_param,L = l_param,equipartition.type = eqp_param,equipartitions.nr.bins = nr.atoms)
    stat.sl = res$sum.lr
    stat.sc = res$sum.chisq
    
  }else if(!skip.ADP.only.MAX){
    ms_vector=(mmin:mmax)
   for(i in 1:length(ms_vector)){
     m=ms_vector[i]
     #call by different variant an m to specific functions (some of which can compute max variants or more efficiently for small m)
     do.nonregular = F
     if (variant == 'DDP'){
       if (m != as.integer(m) || m < 2){
         stop('m, number of partitions,  must be an integer greater than 1')
       } else if (m == 2){
         if(w.sum==0 & w.max==2){
           x.variant = NA
           do.nonregular = T
           do.nonregular.type = 'DDP_2X2_OPTIMAL'
         }else{
           x.variant = 'spr.obs'  
         }
         
       } else if (m == 3){
         x.variant = 'ppr.33.obs'
       } else if (m == 4){
         x.variant = 'tpr.obs'
       } else {
         x.variant = 'ddp.obs'
       }
     } else if (variant == 'ADP'){
       if (m != as.integer(m) || m < 2){
         stop('m, number of partitions,  must be an integer greater than 1')
       } else if (m == 2) {
         x.variant = 'spr.all'
       } else if (m == 3) {
         # One could use 'ppr.33.all' that has the same complexity, but in practice it is much slower
         x.variant = 'ddp.all'
       } else if (m == 4) {
         # tpr.all is too time consuming
         x.variant = 'ddp.all'
       } else {
         x.variant = 'ddp.all'
       }
     }else{stop("Unkown Variant, should be \'ADP\' or \'DDP\' " )}
     
     if(!do.nonregular){
       ret_raw = .hhg.test.udfree(x = x, y = y, variant = x.variant, K = m, correct.mi.bias = correct.mi.bias,w.sum = w.sum,w.max = w.max)
       if (max_not_available) {
         ret_raw$max.chisq = NA
         ret_raw$max.lr    = NA
       }
       
       stat.sl[i] = ret_raw$sum.lr
       stat.sc[i] = ret_raw$sum.chisq
       stat.ml[i] = ret_raw$max.lr
       stat.mc[i] = ret_raw$max.chisq  
     }else{
       if(do.nonregular.type == 'DDP_2X2_OPTIMAL'){
         ret_raw = .hhg.ddp2.inversions(x,y)
         stat.sl[i] = ret_raw$sum.lr
         stat.sc[i] = ret_raw$sum.chisq
         stat.ml[i] = ret_raw$max.lr
         stat.mc[i] = ret_raw$max.chisq  
           
       }
     }
     
   } 
   titles = paste0('m.',as.character(mmin:mmax),'X',as.character(mmin:mmax))
  }
  
  if(do.HHGRcpp.maximum.variants){
    m_param = NULL
    l_param = NULL
    eqp_param=NULL
    create_params = .ml_by_variant(variant,mmin,mmax)
    m_param = create_params$m_param
    l_param = create_params$l_param
    titles = create_params$titles
    max_variant_nr_atoms = NA
    if(variant %in% c('ADP','ADP-ML')){
      max_variant_nr_atoms = length(x)
    }else if(variant %in% c('ADP-EQP','ADP-EQP-ML')){
      max_variant_nr_atoms = nr.atoms
    }
    for(i in 1:length(m_param)){
      res = NA
      if(m_param[i] == 2 && l_param[i] == 2){
        res = ADP_MAX_2X2_statistic(x,y,max_variant_nr_atoms,w.max = w.max)   
      }
      if(m_param[i] == 3 && l_param[i] == 2){
        res = ADP_MAX_3X2_statistic(x,y,max_variant_nr_atoms,w.max = w.max)   
      }
      if(m_param[i] == 2 && l_param[i] == 3){
        res = ADP_MAX_3X2_statistic(y,x,max_variant_nr_atoms,w.max = w.max)    
      }
      if(m_param[i] == 3 && l_param[i] == 3){
        res = ADP_MAX_3X3_statistic(x,y,max_variant_nr_atoms,w.max = w.max)   
      }
      
      stat.ml[i] = res$loglik.max
      stat.mc[i] = res$chisq.max
      
    }
    
  }
  
  flag_SL = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_SC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_ML = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  flag_MC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  
  #arrange output:
  if(flag_SL){names(stat.sl) = titles}
  if(flag_SC){names(stat.sc) = titles}
  if(flag_ML){names(stat.ml) = titles}
  if(flag_MC){names(stat.mc) = titles}

  ret=list()
  if(flag_SL + flag_SC + flag_ML + flag_MC==1){
      if(flag_SL){
        ret$statistic = stat.sl
      }else if(flag_SC){
        ret$statistic = stat.sc
      }else if(flag_ML){
        ret$statistic = stat.ml
      }else if(flag_MC){
        ret$statistic = stat.mc
      }
  }else{
    if(flag_SL){
      ret$sum.lr = stat.sl
    }
    if(flag_SC){
      ret$sum.chisq = stat.sc
    }
    if(flag_ML){
      ret$max.lr = stat.ml
    }
    if(flag_MC){
      ret$max.chisq = stat.mc
    }
  }
  ret$type = type
  ret$stat.type = 'Independence-Stat'
  ret$size = length(x)
  ret$variant = variant
  ret$score.type = score.type
  ret$aggregation.type = aggregation.type
  ret$mmin = mmin
  ret$mmax = mmax
  ret$nr.atoms = nr.atoms
  ret$additional = c(correct.mi.bias , w.sum,w.max)
  class(ret)='UnivariateStatistic'
  return (ret)
}

#function for performing validation checks on inputs. the logic of this function is being shared by many other functions,
#but other functions also have specific checks
.hhg.univariate.check.inputs=function(type,variant,size,score.type,aggregation.type,nr.atoms=NULL){
  if(is.null(variant)){
    stop("unknown variant - for independence should be 'ADP', 'DDP', 'ADP-EQP','ADP-ML' or 'ADP-EQP-ML'.\n For K-Sample should be 'KSample-Variant' or 'KSample-Equipartition'.")}
  
  if(is.null(score.type)){
    stop("score.type should be \'LikelihoodRatio\' , \'Pearson\'  or \'both\' ")
    
  }
  if(is.null(aggregation.type)){
    stop("aggregation.type should be \'sum\' , \'max\'  or \'both\' ")
  }
  if(type=='Independence'){
    if( length(size)>1  ){stop("For independence - group size should be single number")}
    
    if(!is.element(score.type,c('LikelihoodRatio','Pearson','both'))){
      stop("score.type should be \'LikelihoodRatio\' , \'Pearson\'  or \'both\' ")
    }
    
    if(!is.element(aggregation.type,c('sum','max','both'))){
      stop("aggregation.type should be \'sum\' , \'max\'  or \'both\' ")
    }
    
    if(!is.element(variant,c('DDP','ADP','ADP-EQP','ADP-ML','ADP-EQP-ML'))){
      stop("unknown variant - for independence should be 'ADP', 'DDP', 'ADP-EQP','ADP-ML' or 'ADP-EQP-ML'.")}
    
    if(is.element(variant,c('ADP-EQP','ADP-EQP-ML'))){
      if(is.null(nr.atoms)){stop('Nr.atoms must be an integer smaller than or equal to m.max')}
      if(is.na(nr.atoms)){stop('Nr.atoms must be an integer smaller than or equal to m.max')}
      if(nr.atoms != floor(nr.atoms)){stop('Nr.atoms must be an integer smaller than or equal to m.max')}
    }
    
  }else if(type=='KSample'){
    if(!is.element(variant,c('KSample-Variant','KSample-Equipartition'))){
      stop("unknown variant - for independence should be 'KSample-Variant' or 'KSample-Equipartition'.")}
    if(length(size)<2 ){
      stop("For k-sample - data should have at least two samples")}
    if(!is.element(score.type,c('LikelihoodRatio','Pearson','both'))){
      stop("score.type should be \'LikelihoodRatio\' , \'Pearson\'  or \'both\' ")
    }
    if(!is.element(aggregation.type,c('sum','max','both'))){
      stop("aggregation.type should be \'sum\' , \'max\'  or \'both\' ")
    }
    
    if(is.element(variant,c('KSample-Equipartition'))){
      if(is.null(nr.atoms)){stop('Nr.atoms must be an integer smaller than m.max')}
      if(is.na(nr.atoms)){stop('Nr.atoms must be an integer smaller than m.max')}
      if(nr.atoms != floor(nr.atoms)){stop('Nr.atoms must be an integer smaller than m.max')}
    }
    
  }else{
    stop('Unkown type, should be either \'Independence\' or \'KSample\' ')
  }
}

# MinP Object constructor from m.stats. Function constructs a univariate object, with ecdfs over the specific m distributions,
#but also ecdf for minp and fisher statistics with the specific m
#in order to handle ties in distributions, pvalues are computed using ecdf over the negative values (of single m stattistics, and Fisher)
# with an addition of a value of -Inf. MinP ecdf is attached a value of 0.
.hhg.univariate.object.constructor.from.mstats=function(m.stats,minm,maxm,type,variant,size,score.type,aggregation.type,nr.atoms,additional=NULL,ml_variant = F,compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001){
  .hhg.univariate.check.inputs(type,variant,size,score.type,aggregation.type,nr.atoms = nr.atoms)
  MESSAGES=F #flag used to show messages while computing the null distribution. currently has a value of F hardcoded.
  if(min(maxm)-minm+1 > dim(m.stats)[2] | min(maxm) <2  | min(maxm)<minm){
    stop("minm should be at least 2,  columns in m.stats should be the statistics for null table from, minm to maxm")
  }
  if(max(maxm)> nr.atoms){
    stop('Nr.atoms should be >= m.max')
  }
  if(type =='Independence'){
     #additional checks required? - not currently, all implemented in the checks function
    }
  
  if(type =='KSample'){
    if(sum(size) <min(maxm)){
      stop("size should be larger than maximum number of partitions")
    }
  }
  
  #mode switches - shouldn't be touched
  
  lookup.table=m.stats
  ret=list()
  if(ml_variant & length(maxm)>1){
    stop('Null table object does not support ADP-ML and ADP-EQP-ML with more than one mmax')
  }
  if(ml_variant & minm>2){
    stop('Null table object does not support ADP-ML and ADP-EQP-ML with mmin >2')
  }
  m.bound=max(maxm)
  nr.mmax=length(maxm)
  
  if(ml_variant){
    m.bound = (maxm-minm + 1)^2
  }
  
  # used to sort and count the ecdf by m (number of partitions)
  m.stat.negative.ecdf=list()
  m.stat.size = list()
  
  #These are matrices for holding the results for many computed m!
  MinP=matrix(rep(NA,nr.mmax*nrow(lookup.table)),nrow =nrow(lookup.table),ncol = nr.mmax)
  m_chosen=matrix(rep(NA,nr.mmax*nrow(lookup.table)),nrow =nrow(lookup.table),ncol = nr.mmax)
  Fisher=matrix(rep(0,nr.mmax*nrow(lookup.table)),nrow =nrow(lookup.table),ncol = nr.mmax)
  
  MinP_ECDF_LIST=list()
  FISHER_NEGATIVE_ECDF_LIST=list() #for higher alternative , one should use an ecdf on the negative to handle ties (conservativly).
  
  nr_cols = m.bound - minm + 1
  if(ml_variant){
    nr_cols = m.bound
  }
  p.val=matrix(rep(NA,nrow(lookup.table)*(nr_cols)),nrow = nrow(lookup.table),ncol = nr_cols)
  
  for(ci in 1:nr_cols){ 
    if(MESSAGES){
      print(paste("Sorting statistics, m-wise:",as.character(ci+minm-1)))
    }
    #computing the ecdf over the unique values
    m.stat.negative.ecdf[[ci]] = ecdf(-1*c(lookup.table[,ci],Inf))
    m.stat.size[[ci]] = length(lookup.table[,ci] +1)
    
    #compte the null ranks
    p.val[,ci]=m.stat.negative.ecdf[[ci]](-1* lookup.table[,ci])
    
  }
  #Compute the actual statistics
  for(mi in 1:nr.mmax){
      mim=maxm[mi] # the current k
      if(MESSAGES){      print(paste("Computing Statistic for M.max:",mim))    }
      
      if(MESSAGES){     print("Finding fields with NA")   }
      ind.na=is.na(lookup.table[,(mim-minm+1)]) #find N.A.'s
      if(MESSAGES){     print("Calling quick compute function")   }
      mim.param = mim
      if(variant %in% c('ADP-ML','ADP-EQP-ML')){
        mim.param = ncol(p.val)+1
      }
      qc=.hhg.univariate.object.fast.compute.statistic(p.val,mmax=mim.param,mmin=minm,ind.na,MESSAGES)
      MinP[,mi]=qc$MinP
      m_chosen[,mi]=qc$m_chosen
      Fisher[,mi]=qc$Fisher
      MinP_ECDF_LIST[[mi]]=ecdf(c(qc$MinP,0))
      FISHER_NEGATIVE_ECDF_LIST[[mi]]=ecdf(c(-qc$Fisher,-Inf))
  }    
  
  #create the object
  ret$m.stat.negative.ecdf = m.stat.negative.ecdf
  ret$m.stat.size = m.stat.size
  
  
  ret$MinP=MinP
  ret$m_chosen=m_chosen
  if(variant %in% c('ADP-ML','ADP-EQP-ML')){
    temp = .ml_index_to_chosen_vector(m_chosen[,1]-1,maxm,minm)#-1 is for going back from m to index
    ret$m_chosen = temp$m
    ret$l_chosen = temp$l
  }
  ret$Fisher=Fisher
  ret$MinP_ECDF_LIST=MinP_ECDF_LIST
  ret$FISHER_NEGATIVE_ECDF_LIST=FISHER_NEGATIVE_ECDF_LIST
    
  ret$type=type
  ret$variant=variant
  ret$size=size
  ret$score.type=score.type
  ret$aggregation.type=aggregation.type
  ret$minm = minm 
  ret$maxm = maxm
  ret$nr.atoms=nr.atoms
  ret$additional=additional
  ret$compress = compress
  ret$compress.p0 = compress.p0
  ret$compress.p = compress.p
  ret$compress.p1 = compress.p1
  if(ret$compress){ #compression is needed
    ret$m.stat.lean.null = list()
    for(ci in 1:ncol(lookup.table)){ 
      if(MESSAGES){
        print(paste("Compressing statistics, column-wise:",as.character(ci+minm-1)))
      }
      ret$m.stat.lean.null[[ci]] = lean_null_constructor(lookup.table[,ci],compress.p0,compress.p,compress.p1)
    }
    ret$m.stat.negative.ecdf=NULL
  }
  class(ret)='UnivariateObject'
  return(ret)
}

#finding pval in null table for specifiv m. this is the overloaded method for compressed and not compressed null tables (lean null tables)
UnivariateObject.findpval=function(UO,value,mindex){
  compress_flag = FALSE
  if('compress' %in% names(UO)){
    if(UO$compress){
      compress_flag=TRUE
    }
    else{
      compress_flag=FALSE
    }
  }else{
    compress_flag=FALSE
  }
  if(compress_flag){
    return(lean_null_p.val(value,UO$m.stat.lean.null[[mindex]]))
  }else{
    return(UO$m.stat.negative.ecdf[[mindex]](-1 * value))
  }
}

#function for printing the univariate object
print.UnivariateObject = function(x,...){
  cat(paste0('Univariate Null Table Object \n'))
  cat(paste0('Type: ',x$type,'\n'))
  cat(paste0('Variant: ',x$variant,'\n'))
  if(x$type == 'Independence'){
    cat(paste0('Sample Size: ',x$size,'\n'))  
  }else if(x$type == 'KSample'){
    cat(paste0('Group Sizes: ','\n'))  
    cat(paste0(x$size))
    cat('\n')
  }
  cat(paste0('Score Type: ',x$score.type,'\n'))
  cat(paste0('Aggregation Type: ',x$aggregation.type,'\n'))
  cat(paste0('mmin: ',x$minm,'\n'))
  cat(paste0('mmax: ',x$maxm,'\n'))
  if(is.element(x$variant,c('ADP-EQP','ADP-EQP-ML','KSample-Equipartition'))){
    cat(paste0('Equipartition nr.atoms: ',x$nr.atoms,'\n'))
  }
  cat(paste0('additional: ','\n'))
  cat((x$additional))
}

#function for printing the univariate statistic object.
print.UnivariateStatistic = function(x,...){
  if(x$stat.type == 'Independence-Stat'){
    .print.ind.stat(x)
  }else if(x$stat.type == 'KSample-Stat'){
    .print.ks.stat(x)
  }else if(x$stat.type == 'Independence-Combined'){
    .print.ind.combined(x)
  }else if(x$stat.type == 'KSample-Combined'){
    .print.ks.combined(x)
  }

}

.print.ml.pretty=function(vector,mmin,mmax){
  temp_mat = matrix(vector,nrow = (mmax-mmin+1),byrow = F)
  colnames(temp_mat) = paste0('m=',(mmin):(mmax))
  rownames(temp_mat) = paste0('l=',(mmin):(mmax))
  print(temp_mat,digits = 3)
}

#auxliary function used for printing hhg.univariate.ind.stat
.print.ind.stat=function(x){
  cat(paste0('HHG univariate independence statistic of type: \n'))
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  
  cat(paste0(agg.text,' of ',x$variant, ' on ',score.text, ' scores.\n\n'))  
  if(x$aggregation.type == 'sum'){
    cat('Single m statistics are the sum of scores over All Derived Partitions (ADP) of the data.\nStatistics are normalized by the number of possible partitions and sample size.\n') 
  }else if(x$aggregation.type == 'max'){
    cat('Single m statistics are the maximum score acheived by a partition of size m.\n') 
  }
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size: ',x$size,' \n\n'))
  if(is.element(x$variant,c('ADP-EQP','ADP-EQP-ML'))){
    cat(paste0('Equipartition nr.atoms: ',x$nr.atoms,'\n'))
  }
  if('statistic' %in% names(x)){
    cat(paste0('Statistics, by partition size: \n'))
    if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
      .print.ml.pretty(x$statistic,x$mmin,x$mmax)
    }else{
      print(x$statistic,digits = 3)  
    }
    
  }
  if('sum.lr' %in% names(x)){
    cat(paste0('Sum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
      .print.ml.pretty(x$sum.lr,x$mmin,x$mmax)
    }else{
      print(round(x$sum.lr,digits = 3))
    }
    cat('\n')
    
  }
  if('max.lr' %in% names(x)){
    cat(paste0('Maximum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
      .print.ml.pretty(x$max.lr,x$mmin,x$mmax)
    }else{
      print(round(x$max.lr,digits = 3))
    }
    cat('\n')
    
  }
  if('sum.chisq' %in% names(x)){
    cat(paste0('Sum over Pearson Chi Square scores of Partitions, by partition size: \n'))
    if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
      .print.ml.pretty(x$sum.chisq,x$mmin,x$mmax)
    }else{
      print(round(x$sum.chisq,digits = 3))
    }
    cat('\n')
  }
  if('max.chisq' %in% names(x)){
    cat(paste0('Maximum over Pearson Chi Square scores of Partitions, by partition size:  \n'))
    if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
      .print.ml.pretty(x$max.chisq,x$mmin,x$mmax)
    }else{
      print(round(x$max.chisq,digits = 3))
    }
    cat('\n')
    
  }
}

#auxliary function used for printing hhg.univariate.ks.stat
.print.ks.stat = function(x){
  cat(paste0('HHG univariate ksample statistic of type: \n'))
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',score.text, ' scores over possible partitions.\n\n'))  
  if(x$aggregation.type == 'sum'){
    cat('Single m statistics are the sum of scores over all possible partitions of the data.\nStatistics are normalized by the number of possible partitions and sample size.\n') 
  }else if(x$aggregation.type == 'max'){
    cat('Single m statistics are the maximum score acheived by a partition of size m.\n') 
  }
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size, by groups: ',' \n'))
  cat(paste0(x$size))
  cat('\n\n')
  if(is.element(x$variant,c('KSample-Equipartition'))){
    cat(paste0('Equipartition nr.atoms: ',x$nr.atoms,'\n'))
  }
  if('statistic' %in% names(x)){
    cat(paste0('Statistics, by partition size: \n'))
    print(x$statistic)
  }
  if('sum.lr' %in% names(x)){
    cat(paste0('Sum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    print(round(x$sum.lr,digits = 3))
    cat('\n')
    
  }
  if('max.lr' %in% names(x)){
    cat(paste0('Maximum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    print(round(x$max.lr,digits = 3))
    cat('\n')
    
  }
  if('sum.chisq' %in% names(x)){
    cat(paste0('Sum over Pearson Chi Square scores of Partitions, by partition size: \n'))
    print(round(x$sum.chisq,digits = 3))
    cat('\n')
    
  }
  if('max.chisq' %in% names(x)){
    cat(paste0('Maximum over Pearson Chi Square scores of Partitions, by partition size:  \n'))
    print(round(x$max.chisq,digits = 3))
    cat('\n')
    
  }
}

#auxliary function used for printing hhg.univariate.ind.combined.test
.print.ind.combined = function(x){
  cat(paste0('HHG univariate combined independence statistic\n'))
  cat('Statistics type combined: \n')
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',x$variant, ' on ',score.text, ' scores.\n\n'))
  if(x$aggregation.type == 'sum'){
    cat('Single m statistics are the sum of scores over All Derived Partitions (ADP) of the data.\nStatistics are normalized by the number of possible partitions and sample size.\n') 
  }else if(x$aggregation.type == 'max'){
    cat('Single m statistics are the maximum score acheived by a partition of size m.\n') 
  }
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size: ',x$size,' \n\n'))
  if(is.element(x$variant,c('ADP-EQP','ADP-EQP-ML'))){
    cat(paste0('Equipartition nr.atoms: ',x$nr.atoms,'\n'))
  }
  cat('Single m (partition size) statistics:\n')
  if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
    .print.ml.pretty(x$m.stats,x$mmin,x$mmax)
  }else{
    print(x$m.stats,digits = 3)  
  }
  cat('\n')
  
  cat('Single m (partition size) pvalues:\n')
  if(x$variant %in% c('ADP-ML','ADP-EQP-ML')){
    .print.ml.pretty(x$pvalues.of.single.m,x$mmin,x$mmax)
  }else{
    print(x$pvalues.of.single.m,digits = 3)  
  }
  cat('\n\n')
  
  if('MinP' %in% names(x)){
    cat('MinP Statistic - Test statistic is minimum of above single m pvalues: ')
    cat(round(x$MinP,digits = 4))
    cat('\n\n')
    
    cat('Partition size with minimum p-value:\n')
    cat(round(x$MinP.m.chosen,digits = 1))
    cat('X')
    if('MinP.l.chosen' %in% names(x)){
      cat(round(x$MinP.l.chosen,digits = 1))  
    }else{
      cat(round(x$MinP.m.chosen,digits = 1))  
    }
    
    cat('\n\n')
    
    cat('p-value for MinP test:')
    cat(round(x$MinP.pvalue,digits = 4))
    cat('\n\n')
  }
  
  if('Fisher' %in% names(x)){
    cat('Fisher Statistic - Test statistic is -sum(log(pval)), over  single m pvalues: ')
    cat(round(x$Fisher,digits = 4))
    cat('\n\n')
    
    cat('p-value for Fisher test:')
    cat(round(x$Fisher.pvalue,digits = 4))
    cat('\n\n')
  }
}

#auxliary function used for printing hhg.univariate.ind.ks.test
.print.ks.combined = function(x){
  cat(paste0('HHG univariate combined K-sample statistic\n'))
  cat('Statistics type combined: \n')
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',score.text, ' scores over possible partitions.\n\n'))
  if(x$aggregation.type == 'sum'){
    cat('Single m statistics are the sum of scores over all possible partitions of the data.\nStatistics are normalized by the number of possible partitions and sample size.\n') 
  }else if(x$aggregation.type == 'max'){
    cat('Single m statistics are the maximum score acheived by a partition of size m.\n') 
  }
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size, by groups: ',' \n'))
  cat(paste0(x$size))
  cat('\n\n')
  if(is.element(x$variant,c('KSample-Equipartition'))){
    cat(paste0('Equipartition nr.atoms: ',x$nr.atoms,'\n'))
  }
  
  cat('Single m (partition size) statistics:\n')
  print(round(x$m.stats,digits = 4))
  cat('\n')
  
  cat('Single m (partition size) pvalues:\n')
  print(round(x$pvalues.of.single.m,digits = 4))
  cat('\n\n')
  
  if('MinP' %in% names(x)){
    cat('MinP Statistic - Test statistic is minimum of above single m pvalues: ')
    cat(round(x$MinP,digits = 4))
    cat('\n\n')
    
    cat('Partition size with minimum p-value(# of cells):\n')
    cat(round(x$MinP.m.chosen,digits = 4))
    cat('\n\n')
    
    cat('p-value for MinP test:')
    cat(round(x$MinP.pvalue,digits = 4))
    cat('\n\n')
  }
  
  if('Fisher' %in% names(x)){
    cat('Fisher Statistic - Test statistic is -sum(log(pval)), over  single m pvalues: ')
    cat(round(x$Fisher,digits = 4))
    cat('\n\n')
    
    cat('p-value for Fisher test:')
    cat(round(x$Fisher.pvalue,digits = 4))
    cat('\n\n')
  }
}

#Internal Function for fast computing MinP & Fisher statistics by p.values (works on a matrix of values, not just a single observations)
.hhg.univariate.object.fast.compute.statistic=function(p.val,mmax,mmin=2,ind.na,MESSAGES=F){
  nr.row=nrow(p.val)
  if(MESSAGES){ print(paste("Fast Compute, # of rows is: ",nr.row,sep=''))}
  m.bound=max(mmax)
  nr.kmax=length(mmax)
  
  m_chosen=rep(NA,nr.row)
  Statistic=rep(NA,nr.row)
  Fisher=rep(NA,nr.row)
  
  if(MESSAGES){      print("Generating : finding indexes of mins")    }
  if(mmax-mmin>0){
    ind=apply(p.val[,1:(mmax-mmin+1)],1,which.min) # before na check
  }else{
    ind=rep(2,nrow(p.val))
  }
  if(MESSAGES){     print("Generating : writing results")   }
  for(ri in 1:nr.row){
    
    if(mmax>2){
      Statistic[ri]=p.val[ri,ind[ri]]
      Fisher[ri]=(-1)*sum(log(p.val[ri,(1 ):(mmax-mmin+1)]))
    }else{
      Statistic[ri]=p.val[ri]
      Fisher[ri]=(-1)*sum(log(p.val[ri]))
    }
  }
  m_chosen=ind+1
  Statistic[ind.na]=NA #handle na's
  m_chosen[ind.na]=NA
  Fisher[ind.na]=NA
  
  ret=list()
  ret$MinP=Statistic
  ret$Fisher=Fisher
  ret$m_chosen=m_chosen
  return(ret)
}

#function for constructing a null table object for independece statistics
hhg.univariate.ind.nulltable=function(size,mmin=2,mmax = max(floor(sqrt(size)/2),2),variant = 'ADP',aggregation.type = 'sum',score.type='LikelihoodRatio', w.sum = 0, w.max = 2,nr.replicates=1000,keep.simulation.data=F, nr.atoms = nr_bins_equipartition(size),compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001){  
  correct.mi.bias = F #mi estimation correction not supported yet
  type='Independence'
  
  #input checks:
  
  if(is.null(size) |any(is.null(mmax)) | is.null(nr.replicates) | is.null(variant) ){
    stop("size, mmax, nr.replicates should be integers. mmax and mmin should be less than size and >=2")
  }
  
  if(is.na(size) | is.na(nr.replicates) |is.na(variant) |any(is.na(mmax))){
    stop("size, Max.m, nr.replicates should be integers. mmax and mmin should be less than size and >=2")
  }
  
  
  for(i in 1:length(mmax)){
    if(mmax[i]!=as.integer(mmax[i]) & !is.na(mmax[i])){
      stop("Max.m must be vector of integers")
    }
  }  
  
  .hhg.univariate.check.inputs(type,variant,size,score.type,aggregation.type,nr.atoms = nr.atoms)
  
  if(size!= as.integer(size) |  nr.replicates!=as.integer(nr.replicates) | size<1 | max(mmax) > size | min(mmax) <2 | mmin<2| nr.replicates < 1){
    stop("size, mmax, nr.replicates should be integers. mmax shoudld be less than size and >=2")
  }
  
  if(score.type=='both' | aggregation.type=='both'){
    stop('null table can only be constructed of a single statistic, \'both\' not allowed in score.type or aggregation.type')
  }
  .check_compress_args(compress.p,compress.p0,compress.p1)
  #create null table statistics
  mat_size = max(mmax)-mmin+1
  if(is.element(variant,c('ADP-ML','ADP-EQP-ML'))){
    mat_size = mat_size*mat_size
  }
  m.stats=matrix(rep(NA,(mat_size)*nr.replicates),nrow = nr.replicates,ncol = mat_size)
  x=(1:size)
  y_samp=(1:size)
  titles = NULL
  for(r in 1:nr.replicates){
    y=sample(y_samp)
    current_test=hhg.univariate.ind.stat(x,y,variant,aggregation.type,score.type,mmax = max(mmax),mmin = mmin,w.sum = w.sum,w.max = w.max,nr.atoms = nr.atoms)
    m.stats[r,]=current_test$statistic
    if(r ==1){titles = names(current_test$statistic)}
  }
  
  #create objects
  colnames(m.stats)=titles
  univariate.object = .hhg.univariate.object.constructor.from.mstats(m.stats,mmin,mmax,type,variant,size,score.type,aggregation.type,nr.atoms,c(correct.mi.bias,w.sum,w.max),ml_variant = is.element(variant,c('ADP-ML','ADP-EQP-ML')),compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1)
  ret=list()
  if(keep.simulation.data){
    ret$m.stats = m.stats  
  }
  ret$univariate.object = univariate.object
  return(ret)
}

# Perform MinP/Fisher Test from m.stats
.hhg.univariate.compute.from.mstats=function(m.stats,univariate.object,minm,maxm,ml_variant = F){
  
  #inputs checks:
  
  if( minm!=as.integer(minm) | maxm!=as.integer(maxm)){
    stop("maxm and minm must be integers must be an integer")
  }
  if(maxm-minm+1>length(m.stats)){
    s=paste("m.stats provided is shorter than required maxm. # m.stats provided:",length(m.stats), "(maybe due to clumps) maxm:",maxm,' minm:',minm,sep='')
    stop(s)
  }
  if(max(univariate.object$maxm)<maxm){
    stop("Univariate object's maximum number of partitions computed is not sufficient")
  }
  if((univariate.object$minm)!=minm){
    stop("Univariate object's minimum number of partitions computed is not same as required in statistic computation")
  }
  
  #compute:
  ret=list()
  mat_size =maxm-minm+1
  if(ml_variant){
    mat_size = mat_size * mat_size
  }
  pm=rep(NA,mat_size)
  if(is.na(m.stats[length(m.stats)])){ #cannot perform test 
    stop("Missing values in m.stats")
  }
  for(pi in 1:(mat_size)){
    #compute the p.vals of the m.stats by the MinP.Object
    pm[pi]=UnivariateObject.findpval(univariate.object,m.stats[pi],pi) #univariate.object$m.stat.negative.ecdf[[pi]](-1*m.stats[pi])
  }
  
  #create object and return
  ret$MinP.m.chosen=which.min(pm)+1 #find the minimum until the last if, chosen m is index, not ML variant output
  ret$MinP=pm[ret$MinP.m.chosen-1]
  ret$pvalues.of.single.m=pm
  ret$Fisher=(-1)*sum(log(pm[1:(mat_size)]))
  if(sum(univariate.object$maxm==maxm)>0){
    current_m=which(univariate.object$maxm==maxm)
    ret$MinP.pvalue=univariate.object$MinP_ECDF_LIST[[current_m]](ret$MinP)
    ret$Fisher.pvalue=univariate.object$FISHER_NEGATIVE_ECDF_LIST[[current_m]](-ret$Fisher)
  }
  if(ml_variant){#this is for the output
    temp = .ml_index_to_chosen(which.min(pm),maxm,minm)
    ret$MinP.m.chosen = temp$m
    ret$MinP.l.chosen = temp$l
  }
  return(ret)
}

#function for compiling a null table object from existing m.stats
#(like the ones given on the site, or the ones used in the code example for function and vignette, for computing a null table using multiple cores)
hhg.univariate.nulltable.from.mstats=function(m.stats,minm,maxm,type,variant,size,score.type,aggregation.type, w.sum = 0, w.max = 2, keep.simulation.data=F,nr.atoms = nr_bins_equipartition(sum(size)),compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001){
  correct.mi.bias = F #not supported yet
  .check_compress_args(compress.p,compress.p0,compress.p1)
  ret=list()
  #call the private function to create and object, and return it. all inputs checks are to be dont inside
  ret$univariate.object = .hhg.univariate.object.constructor.from.mstats(m.stats,minm,maxm,type,variant,size,score.type,aggregation.type,additional =c(correct.mi.bias,w.sum,w.max), nr.atoms = nr.atoms,ml_variant = is.element(variant,c('ADP-ML','ADP-EQP-ML')),compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1)
  if(keep.simulation.data){
    ret$m.stats=m.stats
  }
  return(ret)
}

#function for performing the combined independence test, combining different partition sizes.
hhg.univariate.ind.combined.test=function(X,Y=NULL,NullTable=NULL,mmin=2,mmax=max(floor(sqrt(length(X))/2),2),variant='ADP',aggregation.type='sum',score.type='LikelihoodRatio', w.sum = 0, w.max = 2 ,combining.type='MinP',nr.perm=100,nr.atoms = nr_bins_equipartition(length(X)),compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001,keep.simulation.data=T){
  
  correct.mi.bias = F #not supported yet
  # function gets either a vector at X and Y, or the result of hhg.univariate.ind.stat at X and null Y. function
  # may also receive null table object or create null table by parameters. since exiting data and null table need to work toghether,
  # or null table and hhg.univariate.ind.stat result, the function has the following priorities:
  
  #  if no null table
  #     if statisticc is found
  #         take statistic parameters
  #     else
  #         take parameters from function inputs
  #  else (there is a null table)
  #     if there is a statistic, 
  #       take its parameters.
  #       try to find an mmax in the null table, that fits statistic.
  #       if not found, take the min mmax (it could be the statistic will still not have enough terms, or it wont use all   		terms)
  #     else
  #       take null table parameters
    #these variables house the parameters chosen by the above logic.
    current_null_table=NULL
    current_mmax=NULL
    current_mmin=NULL
    current_variant=NULL
    current_aggregation.type=NULL
    current_score.type=NULL
    current_correct.mi.bias = NULL
    current_size = NULL
    current_w.sum = NULL
    current_w.max = NULL
    current_nr.atoms = NULL
    XY_is_stat=F
    m.stats=rep(NA,mmax-mmin+1) #if Y is null and X is statistic, this doesn't matter, since it is overridden (it is of a false dimensionality)
    
    #input checks
    .check_compress_args(compress.p,compress.p0,compress.p1)
    
    if(!is.element(combining.type,c('MinP','Fisher','Both'))){
      stop('combining.type should be \'MinP\' , \'Fisher\'  or \'Both\' ')  
    }
    if(is.null(Y)){ if(is.list(X)){ if( class(X) =='UnivariateStatistic'){ if(X$stat.type=='Independence-Stat'){
      XY_is_stat=T #we already acctually have the statistics
    }}}}
    if(XY_is_stat == F & is.null(Y)){
      stop('Y not supplied or X is not valid result of hhg.univariate.ind.stat.')
    }
    if(is.null(NullTable)){ #handle the case of no null table found, generate, according to parameters. remember to take parameters from object if found.
      current_mmax = mmax
      current_mmin = mmin
      current_size = length(X)
      current_variant = variant
      current_aggregation.type = aggregation.type
      current_score.type = score.type
      current_correct.mi.bias = correct.mi.bias
      current_w.sum = w.sum
      current_w.max = w.max
      current_nr.atoms = nr.atoms
      current_size = length(X)
      if(XY_is_stat){ #handle object
        current_mmax=X$mmax
        current_mmin=X$mmin
        current_variant=X$variant
        if("nr.atoms" %in% names(X)){
          current_nr.atoms = X$nr.atoms
        }
        current_size = X$size
        current_aggregation.type = X$aggregation.type
        current_score.type = X$score.type
        current_correct.mi.bias = X$additional[1]
        current_w.sum = X$additional[2]
        current_w.max = X$additional[3]
      }else{
        if(length(X) != length(Y)){
          stop('X and Y not of the same length')
        }
      }
      .hhg.univariate.check.inputs('Independence',current_variant,current_size,current_score.type,current_aggregation.type,nr.atoms = current_nr.atoms)
      
      current_null_table = hhg.univariate.ind.nulltable(size= current_size,mmin,mmax,variant,aggregation.type,score.type,nr.replicates = nr.perm, w.sum = w.sum, w.max = w.max,keep.simulation.data = keep.simulation.data,nr.atoms = current_nr.atoms,compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1)
      
    }else{ #null table is found
      current_null_table = NullTable
      #checking its of the right type
      if(!is.list(current_null_table) | !('univariate.object' %in% names(current_null_table))){
        stop('NullTable supplied is not valid. construct null table using hhg.univariate.ind.nulltable')
      }
      if(current_null_table$univariate.object$type!='Independence'){
        stop('null table type is not \'Independence\' ')
      }
      
      current_variant = NullTable$univariate.object$variant
      current_aggregation.type = NullTable$univariate.object$aggregation.type
      current_score.type = NullTable$univariate.object$score.type
      current_mmin = NullTable$univariate.object$minm
      if("nr.atoms" %in% names(NullTable$univariate.object)){
        current_nr.atoms = NullTable$univariate.object$nr.atoms
      }# if it is not found, it is a pre 1.6 null table, and there is no problem.... value of nr.atoms doesnt matter
        
      for(m in current_null_table$univariate.object$maxm){
        if(m == mmax){current_mmax=m}
      }
      if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}

      current_correct.mi.bias = NullTable$univariate.object$additional[1]
      current_w.sum = NullTable$univariate.object$additional[2]
      current_w.max = NullTable$univariate.object$additional[3]
      current_size = NullTable$univariate.object$size
      if(XY_is_stat){ #handle the case we have both a null table and an object
        
        current_mmax=NULL
        for(m in current_null_table$univariate.object$maxm){
          if(m == mmax){current_mmax=m}
        }
        if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}
        
        
        
        if(current_mmax > X$mmax){
          stop('Null Table mmax bigger than computed stat m')
        }
        if(current_mmin != X$mmin){
          stop('Null Table mmin is different than stat m')
        }
        if(current_variant != X$variant){
          stop('Null Table and Statistic not of same variant')
        }
        if(current_size != X$size){
          stop('Null Table and Statistic not of same size')
        }
        if(current_aggregation.type != X$aggregation.type){
          stop('Null Table and Statistic not of same aggregation type')
        }
        if(current_score.type != X$score.type){
          stop('Null Table and Statistic not of same score type')
        }
        if(current_correct.mi.bias != X$additional[1]){
          stop('Null Table and Statistic not of same bias correction')
        }
        if(current_w.sum != X$additional[2]){
          stop('Null Table and Statistic not of same w.sum')
        }
        if(current_w.max != X$additional[3]){
          stop('Null Table and Statistic not of same w.max')
        }
        if(is.element(current_variant,c('ADP-EQP','ADP-EQP-ML'))){
          if(current_nr.atoms != X$nr.atoms){
            stop(paste0('Null Table and Statistic not of same nr.atoms of ADP-EQP/ADP-EQP-ML variant. Null table is: ',current_nr.atoms,' while statistic object is: ', X$nr.atoms))
          }
        }
      }else{
      }  
        
    }#end of null table found
    if(XY_is_stat){
      mat_size = current_mmax - current_mmin +1
      if(is.element(X$variant,c('ADP-ML','ADP-EQP-ML'))){
        mat_size = mat_size * mat_size
      }
      m.stats[1:(mat_size)]  = X$statistic[1:mat_size]
    }else{
      if(length(X)!=current_null_table$univariate.object$size){
        stop(paste0('Sample size is ',length(X),' while null table was constructed on sample size ',NullTable$univariate.object$size))
      }        
      m.stats=hhg.univariate.ind.stat(X,Y,current_variant,current_aggregation.type,current_score.type,mmin = current_mmin , mmax = current_mmax,w.sum = current_w.sum , w.max = current_w.max,nr.atoms = current_nr.atoms)$statistic        
    }
    #compute statistics and pvalues.
    computed = .hhg.univariate.compute.from.mstats(m.stats,current_null_table$univariate.object,current_mmin,current_mmax,ml_variant = is.element(current_variant,c('ADP-ML','ADP-EQP-ML')))
  
    #pack object
    ret=computed
    ret$m.stats=m.stats
    create_params = .ml_by_variant(current_variant,current_mmin,current_mmax)
    m_param = create_params$m_param
    l_param = create_params$l_param
    titles = create_params$titles
    names(ret$m.stats)= titles#paste0('m.',(current_mmin:current_mmax))
    names(ret$pvalues.of.single.m) = paste0('pval.',titles)
    if(is.null(NullTable)){
      ret$generated_null_table=current_null_table
    }
    if(combining.type == 'MinP'){
      ind.to.remove = which(names(ret) %in% c('Fisher','Fisher.pvalue'))
      ret=ret[-ind.to.remove]
    }else if(combining.type =='Fisher'){
      ind.to.remove = which(names(ret) %in% c('MinP','MinP.pvalue','MinP.m.chosen'))
      ret=ret[-ind.to.remove]
    }else if(combining.type =='Both'){
      #do nothing, no need to remove
    }
    ret$stat.type = 'Independence-Combined'
    ret$mmin = current_mmin
    ret$mmax = current_mmax
    ret$score.type = current_score.type
    ret$aggregation.type = current_aggregation.type
    ret$variant = current_variant
    ret$w.sum = current_w.sum
    ret$w.max = current_w.max
    ret$size = current_size
    ret$nr.atoms = current_nr.atoms
    class(ret)='UnivariateStatistic'
    return(ret)
}

#function for computing the pvalue of an independence statistic
hhg.univariate.ind.pvalue=function(statistic, NullTable, m=min(statistic$mmax,4),l=m){
  #input checks
  if(!is.list(statistic) | class(statistic)!='UnivariateStatistic'){
    stop(' statistic is not result of hhg.univariate.ind.stat')
  }
  if(statistic$stat.type != 'Independence-Stat'){
    stop(' statistic is not result of hhg.univariate.ind.stat')
  }
  if(!is.list(NullTable) | !('univariate.object' %in% names(NullTable))){
    stop('NullTable supplied is not valid. construct null table using hhg.univariate.ind.nulltable')
  }
  
  uvo=NullTable$univariate.object
  if(class(uvo)!= 'UnivariateObject'){
    stop('Null table supplied does not contain valid univariate object')
  }
  if(length(statistic$size) != length(uvo$size)){
    stop(paste0('statistic group sizes are ',statistic$size,' while null table group sizes are:',uvo$size ))
  }
  if(any(statistic$size != uvo$size)){
    stop(paste0('statistic size is ',statistic$size,' while null table size is ', uvo$size))
  }
  if(statistic$variant != uvo$variant){
    stop(paste0('statistic variant is ',statistic$variant,' while null table variant is ', uvo$variant))
  }
  if(statistic$aggregation.type != uvo$aggregation.type){
    stop(paste0('statistic aggregation.type is ',statistic$aggregation.type,' while null table aggregation.type is ', uvo$aggregation.type))
  }
  if(statistic$score.type != uvo$score.type){
    stop(paste0('statistic score.type is ',statistic$score.type,' while null table score.type is ', uvo$score.type))
  }
  if(statistic$type != uvo$type){
    stop(paste0('statistic type is ',statistic$type,' while null table type is ', uvo$type))
  }
  
  
  if(statistic$additional[1] != uvo$additional[1]){
    stop(paste0('MI bias correction in test sample is : ',statistic$additional[1] , ' while in null table is : ',uvo$additional[1]))
  }
  if(TRUE){#if(statistic$additional[1]){
    if(statistic$aggregation.type=='max' & statistic$additional[3] != uvo$additional[3]){
      stop(paste0('w.max in test sample: ', statistic$additional[3] , ' while in null table: ', uvo$additional[3]))
    }
    if(statistic$aggregation.type=='sum' & statistic$additional[2] != uvo$additional[2]){
      stop(paste0('w.sum in test sample: ', statistic$additional[2] , ' while in null table: ', uvo$additional[2]))
    }
  }
  
  if(is.element(statistic$variant,c('ADP-EQP','ADP-EQP-ML'))){
    if(statistic$nr.atoms != uvo$nr.atoms){
      stop(paste0('Null Table and Statistic not of same nr.atoms of ADP-EQP/ADP-EQP-ML variant. Null table is: ',uvo$nr.atoms,' while statistic object is: ', statistic$nr.atoms))
    }
  }
  if(is.element(statistic$variant,c('ADP-ML','ADP-EQP-ML'))){
    flag = F
   if(is.null(l)){ flag = T }else if(!is.numeric(l)){flag = T}else if(l != round(l)){flag = T}
   if(flag ==F){if(l<2 ){flag=T}}
   if(flag == T){
      stop("for ADP-ML and ADP-EQP-ML variants, l should also be supplied, being an integer between m.min and m.max.")   
   }
  }
   
  
  m.ind = NULL
  m.stat.ind = NULL
  ret=list()
  if(statistic$stat.type == 'Independence-Stat'){
    if(uvo$minm >m | max(uvo$maxm) < m){
      stop(paste0('null table does not include number of partitions m: ',m))
    }
    
    if(statistic$mmin >m | max(statistic$mmax) < m){
      stop(paste0('statistic  does not include number of partitions m: ',m))
    }
    if(is.element(statistic$variant,c('ADP-ML','ADP-EQP-ML'))){
      if(uvo$minm >l | max(uvo$maxm) < l){
        stop(paste0('null table does not include number of partitions m: ',m))
      }
      
      if(statistic$mmin >l | max(statistic$mmax) < l){
        stop(paste0('statistic  does not include number of partitions l: ',l))
      }
    }
    #pvalue can be computed
    m.ind = m - uvo$minm +1
    m.stat.ind = m - statistic$mmin +1
    if(is.element(statistic$variant,c('ADP-ML','ADP-EQP-ML'))){
      m.ind = (uvo$maxm - uvo$minm +1)*(m - uvo$minm ) + (l - uvo$minm +1)
      m.stat.ind = (statistic$mmax - statistic$mmin +1)*(m - statistic$mmin ) + (l - statistic$mmin +1)
    }
    stat=statistic$statistic[m.stat.ind] #the single m statistic value under the statistic object
    return(UnivariateObject.findpval(uvo,stat,m.ind))#return(uvo$m.stat.negative.ecdf[[m.ind]](-1 * stat))

  }else{
    stop('incorrect statistic type: statistic should be the result of hhg.univariate.ind.stat')
  }  
}

#function for computing K sample statistic
hhg.univariate.ks.stat=function(x, y,variant = 'KSample-Variant',aggregation.type='sum',score.type='LikelihoodRatio', mmax = max(4,round(min(table(y))/3)),mmin=2,nr.atoms= nr_bins_equipartition(length(x))){
  type='KSample'
  #input checks
  if(is.null(x) |is.null(y)){
    stop("x & y should be vectors of doubles.")
  }
  if(length(x)!= length(y)){
    stop("X & Y not the same length")
  }
  if(!is.numeric(mmin) | !is.numeric(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin !=round(mmin) | mmax !=round(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin<2 | mmax<mmin){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(is.null(nr.atoms) | is.na(nr.atoms)){
    stop('For KSample Equipartition nr.atoms should be between 2 and sample size, and bigger than mmax')
  }
  if(variant =='KSample-Equipartition' & (nr.atoms > length(x) | nr.atoms <2 | mmax>nr.atoms ) ){
    stop('For KSample Equipartition nr.atoms should be between 2 and sample size, and bigger than mmax')
  }
  if(nr.atoms!=round(nr.atoms)){
    stop('nr.atoms should be an integer')
  }
  
  
  .hhg.univariate.check.inputs(type,variant,(table(y)),score.type,aggregation.type,nr.atoms)
  
  #check which types should be computed
  flag_SL = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_SC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_ML = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  flag_MC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  results_SL=NULL
  results_SC=NULL
  results_ML=NULL
  results_MC=NULL
  param_dynslice_type = 1
  param_sm_type = 0
  if(variant == 'KSample-Equipartition'){
    param_dynslice_type = 2
    param_sm_type = 1
  }
  if(mmax==mmin & variant == 'KSample-Variant'){
    if(flag_ML| flag_MC){
      test=.hhg.univariate.ks.stat.max(x, y, variant = 'Mm',Max.m =mmax ,m.stats.wanted=T)
      if(length(test$m.stats)< 2*(mmax-1)){ #this is due to use of clumps in the original max computation algorithem
        stop(paste0('Too many clumps in data, cannot compute ', mmax ,' partitions'))
      }
      results_ML=test$m.stats[(mmin-1)] 
      results_MC = test$m.stats[(mmax-1)+(mmin-1)]
    }
    if(flag_SC | flag_SL){
      test=.xdp.test.k.sample(x,y,mmax)
      results_SL = test$sum.lr
      results_SC = test$sum.chisq
    }
  }else{
    if(flag_ML | flag_MC){
      test=.hhg.univariate.ks.stat.max(x, y, variant = 'Mm',Max.m =mmax ,m.stats.wanted=T,DS_type = param_dynslice_type,nr_equipartition_bins = nr.atoms) #computing the Mm statistic.
      results_ML=test$m.stats[(mmin-1):(mmax-1)] 
      names(results_ML)=paste('M.m_',(mmin:(mmax)),sep = '')  
      results_ML=results_ML[!is.na(results_ML)]
      results_MC=test$m.stats[(mmax-1) + (mmin-1):(mmax-1)] 
      names(results_MC)=paste('M.m_',(mmin:(mmax)),sep = '')  
      results_MC=results_MC[!is.na(results_MC)]
    }
    if(flag_SC | flag_SL){
      test=.xdp.test.k.sample.mk(x,y,mmax,equipartition.type = param_sm_type,equipartition.nr.bins = nr.atoms) #computing the Sm statistic
      results_SC = test[(mmin-1):(mmax-1),1]
      names(results_SC)=paste('S.m_',(mmin:mmax),sep = '')  
      results_SL = test[(mmin-1):(mmax-1),2]
      names(results_SL)=paste('S.m_',(mmin:mmax),sep = '')  
    }
  }
  #writing to output object
  ret=list()
  if(flag_SL + flag_SC + flag_ML + flag_MC==1){
    if(flag_SL){
      ret$statistic = results_SL
    }else if(flag_SC){
      ret$statistic = results_SC
    }else if(flag_ML){
      ret$statistic = results_ML
    }else if(flag_MC){
      ret$statistic = results_MC
    }
  }else{
    if(flag_SL){
      ret$sum.lr = results_SL
    }
    if(flag_SC){
      ret$sum.chisq = results_SC
    }
    if(flag_ML){
      ret$max.lr = results_ML
    }
    if(flag_MC){
      ret$max.chisq = results_MC
    }
  }
  ret$type = type
  ret$variant=variant
  ret$stat.type = 'KSample-Stat'
  ret$size = table(y)
  ret$score.type = score.type
  ret$aggregation.type = aggregation.type
  ret$mmin = mmin
  ret$mmax = mmax
  ret$nr.atoms=nr.atoms
  class(ret)='UnivariateStatistic'
  return (ret)
}

#function for creating K sample null table object
hhg.univariate.ks.nulltable=function(group.sizes,mmin=2,mmax=max(4,round(min(group.sizes)/3)),variant = 'KSample-Variant',aggregation.type='sum',score.type='LikelihoodRatio',nr.replicates=1000,keep.simulation.data=F,nr.atoms = nr_bins_equipartition(sum(group.sizes)),compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001){
  type='KSample'
  #check inputs
  if(any(is.null(group.sizes)) |any(is.null(mmax)) | is.null(nr.replicates) | is.null(variant) ){
    stop("Size should be a vector of valid group sizes (more than two groups). mmax,mmin , nr.replicates should be integers. mmax and mmin should be less than sum(group.sizes) and >=2")
  }
  
  if(any(is.na(group.sizes)) | is.na(nr.replicates) |is.na(variant) |any(is.na(mmax))){
    stop("Size should be a vector of valid group sizes (more than two groups). mmax,mmin, nr.replicates should be integers. mmax and mmin should be less than sum(group.sizes) and >=2")
  }
  
  if(any(mmax!=as.integer(mmax)) & any(is.na(mmax))){
    stop("mmax must be vector of integers")
  }
  
  .hhg.univariate.check.inputs(type,variant,group.sizes,score.type,aggregation.type,nr.atoms)
  
  if(any(group.sizes!= as.integer(group.sizes)) |  nr.replicates!=as.integer(nr.replicates) | any(group.sizes<1) | max(mmax) > sum(group.sizes) | mmin<2 |min(mmax) <2 | nr.replicates < 1){
    stop("Size should be a vector of valid group sizes (more than two groups). mmax,mmin nr.replicates should be integers. mmax and mmin should be less than sum(group.sizes) and >=2")
  }
  
  if(score.type=='both' | aggregation.type=='both'){
    stop('null table can only be constructed of a single statistic, \'both\' not allowed in score.type or aggregation.type')
  }
  
  .check_compress_args(compress.p,compress.p0,compress.p1)
  
  #perform permutations on ranks
  m.stats=matrix(rep(NA,(max(mmax)-mmin+1)*nr.replicates),nrow = nr.replicates,ncol = (max(mmax)-mmin+1))
  x=(1:sum(group.sizes))
  y_samp=NULL
  
  for(i in 1:length(group.sizes)){
    y_samp=c(y_samp,rep(i-1,group.sizes[i]))
  }
  for(r in 1:nr.replicates){
    y=sample(y_samp)
    current_test=hhg.univariate.ks.stat(x,y,variant,aggregation.type,score.type,mmax = max(mmax),mmin=mmin,nr.atoms = nr.atoms)
    m.stats[r,]=current_test$statistic
  }
  
  #compile statistics matrix to object
  colnames(m.stats)=paste0('m.',(mmin:max(mmax)))
  univariate.object = .hhg.univariate.object.constructor.from.mstats(m.stats,mmin,mmax,type,variant,group.sizes,score.type,aggregation.type,additional = NULL,nr.atoms = nr.atoms,compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1)
  ret=list()
  if(keep.simulation.data){
    ret$m.stats = m.stats  
  }
  ret$univariate.object = univariate.object
  return(ret)
}

#function for performing the combined K sample test, combining different partition sizes.
hhg.univariate.ks.combined.test=function(X,Y=NULL,NullTable=NULL,mmin=2,mmax=ifelse(is.null(Y),4,max(4,round(min(table(Y))/3))), aggregation.type='sum',score.type='LikelihoodRatio' ,combining.type='MinP',nr.perm=1000,variant='KSample-Variant', nr.atoms = nr_bins_equipartition(length(X)),compress=F,compress.p0=0.001,compress.p=0.99,compress.p1=0.000001,keep.simulation.data=T){
  # function gets either a vector at X and Y, or the result of hhg.univariate.ks.stat at X and null Y. function
  # may also receive null table object or create null table by parameters. since exiting data and null table need to work toghether,
  # or null table and hhg.univariate.ind.stat result, the function has the following priorities:
  
#  if no null table
#     if statistic is found
#         take statistic parameters
#     else
#         take parameters from function inputs
#  else (there is a null table)
#     if there is a statistic, 
#       take its parameters.
#       try to find an mmax in the null table, that fits statistic.
#       if not found, take the min mmax (it could be the statistic will still not have enough terms, or it wont use all 			terms)
#     else
#       take null table parameters
  
  #these hold the final parameters for the test, according to the above decision logic.
  current_null_table=NULL
  current_mmax=NULL
  current_mmin=NULL
  current_variant = NULL
  current_aggregation.type=NULL
  current_score.type=NULL
  current_size = NULL
  current_nr.atoms = NULL
  XY_is_stat=F # check if X is an actual UnivariateStatistic object of the required type.
  
  if(!is.numeric(mmax) | is.infinite(mmax)){ #this is to handle a case were one uses Y as null, and X as statistic
    mmax=mmin
  }
  
  m.stats=rep(NA,mmax-mmin+1)
  if(!is.element(combining.type,c('MinP','Fisher','Both'))){
    stop('combining.type should be \'MinP\' , \'Fisher\'  or \'Both\' ')  
  }
  if(is.null(Y)){ if(is.list(X)){ if( class(X) == 'UnivariateStatistic'){ if(X$stat.type=='KSample-Stat'){
    XY_is_stat=T #we already acctually have the statistics
  }}}}
  if(XY_is_stat == F & is.null(Y)){
    stop('Y not supplied or X is not valid result of hhg.univariate.ks.stat.')
  }
  
  .check_compress_args(compress.p,compress.p0,compress.p1)
  
  if((is.null(NullTable))){ #no null table is given, we have to generate by user specific parameters.
    current_mmax=mmax
    current_mmin=mmin
    current_aggregation.type = aggregation.type
    current_score.type = score.type
    current_nr.atoms = nr.atoms
    current_variant = variant
    if(XY_is_stat){ #handle the case no null table is given, but we have an object in X
      current_mmax = X$mmax
      current_mmin = X$mmin
      current_variant = X$variant
      current_size = X$size
      current_aggregation.type = X$aggregation.type
      current_score.type = X$score.type
      current_nr.atoms = X$nr.atoms
    }else{
      current_size = table(Y)
      if(length(X) != length(Y)){
          stop('X and Y not of the same length')
      }
    }
    .hhg.univariate.check.inputs('KSample',current_variant,current_size,current_score.type,current_aggregation.type,nr.atoms = current_nr.atoms)
    #generate null table
    current_null_table = hhg.univariate.ks.nulltable(current_size,current_mmin,current_mmax,current_variant,current_aggregation.type,current_score.type,nr.replicates = nr.perm,nr.atoms = current_nr.atoms,compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1,keep.simulation.data = keep.simulation.data)
  }else{ #null table is given, we have to check if its of the correct type
    current_null_table = NullTable
    if(!is.list(current_null_table) | !('univariate.object' %in% names(current_null_table))){
      stop('NullTable supplied is not valid. construct null table using hhg.univariate.ks.nulltable')
    }
    if(current_null_table$univariate.object$type != 'KSample'){
      stop('null table type is not \'KSample\' ')
    }
    current_variant = NullTable$univariate.object$variant
    current_aggregation.type = NullTable$univariate.object$aggregation.type
    current_score.type = NullTable$univariate.object$score.type
    current_mmin = NullTable$univariate.object$minm
    current_nr.atoms = NullTable$univariate.object$nr.atoms
    for(m in current_null_table$univariate.object$maxm){
      if(m == mmax){current_mmax=m
      }
    }
    if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}
    current_size = NullTable$univariate.object$size

    if(XY_is_stat){ #handle a case where we have both a null table and object
      current_mmax = NULL
      for(m in current_null_table$univariate.object$maxm){
        if(m == X$mmax){current_mmax=m
        }
      }
      if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}
      if(current_mmax > X$mmax){
        stop('Null Table mmax bigger than computed stat m')
      }
      if(current_mmin != X$mmin){
        stop('Null Table mmin is different than stat m')
      }
      if(current_variant != X$variant){
        stop('Null Table and Statistic not of same variant')
      }
      if(length(current_size) != length(X$size)){
        stop('Null Table and Statistic not of same number of groups')
      }
      if(any(sort(current_size) != sort(X$size))){
        stop('Null Table and Statistic not of same group sizes')
      }
      if(current_aggregation.type != X$aggregation.type){
        stop('Null Table and Statistic not of same aggregation type')
      }
      if(current_score.type != X$score.type){
        stop('Null Table and Statistic not of same score type')
      }
      if(is.element(current_variant,c('KSample-Equipartition'))){
        if(current_nr.atoms != X$nr.atoms){
          stop(paste0('Null Table and Statistic not of same nr.atoms of KSample-Equipartition variant. Null table is: ',current_nr.atoms,' while statistic object is: ', X$nr.atoms))
        }
      }
    }else{
    }  
    
  }#end of null table found
  if(XY_is_stat){
    m.stats[1:(current_mmax - current_mmin +1)]  = X$statistic[1:(current_mmax - current_mmin +1)]
  }else{
    tab_y=table(Y)
    if(length(current_size)!=length(tab_y)){
      stop('Test sample group sizes and null table group sizes do not match - Not same number of groups')
    }
    if(any(sort(current_size)!=sort(tab_y))){
      stop('Test sample group sizes and null table group sizes do not match!')
    }
    m.stats = hhg.univariate.ks.stat(X,Y,current_variant,current_aggregation.type,current_score.type,current_mmax,current_mmin,nr.atoms=current_nr.atoms)$statistic
    
    actual_mmax = current_mmin + length(m.stats)-1 #might we shorter than wanted for dynamic slicing
    if(actual_mmax != current_mmax){# Mm had too few clumps
      ret=list()
      ret$m.stats=m.stats
      return(ret)
      warning('Too few clumps in data, lower mmax or use \'sum\' in aggregation.type')
    }
  }
  
  # compute pvalues and MinP,Fisher statistics
  computed = .hhg.univariate.compute.from.mstats(m.stats,current_null_table$univariate.object,current_mmin,current_mmax)
  
  #pack the results
  ret=computed
  ret$m.stats=m.stats
  if(aggregation.type=='sum'){
    names(ret$m.stats)=paste0('Sm.',current_mmin:current_mmax)
  }else if(aggregation.type=='max'){
    names(ret$m.stats)=paste0('Mm.',current_mmin:current_mmax)
  }
  if(is.null(NullTable)){
    ret$generated_null_table=current_null_table
  }
  if(combining.type == 'MinP'){
    ind.to.remove = which(names(ret) %in% c('Fisher','Fisher.pvalue'))
    ret=ret[-ind.to.remove]
  }
  else if(combining.type =='Fisher'){
    ind.to.remove = which(names(ret) %in% c('MinP','MinP.pvalue','MinP.m.chosen'))
    ret=ret[-ind.to.remove]
  }else if(combining.type =='Both'){
    #do nothing, no need to remove
  }
  ret$stat.type = 'KSample-Combined'
  ret$mmin = current_mmin
  ret$mmax = current_mmax
  ret$score.type = current_score.type
  ret$aggregation.type = current_aggregation.type
  ret$variant = current_variant
  ret$size = current_size
  ret$nr.atoms = current_nr.atoms
  class(ret)='UnivariateStatistic'
  return(ret)
}

#function for computing p-values of the K sample statistic
hhg.univariate.ks.pvalue=function(statistic, NullTable,m = min(statistic$mmax,4)){
  #input checks
  if(!is.list(statistic) | class(statistic)!='UnivariateStatistic'){
    stop(' statistic is not result of hhg.univariate.ks.stat')
  }
  if(statistic$stat.type != 'KSample-Stat'){
    stop(' statistic is not result of hhg.univariate.ks.stat')
  }
  if(!is.list(NullTable) | !('univariate.object' %in% names(NullTable))){
    stop('NullTable supplied is not valid. construct null table using hhg.univariate.ks.nulltable')
  }
  uvo=NullTable$univariate.object
  if(class(uvo)!= 'UnivariateObject'){
    stop('Null table supplied does not contain valid univariate object')
  }
  if(length(statistic$size)!=length(uvo$size)){
    stop('Test sample group sizes and null table group sizes do not match - not same number of groups')
  }
  if(any(sort(statistic$size)!=sort(uvo$size))){
    stop('Test sample group sizes and null table group sizes do not match!')
  }
  if(statistic$variant != uvo$variant){
    stop(paste0('statistic variant is ',statistic$variant,' while null table variant is ', uvo$variant))
  }
  if(statistic$aggregation.type != uvo$aggregation.type){
    stop(paste0('statistic aggregation.type is ',statistic$aggregation.type,' while null table aggregation.type is ', uvo$aggregation.type))
  }
  if(statistic$score.type != uvo$score.type){
    stop(paste0('statistic score.type is ',statistic$score.type,' while null table score.type is ', uvo$score.type))
  }
  if(statistic$type != uvo$type){
    stop(paste0('statistic type is ',statistic$type,' while null table type is ', uvo$type))
  }
  
  if(is.element(statistic$variant,c('KSample-Equipartition'))){
    if(statistic$nr.atoms != uvo$nr.atoms){
      stop(paste0('Null Table and Statistic not of same nr.atoms of KSample-Equipartition variant. Null table is: ',uvo$nr.atoms,' while statistic object is: ', statistic$nr.atoms))
    }
  }

  ret=list()
  if(statistic$stat.type == 'KSample-Stat'){
    if(uvo$minm >m | max(uvo$maxm) < m){
      stop(paste0('null table does not include number of partitions m: ',m))
    }
    if(statistic$mmin >m | statistic$mmax < m){
      stop(paste0('statistic does not include number of partitions m: ',m))
    }
    #compute pvalue
    m.ind = m - uvo$minm +1
    m.stat.ind = m - statistic$mmin +1
    stat=statistic$statistic[m.stat.ind] #the single m statistic value under the statistic object
    return(UnivariateObject.findpval(uvo,stat,m.ind))#return(uvo$m.stat.negative.ecdf[[m.ind]](-1 * stat))
  }else{
    stop('incorrect statistic type: statistic should be the result of hhg.univariate.ks.stat')
  }  
}

#'@title k-sample DS and Mm test. 
#'@description function for performin k sample DS and Mm test. this function is called by the M.m Test Function
#'@param x vector of sample values
#'@param y vector of sample numbers
#'@param variant should be equal 'ds' or 'mds'. if 'ds' then lambda is valid. on 'mds' prior is valid 
#'@param lambda for 'ds'
#'@param Max.m The maximum m partition to be used.
#'@param m.stats.wanted Boolean value indicating whether to return the list of Mm statistics.
#'@param prior vector of prior/penalty function to be used
.hhg.univariate.ks.stat.max = function(x, y, variant = 'Mm', lambda = 1,Max.m=2,m.stats.wanted=F,prior=NA,DS_type=1,nr_equipartition_bins = nr_bins_equipartition(length(x))){
  w.max = 0
  w.sum = 2
  tables.wanted = F
  perm.stats.wanted = F
  nr.threads = 1
  
  if (variant == 'ds') {
    test_type = .UV_KS_DS
  } else if (variant == 'Mm') {
    test_type = .UV_KS_MDS
  } else {
    stop('Unexpected variant specified')
  }
  
  if(DS_type == 2){
    Max.m = min(Max.m, nr_equipartition_bins)
  }
  
  nr.perm = 0
  is.sequential = F
  total.nr.tests = 1
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  is_sequential = as.integer(is.sequential)  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  # y is passed as numbers in 0:(K - 1), sorted according to the order of x
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  }
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  y = as.matrix(y[order(x)])
  
  # Dx and Dy are not used
  Dx = 0
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  #fixme:  this assumes that prior length is at least maxk in length
  prior_param=F
  prior_length_param=1
  if(!is.na(prior[1])){
    prior_param=  prior
    prior_length_param=length(prior)
  }
  else{
    prior_param=rep(0,Max.m)  
    prior_length_param = Max.m
  }
  if(prior_length_param <Max.m-2){
    stop('prior supplied is not of sufficient length')
  }
  
  extra_params = c(as.double(DS_type),as.double(nr_equipartition_bins), as.double(lambda),as.double(Max.m+0.1),as.double(prior_length_param+0.1),as.double(prior_param))#0.1 is in order to make sure after int casting we get the right number
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(y), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F,m.stats.wanted =m.stats.wanted ,m.stats.size = 2*(Max.m-1)) #changed size of n from nrow(Dx) to y
  
  return (ret)
}

# The k-sample version of the XDP test
.xdp.test.k.sample.mk = function(x, y, ddp.K = 3, w.sum = 0, w.max = 2,equipartition.type = 0, equipartition.nr.bins = nr_bins_equipartition(length(x))) 
{
  # The interface may need more work..
  
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (ddp.K != as.integer(ddp.K) || ddp.K < 2 || ddp.K > length(x)) {
    stop('K must be an integer between 2 and length(x)')
  } 
  
  
  if(!(equipartition.type ==0 || equipartition.type ==1)){
    stop('equipartition type for XDP KS should be 0 or 1')
  }
  
  if(!is.numeric(equipartition.nr.bins)){
    stop('equipartition.nr.bins should be a numeric greater than 4 and smaller than sample size')
  }
  
  if(equipartition.nr.bins<4 || equipartition.nr.bins>length(x)){
    stop('equipartition.nr.bins should be a numeric greater than 4 and smaller than sample size')
  }
  
  
  test_type = .UV_KS_XDP_MK
  
  
  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  
  # y is passed as numbers in 0:(K - 1)
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  # Can make these parameter at some point
  nr.perm = 0
  total.nr.tests = 1
  is.sequential = F
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  nr.threads = 1
  tables.wanted = F
  perm.stats.wanted = F  
  
  extra_params = as.double(c(ddp.K,equipartition.type,equipartition.nr.bins))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0,m.stats.wanted = T,m.stats.size=ddp.K*2)
  
  if (ddp.K > 3) {
    ret$max.chisq = NA
    ret$max.lr    = NA
    
    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
  ret_mat=matrix(NA,nrow = ddp.K-1,ncol = 2)
  ret_mat[,1]=ret$m.stats[1:(ddp.K-1)]
  ret_mat[,2]=ret$m.stats[(ddp.K-1) +1:(ddp.K-1)]
  colnames(ret_mat)=c('SC','SL')
  rownames(ret_mat)=paste('m=',as.character(2:(ddp.K)),sep = '')
  return (ret_mat)
}

# The k-sample version of the XDP test
.xdp.test.k.sample = function(x, y, ddp.K = 3, w.sum = 0, w.max = 2) 
{
  # The interface may need more work..
  
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (ddp.K != as.integer(ddp.K) || ddp.K < 2 || ddp.K > length(x)) {
    stop('K must be an integer between 2 and length(x)')
  } 
  
  if (ddp.K == 2) {
    test_type = .UV_KS_XDP2
  } else if (ddp.K == 3) {
    test_type = .UV_KS_XDP3
  } else {
    test_type = .UV_KS_XDP
  }
  
  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  
  # y is passed as numbers in 0:(K - 1)
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  # Can make these parameter at some point
  nr.perm = 0
  total.nr.tests = 1
  is.sequential = F
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  nr.threads = 1
  tables.wanted = F
  perm.stats.wanted = F  
  
  extra_params = as.double(c(ddp.K))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)

  
  if (ddp.K > 3) {
    ret$max.chisq = NA
    ret$max.lr    = NA
    
    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
  
  return (ret)
}

#function for computing ADP for multiple partitions at once
.hhg.test.adp_mk = function(x, y, w.sum = 0, w.max = 2,
                            nr.perm = 0, M = 3,L=M, correct.mi.bias = F, total.nr.tests = 1, 
                            is.sequential = F, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
                            nr.threads = 1, tables.wanted = F, perm.stats.wanted = F,equipartition.type=0 ,equipartitions.nr.bins = nr_bins_equipartition(length(x)))
{
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.ordered(y)) {
    stop('y is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  
  test_type = .UV_IND_ADP_MK
  
  
  # Dx is used to store x
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  y  = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  # For historical reasons, the high-k DDP and ADP variants work on 1-based ranks, while the small-k
  # variants work on 0-based ranks.
  #if (!(variant == 'ddp.obs') && !(variant == 'ddp.all')) {
  #  Dx = Dx - 1
  #  y = y - 1
  #}
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(c(correct.mi.bias,equipartition.type,equipartitions.nr.bins, length(M),M,L ))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret.raw = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0,m.stats.wanted = T,m.stats.size = length(M)*2,debug.vec.needed = F)
  ret=list()
  ret$m=M
  ret$l=L
  ret$sum.lr = ret.raw$m.stats[1:(length(ret.raw$m.stats)/2)+ (length(ret.raw$m.stats)/2)]
  ret$sum.chisq = ret.raw$m.stats[1:(length(ret.raw$m.stats)/2)]
  if("debug.vec" %in% names(ret.raw)){
    ret$debug.vec = ret.raw$debug.vec
  }
  return (ret)
}


#general xdp independence wrapper
.hhg.test.udfree = function(x, y, variant = 'ppr.33.obs', w.sum = 0, w.max = 2,
                            nr.perm = 0, K = 3, correct.mi.bias = F, total.nr.tests = 1, 
                            is.sequential = F, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
                            nr.threads = 1, tables.wanted = F, perm.stats.wanted = F)
{
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.ordered(y)) {
    stop('y is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (K < 2 || K > length(x)) {
    stop('m, number of partitions,  must be an integer greater than 1')
  }
  
  if (variant == 'spr.obs') {
    test_type = .UV_IND_DDP2
  } else if (variant == 'spr.all') {
    test_type = .UV_IND_ADP2
  } else if (variant == 'ppr.22.obs') {
    test_type = .UV_IND_DDP3_C
  } else if (variant == 'ppr.22.all') {
    test_type = .UV_IND_ADP3_C
  } else if (variant == 'ppr.33.obs') {
    test_type = .UV_IND_DDP3
  } else if (variant == 'ppr.33.all') {
    test_type = .UV_IND_ADP3
  } else if (variant == 'tpr.obs') {
    test_type = .UV_IND_DDP4
  } else if (variant == 'tpr.all') {
    test_type = .UV_IND_ADP4
  } else if (variant == 'ddp.obs') {
    test_type = .UV_IND_DDP
  } else if (variant == 'ddp.all') {
    test_type = .UV_IND_ADP
  } else {
    stop('Unexpected variant specified.')
  }
  
  # Dx is used to store x
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  y  = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  # For historical reasons, the high-k DDP and ADP variants work on 1-based ranks, while the small-k
  # variants work on 0-based ranks.
  if (!(variant == 'ddp.obs') && !(variant == 'ddp.all')) {
    Dx = Dx - 1
    y = y - 1
  }
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(c(K, correct.mi.bias))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  return (ret)
}


nr_bins_equipartition = function(n){
  wanted_nr_bins = 60+0.5*floor(sqrt(n))
  return(round(min(wanted_nr_bins,n)))
}

.ml_index_to_chosen=function(chosen_ind,mmax,mmin=2){
  if(mmin!=2){
    stop('mmin != 2 not supported')
  }
  pointer=1
  for(i in mmin:mmax){
    for(j in mmin:mmax){
      if(chosen_ind == pointer){
        return(list(m=i,l=j))
      }
      pointer = pointer + 1
    }
  }
}

.ml_index_to_chosen_vector=function(vec,mmax,mmin=2){
  m=rep(NA,length(vec))
  l=rep(NA,length(vec))
  for(i in 1:length(vec)){
    temp = .ml_index_to_chosen(vec[i],mmax,mmin)
    m[i]=temp$m
    l[i]=temp$l
  }
  return(list(m=m,l=l))
}

lean_null_constructor = function(vec, p0=0.001,p=0.995,p1=0.000001){
  if(!is.numeric(p0) | !is.numeric(p) |!is.numeric(p1)){
    stop('For lean null tables, p0,p and p1 must be numerics between 0 and 1.\n p0 is the resolution for holding the null tabulation up to the p percentile. p1 is the resultion for the p% upper values of the null distribution. p0 must be bigger than p1.')
  }
  if(p0<0 | p0>1 | p<0 | p>1 | p1<0 | p1>1 | p0<p1){
    stop('For lean null tables, p0,p and p1 must be numerics between 0 and 1.\n p0 is the resolution for holding the null tabulation up to the p percentile. p1 is the resultion for the p% upper values of the null distribution. p0 must be bigger than p1.')
  }
  ret=list()
  n=length(vec)
  tab = table(vec)
  sorted_values = sort(vec)
  starting.size = 1000
  mult.factor = 2
  values = rep(NA,starting.size)
  p.higher = rep(NA,starting.size)
  values[0] = 0
  p.higher[0] =  1
  flag=TRUE
  pointer=1
  counter = 2
  while(flag){
    if(counter>length(values)){
      values = c(values,rep(NA,length(values)*(mult.factor-1)))
      p.higher = c(p.higher,rep(NA,length(values)*(mult.factor-1)))
    }
    values[counter] = sorted_values[pointer]
    interval = findInterval(sorted_values[pointer],sorted_values)
    value_counter=0
    value_pointer = pointer
    flag2 = TRUE
    while(flag2){
      if(values[counter] == sorted_values[value_pointer]){
        value_counter = value_counter + 1
        value_pointer = value_pointer - 1
      }else{
        flag2=FALSE
      }
      if(value_pointer<1){
        flag2=FALSE
      }
    }
    p.higher[counter] = (n-interval+value_counter)/(n+1)
    counter =counter+1
    pointer = pointer + ceiling(n*p0)
    if(pointer>=n*p){
      flag=FALSE
    }
  }
  flag=TRUE
  while(flag){
    if(counter>length(values)){
      values = c(values,rep(NA,length(values)*(mult.factor-1)))
      p.higher = c(p.higher,rep(NA,length(values)*(mult.factor-1)))
    }
    values[counter] = sorted_values[pointer]
    interval = findInterval(sorted_values[pointer],sorted_values)
    value_counter=0
    value_pointer = pointer
    flag2 = TRUE
    while(flag2){
      if(values[counter] == sorted_values[value_pointer]){
        value_counter = value_counter + 1
        value_pointer = value_pointer - 1
      }else{
        flag2=FALSE
      }
      if(value_pointer<1){
        flag2=FALSE
      }
    }
    p.higher[counter] = (n-interval+value_counter)/(n+1)
    pointer = pointer + ceiling(n*p1)
    counter = counter+1
    if(pointer>n){
      flag=FALSE
    }
  }
  values[counter] = Inf
  p.higher[counter] = 0
  ret$values = c(0,values[!is.na(values)])
  ret$p.higher = c(1,p.higher[!is.na(p.higher)])
  return(ret)
}

lean_null_p.val = function(value,lean_null_object){
  interval = findInterval(value,lean_null_object$values)
  return(lean_null_object$p.higher[interval])
}

Fast.independence.test.nulltable = function(n,mmin=2,mmax=min(10,n),variant = 'ADP-EQP-ML',nr.atoms = min(40,n),score.type='LikelihoodRatio',nr.perm=200,compress=T, compress.p0=0.001, compress.p=0.99, compress.p1=0.000001){
  if(variant %in% c('DDP','ADP','ADP-ML')){
    stop('Fast.independence.test and Fast.independence.test.nulltable are function with parameters m.max and nr.atoms optimized for large n (>100). For Smaller N and the ADP,DDP and ADP-ML variants, use hhg.univarite.combined.test')
  }
  return(hhg.univariate.ind.nulltable(size = n,mmin = mmin,mmax = mmax,variant = variant,aggregation.type = 'sum',score.type = score.type,nr.replicates = nr.perm,keep.simulation.data = (!compress),nr.atoms = nr.atoms,compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1))
}

Fast.independence.test = function(X,Y,NullTable=NULL,mmin=2, mmax=min(10,length(X)), variant='ADP-EQP-ML',nr.atoms = min(40,length(X)),combining.type='MinP',score.type='LikelihoodRatio',nr.perm=200,compress=T, compress.p0=0.001, compress.p=0.99, compress.p1=0.000001){
  if(variant %in% c('DDP','ADP','ADP-ML') & is.null(NullTable)){
    stop('Fast.independence.test and Fast.independence.test.nulltable are function with parameters m.max and nr.atoms optimized for large n (>100). For Smaller N and the ADP,DDP and ADP-ML variants, use hhg.univarite.combined.test')
  }
  return(hhg.univariate.ind.combined.test(X=X, Y=Y, NullTable = NullTable,
                                          mmin = mmin,mmax = mmax,
                                          variant = variant,aggregation.type = 'sum',score.type = score.type,
                                          combining.type = combining.type,
                                          nr.perm = nr.perm,
                                          nr.atoms = nr.atoms,
                                          compress = compress,compress.p0 = compress.p0,compress.p = compress.p,compress.p1 = compress.p1,keep.simulation.data = (!compress)))
}

.hhg.ddp2.inversions = function(x, y, w.sum = 0, w.max = 2)
{
  #Optimized function for computing ddp2 using inversions
  # Argument checking is redundant, since this function can only be called by the HHG Univariate statistic function for independence
  #if (!is.double(Dx) || !is.double(Dy) || !is.matrix(Dx) || !is.matrix(Dy) || 
  #      nrow(Dx) != ncol(Dx) || nrow(Dx) != nrow(Dy) || nrow(Dy) != ncol(Dy)){
  #  stop('Dx and Dy are expected to be square matrices of doubles, and must have the same number of rows/cols')
  #}
  #if (w.sum < 0 || w.max < 0) {
  #  stop('w.sum and w.max should be greater or equal to zero')
  #}
  #if (nr.perm < 0) {
  #  stop('nr.perm should not be negative')
  #}
  if(length(x)> 30000){
    stop('DDP is not supported for n>30000, Use ADP with atoms \n (ADP-EQP,ADP-EQP-ML) instead.')
  }
  nr.perm = 0 
  is.sequential = F
  seq.total.nr.tests = 1
  seq.alpha.hyp = NULL
  seq.alpha0 = NULL 
  seq.beta0 = NULL
  seq.eps = NULL
  nr.threads = 0
  tables.wanted = F
  perm.stats.wanted = F
  Dx=as.numeric(x)
  Dy=as.numeric(y)
  
  test_type = .UV_IND_OPT_DDP2
  
  
  
  
  dummy.y = matrix(0, length(Dy), 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, dummy.y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  
  class(ret) = 'HHG.Optimal.DDP.Test.Result'
  ret$stat.type = 'Optimal.DDP'
  ret$n = ncol(Dx)
  return (ret)
}

