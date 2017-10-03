
.hhg.hoeffding.inversions = function(x, y, w.sum = 0, w.max = 2)
{
  # Argument checking is negligent at this point...
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
  Dx=x
  Dy=y
  
  test_type = .UV_IND_OPT_HOEFFDING
  
  
  
  
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
  
  class(ret) = 'HHG.Optimal.Hoeffding.Test.Result'
  ret$stat.type = 'Optimal.Hoeffding'
  ret$n = ncol(Dx)
  return (ret)
}


.xdp.ks.competitors = function(x, y, nr.perm = 0, total.nr.tests = 1,
                              is.sequential = T, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
                              nr.threads = 1) 
{
  # Can make these parameters at some point
  w.max = 0
  w.sum = 2
  tables.wanted = F
  perm.stats.wanted = F  
  
  test_type = .UV_KS_CVM_KS
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(total.nr.tests, is.sequential, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  
  # y is passed as numbers in 0:(K - 1)
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  }
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(0)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F)
  
  names(ret)[names(ret) == 'sum.chisq'] = 'cvm.chisq'
  names(ret)[names(ret) == 'max.chisq'] = 'ks.chisq'
  names(ret)[names(ret) == 'sum.lr'   ] = 'cvm.lr'
  names(ret)[names(ret) == 'max.lr'   ] = 'ks.lr'
  
  return (ret)
}


