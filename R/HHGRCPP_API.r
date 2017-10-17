ADP_MAX_2X2_statistic = function(x,y,nr_atoms,w.max = 2){
  ranks_x = rank(x)
  ranks_y = rank(y)
  c_res = rcpp_Compute_ADP_MAX_2X2_over_atoms(ranks_x,ranks_y,nr_atoms,w.max)
  ret = list( loglik.max        = c_res[[1]],
              chisq.max         = c_res[[2]],
              loglik.selected.x = c_res[[3]],
              chisq.selected.x  = c_res[[4]],
              loglik.selected.y = c_res[[5]],
              chisq.selected.y  = c_res[[6]]
  )
  return(ret)
}

ADP_MAX_3X3_statistic = function(x,y,nr_atoms,w.max = 2){
  ranks_x = rank(x)
  ranks_y = rank(y)
  c_res = rcpp_Compute_ADP_MAX_3X3_over_atoms(ranks_x,ranks_y,nr_atoms,w.max)
  ret = list( loglik.max        = c_res[[1]],
              chisq.max         = c_res[[2]],
              loglik.selected.low.x = c_res[[3]],
              loglik.selected.high.x = c_res[[4]],
              chisq.selected.low.x  = c_res[[5]],
              chisq.selected.high.x  = c_res[[6]],
              loglik.selected.low.y = c_res[[7]],
              loglik.selected.high.y = c_res[[8]],
              chisq.selected.low.y  = c_res[[9]],
              chisq.selected.high.y  = c_res[[10]]
  )
  return(ret)
}


ADP_MAX_3X2_statistic = function(x,y,nr_atoms,w.max = 2){
  ranks_x = rank(x)
  ranks_y = rank(y)
  c_res = rcpp_Compute_ADP_MAX_3X2_over_atoms(ranks_x,ranks_y,nr_atoms,w.max)
  ret = list( loglik.max        = c_res[[1]],
              chisq.max         = c_res[[2]],
              loglik.selected.low.x = c_res[[3]],
              loglik.selected.high.x = c_res[[4]],
              chisq.selected.low.x  = c_res[[5]],
              chisq.selected.high.x  = c_res[[6]],
              loglik.selected.y = c_res[[7]],
              chisq.selected.y  = c_res[[8]]
            )
  return(ret)
}
