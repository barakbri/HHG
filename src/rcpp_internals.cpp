#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

/***********************************************************
 *    Rcpp Extension of HHG package                        *
 *    Barak Brill, 2017/10/17,                             *
 *    *****************************                        *
 *    This extension is used to implement                  *
 *    3X3 tests over atoms (rank space subsample)          *
 *    Adding functions here requires adding them           *
 *    to the HHG package routine index                     *
 ***********************************************************/

// [[Rcpp::export]]
NumericMatrix ComputeECDF(NumericVector rank_x, NumericVector rank_y, IntegerVector NR_Atoms){
  int nr_atoms = NR_Atoms(0);
  int n = rank_x.length();
  IntegerVector x_atoms_membership(n);
  IntegerVector y_atoms_membership(n);
  
  
  NumericMatrix A(nr_atoms+1,nr_atoms+1);
  NumericMatrix B(nr_atoms+1,nr_atoms+1);
  //zeroize A,B;
  for(int i=0;i<= nr_atoms;i++ ){
    for(int j=0;j<= nr_atoms;j++ ){
      A(i,j) = 0.0;
      B(i,j) = 0.0;
    }  
  }
  //find atoms membership
  for(int i=0;i<n;i++){
    x_atoms_membership(i) = ceil(((double)(nr_atoms * rank_x(i)))/((double)n));
    y_atoms_membership(i) = ceil(((double)(nr_atoms * rank_y(i)))/((double)n));
    if(x_atoms_membership(i)>= nr_atoms+1){
      //Rprintf("Atom Membershipfix x - max\n\r");
      //Rprintf("x_atoms_membership(i): %d \n\r",x_atoms_membership(i));
      x_atoms_membership(i) = nr_atoms;
    }
    if(y_atoms_membership(i)>= nr_atoms+1){
      //Rprintf("Atom Membershipfix y - max \n\r");
      //Rprintf("y_atoms_membership(i): %d \n\r",y_atoms_membership(i));
      y_atoms_membership(i) = nr_atoms;
    }

    B(x_atoms_membership(i),y_atoms_membership(i)) += 1.0;
  }
  for(int s = 1 ; s <= nr_atoms; s++){
    for(int r = 1 ; r <= nr_atoms; r++){
      A(r,s) = B(r,s-1) + B(r-1,s) - B(r-1,s-1) +B(r,s);
      B(r,s) = A(r,s);
    }
  }
  return A;
}

double compute_obs(NumericMatrix ecdf, int xlow,int xhigh, int ylow, int yhigh){
  double obs =0.0;
  bool debug = false;
  if(debug){
    Rprintf("#########\n\r");
    Rprintf("DEBUG OBS\n\r");
    Rprintf("xlow: %d xhigh: %d ylow: %d yhigh: %d \n\r",xlow,xhigh,ylow,yhigh);
    Rprintf("ecdf(xhigh,yhigh): %lf ecdf(xhigh,ylow): %lf\n\r ecdf(xlow,yhigh): %lf ecdf(xlow,ylow): %lf \n\r",
            ecdf(xhigh,yhigh),ecdf(xhigh,ylow),ecdf(xlow,yhigh),ecdf(xlow,ylow));
    Rprintf("#########\n\r");
  }
  obs = ecdf(xhigh,yhigh) - ecdf(xhigh,ylow) - ecdf(xlow,yhigh) + ecdf(xlow,ylow);
  return obs;
}

double compute_expected( int xlow,int xhigh, int ylow, int yhigh, int n, int nr_atoms){
  
  double T_xlow = floor(((float)xlow) * ((float(n)))/(float(nr_atoms)));
  double T_xhigh = floor(((float)xhigh) * ((float(n)))/(float(nr_atoms)));
  double T_ylow = floor(((float)ylow) * ((float(n)))/(float(nr_atoms)));
  double T_yhigh = floor(((float)yhigh) * ((float(n)))/(float(nr_atoms)));
  double res = (T_xhigh - T_xlow)*(T_yhigh - T_ylow)/((float)n);
  //Rprintf("Expected: %lf \n\r",res);
  return res;
}

// [[Rcpp::export]]
List rcpp_Compute_ADP_MAX_3X3_over_atoms(NumericVector rank_x,NumericVector rank_y,
                                          IntegerVector NR_Atoms,IntegerVector Wmax){
  
  //compute ecdf
  int nr_atoms = NR_Atoms(0);
  int n = rank_x.length();
  double wmax = (double)Wmax(0);
  double min_e = 0.0;
  NumericMatrix ecdf = ComputeECDF(rank_x,rank_y,NR_Atoms);
  
  //define counters
  
  NumericVector max_loglik(1);
  NumericVector max_chisq(1);
  
  
  double current_max_loglik = 0.0;
  double current_max_chisq = 0.0;
  
  NumericVector exp_cells(9);
  NumericVector obs_cells(9);
  
  //splitting location of the selected partition
  NumericVector selected_loglik_i(1);
  NumericVector selected_loglik_j(1);
  NumericVector selected_chisq_i(1);
  NumericVector selected_chisq_j(1);
  NumericVector selected_loglik_i2(1);
  NumericVector selected_loglik_j2(1);
  NumericVector selected_chisq_i2(1);
  NumericVector selected_chisq_j2(1);
  
  //go over partitions
  for(int i = 1; i < nr_atoms-1 ; i++){
    for(int j = i+1; j < nr_atoms ; j++){
      for(int i2 = 1; i2 < nr_atoms-1 ; i2++){
        for(int j2 = i2+1; j2 < nr_atoms ; j2++){
        
        //compute obs for partition:
        obs_cells(0) = compute_obs(ecdf, 0 ,i        ,0  ,i2        );
        obs_cells(1) = compute_obs(ecdf, 0 ,i        ,i2 ,j2        );
        obs_cells(2) = compute_obs(ecdf, 0 ,i        ,j2 ,nr_atoms  );
        obs_cells(3) = compute_obs(ecdf, i ,j        ,0  ,i2        );
        obs_cells(4) = compute_obs(ecdf, i ,j        ,i2 ,j2        );
        obs_cells(5) = compute_obs(ecdf, i ,j        ,j2 ,nr_atoms  );
        obs_cells(6) = compute_obs(ecdf, j ,nr_atoms ,0  ,i2        );
        obs_cells(7) = compute_obs(ecdf, j ,nr_atoms ,i2 ,j2        );
        obs_cells(8) = compute_obs(ecdf, j ,nr_atoms ,j2 ,nr_atoms  );
        
        //compute expected for partition:
        exp_cells(0) = compute_expected( 0, i       , 0  ,i2       , n , nr_atoms);
        exp_cells(1) = compute_expected( 0, i       , i2 ,j2       , n , nr_atoms);
        exp_cells(2) = compute_expected( 0, i       , j2 ,nr_atoms , n , nr_atoms);
        exp_cells(3) = compute_expected( i, j       , 0  ,i2       , n , nr_atoms);
        exp_cells(4) = compute_expected( i, j       , i2 ,j2       , n , nr_atoms);
        exp_cells(5) = compute_expected( i, j       , j2 ,nr_atoms , n , nr_atoms);
        exp_cells(6) = compute_expected( j, nr_atoms, 0  ,i2       , n , nr_atoms);
        exp_cells(7) = compute_expected( j, nr_atoms, i2 ,j2       , n , nr_atoms);
        exp_cells(8) = compute_expected( j, nr_atoms, j2 ,nr_atoms , n , nr_atoms);
        
        
        //computing current scores
        current_max_chisq = 0.0;
        current_max_loglik = 0.0;
        for(int k=0;k<9;k++){
          if(k==0)
            min_e = exp_cells(0);
          else
            min_e = std::min(min_e,exp_cells(k));
          
          current_max_chisq += (obs_cells(k) - exp_cells(k))*(obs_cells(k) - exp_cells(k))/exp_cells(k);
          if(obs_cells(k)>0.0){
            current_max_loglik += obs_cells(k) * log(obs_cells(k)/exp_cells(k));
          }
        }
        //check if new max has been attained
        if(current_max_chisq > max_chisq(0) && min_e>wmax){
          max_chisq(0) = current_max_chisq;
          selected_chisq_i(0) = i;
          selected_chisq_j(0) = j;
          selected_chisq_i2(0) = i2;
          selected_chisq_j2(0) = j2;
        }
        if(current_max_loglik > max_loglik(0)){
          max_loglik(0) = current_max_loglik;
          selected_loglik_i(0) = i;
          selected_loglik_j(0) = j;
          selected_loglik_i2(0) = i2;
          selected_loglik_j2(0) = j2;
        }
        
        
        } // end of j2
      } // end of i2      
    } // end of i
  } // end of j
  
  //return results
  List res   = List::create( max_loglik,
                             max_chisq,
                             selected_loglik_i,
                             selected_loglik_j,
                             selected_chisq_i,
                             selected_chisq_j ,
                             selected_loglik_i2,
                             selected_loglik_j2,
                             selected_chisq_i2,
                             selected_chisq_j2);
  return res;
}

// [[Rcpp::export]]
List rcpp_Compute_ADP_MAX_3X2_over_atoms(NumericVector rank_x,NumericVector rank_y,
                                         IntegerVector NR_Atoms,IntegerVector Wmax){
  
  //compute ecdf
  int nr_atoms = NR_Atoms(0);
  int n = rank_x.length();
  double wmax = (double)Wmax(0);
  double min_e = 0.0;
  NumericMatrix ecdf = ComputeECDF(rank_x,rank_y,NR_Atoms);
  
  //define counters
  
  NumericVector max_loglik(1);
  NumericVector max_chisq(1);
  
  
  double current_max_loglik = 0.0;
  double current_max_chisq = 0.0;
  
  NumericVector exp_cells(6);
  NumericVector obs_cells(6);
  
  //splitting location of the selected partition
  NumericVector selected_loglik_i(1);
  NumericVector selected_loglik_j(1);
  NumericVector selected_chisq_i(1);
  NumericVector selected_chisq_j(1);
  NumericVector selected_loglik_i2(1);
  NumericVector selected_chisq_i2(1);
  
  //go over partitions
  for(int i = 1; i < nr_atoms-1 ; i++){
    for(int j = i+1; j < nr_atoms ; j++){
      for(int i2 = 1; i2 < nr_atoms ; i2++){
        
        //compute obs for partition:
        obs_cells(0) = compute_obs(ecdf, 0 ,i        ,0  ,i2        );
        obs_cells(1) = compute_obs(ecdf, 0 ,i        ,i2 ,nr_atoms  );
        obs_cells(2) = compute_obs(ecdf, i ,j        ,0  ,i2        );
        obs_cells(3) = compute_obs(ecdf, i ,j        ,i2 ,nr_atoms  );
        obs_cells(4) = compute_obs(ecdf, j ,nr_atoms ,0  ,i2        );
        obs_cells(5) = compute_obs(ecdf, j ,nr_atoms ,i2 ,nr_atoms  );
        
        //compute expected for partition:
        exp_cells(0) = compute_expected( 0, i       , 0  ,i2       , n , nr_atoms);
        exp_cells(1) = compute_expected( 0, i       , i2 ,nr_atoms , n , nr_atoms);
        exp_cells(2) = compute_expected( i, j       , 0  ,i2       , n , nr_atoms);
        exp_cells(3) = compute_expected( i, j       , i2 ,nr_atoms , n , nr_atoms);
        exp_cells(4) = compute_expected( j, nr_atoms, 0  ,i2       , n , nr_atoms);
        exp_cells(5) = compute_expected( j, nr_atoms, i2 ,nr_atoms , n , nr_atoms);

        //computing current scores
          current_max_chisq = 0.0;
          current_max_loglik = 0.0;
          for(int k=0;k<6;k++){
            if(k==0)
              min_e = exp_cells(0);
            else
              min_e = std::min(min_e,exp_cells(k));
            
            current_max_chisq += (obs_cells(k) - exp_cells(k))*(obs_cells(k) - exp_cells(k))/exp_cells(k);
            if(obs_cells(k)>0.0){
              current_max_loglik += obs_cells(k) * log(obs_cells(k)/exp_cells(k));
            }
          }
          //check if new max has been attained
          if(current_max_chisq > max_chisq(0) && min_e>wmax){
            max_chisq(0) = current_max_chisq;
            selected_chisq_i(0) = i;
            selected_chisq_j(0) = j;
            selected_chisq_i2(0) = i2;
            
          }
          if(current_max_loglik > max_loglik(0)){
            max_loglik(0) = current_max_loglik;
            selected_loglik_i(0) = i;
            selected_loglik_j(0) = j;
            selected_loglik_i2(0) = i2;
            
          }
          
          
        
      } // end of i2      
    } // end of i
  } // end of j
  
  //return results
  List res   = List::create( max_loglik,
                             max_chisq,
                             selected_loglik_i,
                             selected_loglik_j,
                             selected_chisq_i,
                             selected_chisq_j ,
                             selected_loglik_i2,
                             selected_chisq_i2);
  return res;
}

// [[Rcpp::export]]
List rcpp_Compute_ADP_MAX_2X2_over_atoms(NumericVector rank_x,NumericVector rank_y,
                                         IntegerVector NR_Atoms,IntegerVector Wmax){
  //compute ecdf
  int nr_atoms = NR_Atoms(0);
  int n = rank_x.length();
  double wmax = (double)Wmax(0);
  NumericMatrix ecdf = ComputeECDF(rank_x,rank_y,NR_Atoms);
  
  //define counters
  
  NumericVector max_loglik(1);
  NumericVector max_chisq(1);
  
  
  double current_max_loglik = 0.0;
  double current_max_chisq = 0.0;
  double min_e = 0.0;
  NumericVector exp_cells(4);
  NumericVector obs_cells(4);
  
  //splitting location of the selected partition
  NumericVector selected_loglik_i(1);
  NumericVector selected_chisq_i(1);
  NumericVector selected_loglik_i2(1);
  NumericVector selected_chisq_i2(1);
  
  
  //go over partitions
  for(int i = 1; i < nr_atoms ; i++){
      for(int i2 = 1; i2 < nr_atoms ; i2++){
          //compute obs for partition:
          obs_cells(0) = compute_obs(ecdf, 0 ,i        ,0  ,i2        );
          obs_cells(1) = compute_obs(ecdf, 0 ,i        ,i2 ,nr_atoms  );
          obs_cells(2) = compute_obs(ecdf, i ,nr_atoms ,0  ,i2        );
          obs_cells(3) = compute_obs(ecdf, i ,nr_atoms ,i2 ,nr_atoms  );
          
          
          //compute expected for partition:
          exp_cells(0) = compute_expected( 0, i        , 0  ,i2       , n , nr_atoms);
          exp_cells(1) = compute_expected( 0, i        , i2 ,nr_atoms , n , nr_atoms);
          exp_cells(2) = compute_expected( i, nr_atoms , 0  ,i2       , n , nr_atoms);
          exp_cells(3) = compute_expected( i, nr_atoms , i2 ,nr_atoms , n , nr_atoms);
          
          min_e = exp_cells(0);
          min_e = std::min(min_e,exp_cells(1));
          min_e = std::min(min_e,exp_cells(2));
          min_e = std::min(min_e,exp_cells(3));
          
          //computing current scores
          current_max_chisq = 0.0;
          current_max_loglik = 0.0;
          for(int k=0;k<4;k++){
            current_max_chisq += (obs_cells(k) - exp_cells(k))*(obs_cells(k) - exp_cells(k))/exp_cells(k);
            if(obs_cells(k) > 0.5){
              current_max_loglik += obs_cells(k) * log(obs_cells(k)/exp_cells(k));
            }
          }
          
          //check if new max has been attained
          if(current_max_chisq > max_chisq(0) && min_e>wmax){
            max_chisq(0) = current_max_chisq;
            selected_chisq_i(0) = i;
            selected_chisq_i2(0) = i2;
          }
          if(current_max_loglik > max_loglik(0)){
            max_loglik(0) = current_max_loglik;
            selected_loglik_i(0) = i;
            selected_loglik_i2(0) = i2;
          }
      } // end of i2      
    } // end of i
  
  //return results
  List res   = List::create( max_loglik,
                             max_chisq,
                             selected_loglik_i,
                             selected_chisq_i,
                             selected_loglik_i2,
                             selected_chisq_i2
                             );
  return res;
}
