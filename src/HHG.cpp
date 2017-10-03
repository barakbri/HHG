//============================================================================
// Name        : HHG.cpp
// Author      : Barak Brill &Shachar Kaufman
// Version     : 1.6
// Copyright   : TAU
// Description : An implementation of the Heller-Heller-Gorfine test, ADP and DDP tests and the ADP k-sample variant
//============================================================================

#ifdef WIN32
#include <Windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <pthread.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "HHG.h"

#include "SequentialTest.h"

using namespace std;

// Functions called form R
// ================================================================================================

extern "C" {

SEXP HHG_R_C(SEXP R_test_type, SEXP R_dx, SEXP R_dy, SEXP R_y,
		 SEXP R_w_sum, SEXP R_w_max, SEXP R_resampling_test_params,
		 SEXP R_is_sequential, SEXP R_alpha, SEXP R_alpha0, SEXP R_beta0, SEXP R_eps, SEXP R_nr_perm,
		 SEXP R_nr_threads, SEXP R_tables_wanted, SEXP R_perm_stats_wanted)
{
	try {
		//
		// Parse arguments sent by R, and instantiate the necessary object
		//

		ScoreType score_type = (ScoreType)*INTEGER(R_test_type);
		SEXP Rdim = getAttrib(R_y, R_DimSymbol);

		ResamplingTestConfigurable resampling_test_params(score_type, *REAL(R_w_sum), *REAL(R_w_max), REAL(R_resampling_test_params), *INTEGER(R_nr_perm), *INTEGER(R_is_sequential) != 0, *REAL(R_alpha), *REAL(R_alpha0), *REAL(R_beta0), *REAL(R_eps), 0, *INTEGER(R_nr_threads));
		TestIO test_io(INTEGER(Rdim)[0], INTEGER(Rdim)[1], REAL(R_dx), REAL(R_dy), REAL(R_y), *INTEGER(R_tables_wanted) != 0, *INTEGER(R_perm_stats_wanted) != 0, resampling_test_params);
		SequentialTest seq(test_io, resampling_test_params);

		//
		// Showtime
		//

		seq.run();

		//
		// Cleanup
		//

		test_io.release();

		// Hmm... I can't call the recommended pthread_exit(). It doesn't seem to be critical, unless joining encountered a problem.

		return (test_io.R_output);
	} catch (exception& e) {
		Rprintf(e.what());
		SEXP R_res;
		PROTECT(R_res = allocMatrix(REALSXP, 12, 1));
		UNPROTECT(1);
		return (R_res);
	}
}

}

// TestInput
// ================================================================================================

TestIO::TestIO(int xy_nrow, int y_ncol, double* dx, double* dy, double* y, bool tables_wanted, bool perm_stats_wanted, ResamplingTestConfigurable& resampling_test_params) :
		xy_nrow(xy_nrow), y_ncol(y_ncol), dx(dx), dy(dy), y(y), tables_wanted(tables_wanted), perm_stats_wanted(perm_stats_wanted)
{
	null_dist = z = dz = y; // aliases, used only for readability

	nr_groups = 0;
	y_counts = NULL;

	sorted_dx = NULL;
	sorted_dy = NULL;
	sorted_dz = NULL;
	ranked_dx = NULL;
	ranked_dy = NULL;

	adp = adp_l = adp_r = NULL;
	adp_mk=adp_l_mk=adp_r_mk= NULL;
	

	//FIXME: this is a vicious_hack , and it needs to be redesigned. this is used in order to bring max MDS by k as part of the nr_stats framework.
	ScoreType st = resampling_test_params.score_type;
	if(st == UV_KS_MDS || (st == MV_KS_HHG_EXTENDED && resampling_test_params.uv_score_type == UV_KS_MDS) || st == UV_KS_XDP_MK || st == UV_IND_ADP_MK ){
		k_stats_wanted=true;
	}
	debug_vec_wanted = false;
	// this can be used to activate the debug vector. then, one must write to it in the appropriate function, and parse it out in R.
	if(false){
		debug_vec_wanted = true;
	}
	
	allocate_outputs(resampling_test_params);
	preprocess(resampling_test_params);
}

TestIO::~TestIO() {
	// NOTE: release has to be called explicitly (it's a compromise I'm going to have to live with for now)
}

void TestIO::allocate_outputs(ResamplingTestConfigurable& resampling_test_params) {
	
	ScoreType st = resampling_test_params.score_type;

	int res_pvals_offset = 0;
	int res_obs_stats_offest = resampling_test_params.nr_stats;
	int res_tables_offset = resampling_test_params.nr_stats * 2;
	int res_tables_size = 4 * xy_nrow * xy_nrow * tables_wanted;
	int res_perm_stats_offset = res_tables_offset + res_tables_size;
	int res_perm_stats_size = resampling_test_params.nr_stats * resampling_test_params.nr_perm * perm_stats_wanted;
	int res_k_stats_offset = res_perm_stats_offset+res_perm_stats_size;
	int res_k_stats_size =0;
	
	if(st == UV_KS_MDS){
		res_k_stats_size =  2*(resampling_test_params.Mk_Maxk-1) * k_stats_wanted;
	}
	if(st == UV_KS_XDP_MK){  // then we need to return a double sized k-stats_vec for log likelihood and sum chi
		res_k_stats_size = resampling_test_params.K * 2;
	}
	if(st == UV_IND_ADP_MK){
		res_k_stats_size = resampling_test_params.adp_mk_tables_nr * 2;
	}
	
	int res_debug_vec_offset = res_k_stats_offset +res_k_stats_size;
	int res_debug_vec_size = 0;
	if(debug_vec_wanted){
		res_debug_vec_size = DEBUG_VEC_SIZE; // currently this is hardcoded, might be changed later...
	}
	
	
	
	int res_size = resampling_test_params.nr_stats * 2 + res_tables_size + res_perm_stats_size+res_k_stats_size+res_debug_vec_size;

	PROTECT(R_output = allocMatrix(REALSXP, res_size, 1));
	double* res = REAL(R_output);
	
	

	obs_stats = res + res_obs_stats_offest;
	obs_tbls = res + res_tables_offset;
	pvals = res + res_pvals_offset;
	perm_stats = res + res_perm_stats_offset;
	k_stats=res+res_k_stats_offset;
	debug_vec= res + res_debug_vec_offset;
	
	if(debug_vec_wanted){
		for(int i=0;i<DEBUG_VEC_SIZE;i++){
			debug_vec[i] = NA_REAL;
		}
	}
	
	if (tables_wanted) {
		for (int i = 0; i < 4 * xy_nrow * xy_nrow; ++i) {
			obs_tbls[i] = NA_REAL;
		}
	}
}

void TestIO::preprocess(ResamplingTestConfigurable& resampling_test_params) {
	ScoreType st = resampling_test_params.score_type;

	if (IS_KS_TEST(st)) {
		count_unique_y();
	} else {
		nr_groups = 1;
	}

	if (st == MV_TS_HHG || st == MV_TS_EXISTING || st == MV_KS_HHG || st == MV_KS_HHG_EXTENDED || st == MV_IND_HHG_NO_TIES || st == MV_IND_HHG || st == MV_IND_HHG_EXTENDED ) {
		sort_x_distances_per_row();
	}

	if (st == MV_IND_HHG_NO_TIES || st == MV_IND_HHG || st == MV_IND_HHG_EXTENDED) {
		sort_y_distances_per_row();
	}

	if (IS_CI_MVZ_TEST(st)) {
		sort_z_distances_per_row();
	}

	if (st == MV_KS_HHG_EXTENDED || st == MV_IND_HHG_EXTENDED) {
		rank_x_distances_per_row();
	}

	if (st == MV_IND_HHG_EXTENDED) {
		rank_y_distances_per_row();
	}
	
	if(st == UV_IND_OPT_DDP2 || st == UV_IND_OPT_HOEFFDING){
		sort_x_distances_opt();
		sort_y_distances_opt();
	}

	if (st == UV_IND_ADP || st == CI_UDF_ADP_MVZ_NN) {
		declare_adp_independence(xy_nrow, resampling_test_params.K);
		compute_adp_independence(xy_nrow, resampling_test_params.K);
	}else if(st == UV_IND_ADP_MK){
		is_equipartition = resampling_test_params.equipartition_type;
		equipartition_m_nr_bins = resampling_test_params.equipartition_nr_cells_m;
		if(is_equipartition ==1){
			choose_nr_atoms = equipartition_m_nr_bins;
		}else{
			choose_nr_atoms = xy_nrow;
		}
		declare_adp_independence(choose_nr_atoms, resampling_test_params.Mk_Maxk);
		declare_adp_independence_mk(choose_nr_atoms, resampling_test_params.Mk_Maxk);
		compute_adp_independence_mk(choose_nr_atoms, resampling_test_params.Mk_Maxk);
	}else if (st == MV_IND_HHG_EXTENDED && resampling_test_params.uv_score_type == UV_IND_ADP) {
		declare_adp_independence(xy_nrow, resampling_test_params.K);
		compute_adp_independence(xy_nrow - 1, resampling_test_params.K);
	} else if (IS_UV_DF_KS_TEST(st) || IS_UV_DF_GOF_TEST(st)) {
	
		if(st==UV_KS_XDP_MK){
			is_equipartition = resampling_test_params.equipartition_type;
			equipartition_m_nr_bins = resampling_test_params.equipartition_nr_cells_m;
			if(is_equipartition ==1){
				choose_nr_atoms = equipartition_m_nr_bins;
			}else{
				choose_nr_atoms = xy_nrow;
			}
			declare_adp_k_sample(choose_nr_atoms, resampling_test_params.K);
			declare_adp_k_sample_mk(choose_nr_atoms, resampling_test_params.K);
			compute_adp_k_sample_mk(choose_nr_atoms, resampling_test_params.K);
		}else{
			declare_adp_k_sample(xy_nrow, resampling_test_params.K);
			compute_adp_k_sample(xy_nrow, resampling_test_params.K);
		}
		
	} else if (st == MV_KS_HHG_EXTENDED && resampling_test_params.uv_score_type == UV_KS_XDP) {
		declare_adp_k_sample(xy_nrow - 1, resampling_test_params.K);
		compute_adp_k_sample(xy_nrow - 1, resampling_test_params.K);		
	} 
}

void TestIO::release(void) {
	// NOTE: these deletes rely on the fact that the standard "delete" just
	// does nothing when pointer is null
	delete[] y_counts;

	delete sorted_dx;
	delete sorted_dy;
	delete sorted_dz;

	delete[] ranked_dx;
	delete[] ranked_dy;

	delete[] adp;
	delete[] adp_l;
	delete[] adp_r;
	
	delete[] adp_mk;
	delete[] adp_l_mk;
	delete[] adp_r_mk;
	

	UNPROTECT(1);
}

void TestIO::count_unique_y(void) {
	// NOTE: this assumes y are in 0...K-1

	nr_groups = 0;
	for (int i = 0; i < xy_nrow; ++i) {
		nr_groups = max(nr_groups, (int)(y[i]));
	}
	nr_groups = max(2, nr_groups + 1);

	y_counts = new int[nr_groups];
	for (int k = 0; k < nr_groups; ++k) {
		y_counts[k] = 0;
	}

	for (int i = 0; i < xy_nrow; ++i) {
		++y_counts[(int)(y[i])];
	}
}


void TestIO::sort_x_distances_opt(void){
	sorted_dx = new dbl_int_pair_matrix;
	sorted_dx->resize(1);

	int k=0;
	(*sorted_dx)[k].resize(xy_nrow);

	for (int l = 0; l < xy_nrow; ++l) {
		(*sorted_dx)[k][l].first = dx[l];
		(*sorted_dx)[k][l].second = l;
	}

	sort((*sorted_dx)[k].begin(), (*sorted_dx)[k].end(), dbl_int_pair_comparator);
	
}

void TestIO::sort_y_distances_opt(void){
	sorted_dy = new dbl_int_pair_matrix;
	sorted_dy->resize(1);

	// Could be parallelized as well

	int k=0;
	(*sorted_dy)[k].resize(xy_nrow);

	for (int l = 0; l < xy_nrow; ++l){
		(*sorted_dy)[k][l].first = dy[l];
		(*sorted_dy)[k][l].second = l;
	}

	sort((*sorted_dy)[k].begin(), (*sorted_dy)[k].end(), dbl_int_pair_comparator);
	
}

void TestIO::sort_x_distances_per_row(void) {
	sorted_dx = new dbl_int_pair_matrix;
	sorted_dx->resize(xy_nrow);

	// This could be parallelized - trivial to do per row, but currently
	// not necessary since we parallelize whole test calls.

	for (int k = 0; k < xy_nrow; ++k) {
		(*sorted_dx)[k].resize(xy_nrow);

		for (int l = 0; l < xy_nrow; ++l) {
			(*sorted_dx)[k][l].first = dx[l * xy_nrow + k];
			(*sorted_dx)[k][l].second = l;
		}

		sort((*sorted_dx)[k].begin(), (*sorted_dx)[k].end(), dbl_int_pair_comparator);
	}
}

// (no idea why this code is replicated...)

void TestIO::sort_y_distances_per_row(void) {
	sorted_dy = new dbl_int_pair_matrix;
	sorted_dy->resize(xy_nrow);

	// Could be parallelized as well

	for (int k = 0; k < xy_nrow; ++k) {
		(*sorted_dy)[k].resize(xy_nrow);

		for (int l = 0; l < xy_nrow; ++l){
			(*sorted_dy)[k][l].first = dy[l * xy_nrow + k];
			(*sorted_dy)[k][l].second = l;
		}

		sort((*sorted_dy)[k].begin(), (*sorted_dy)[k].end(), dbl_int_pair_comparator);
	}
}

void TestIO::sort_z_distances_per_row(void) {
	// NOTE: We don't actually need to sort dz, only for each row i we list
	// any index j which is in the z-neighborhood of i.
	// The current implementation is O(n*log(n)) per row. It can also be
	// done in O(n*K) with K the neighborhood size, but this can be up to O(n^2)
	// for large K. Is it doable in O(n)?

	sorted_dz = new dbl_int_pair_matrix;
	sorted_dz->resize(xy_nrow);

	// Could be parallelized as well

	for (int k = 0; k < xy_nrow; ++k) {
		(*sorted_dz)[k].resize(xy_nrow);

		for (int l = 0; l < xy_nrow; ++l) {
			(*sorted_dz)[k][l].first = dz[l * xy_nrow + k];
			(*sorted_dz)[k][l].second = l;
		}

		sort((*sorted_dz)[k].begin(), (*sorted_dz)[k].end(), dbl_int_pair_comparator);
	}
}

// FIXME perhaps the ranking I implemented is too simplistic. What if there are tied distances?

void TestIO::rank_x_distances_per_row(void) {
	ranked_dx = new int[xy_nrow * xy_nrow]; // rank(dx) per row
	for (int i = 0; i < xy_nrow; ++i) {
		for (int j = 0; j < xy_nrow; ++j) {
			ranked_dx[(*sorted_dx)[i][j].second * xy_nrow + i] = j + 1;
		}
	}
}

void TestIO::rank_y_distances_per_row(void) {
	ranked_dy = new int[xy_nrow * xy_nrow]; // rank(dy) per row
	for (int i = 0; i < xy_nrow; ++i) {
		for (int j = 0; j < xy_nrow; ++j) {
			ranked_dy[(*sorted_dy)[i][j].second * xy_nrow + i] = j + 1;
		}
	}
}

void TestIO::compute_adp_independence(int n, int K) {

	#if 0 // partition at ranks
	for (int xl = 1; xl <= n; ++xl) {
		for (int xh = xl; xh <= n; ++xh) {
			int idx = (xl - 1) * xy_nrow + xh - 1;

			if (xl == 1) {
				// left anchored interval
				adp[idx] = my_choose(n - xh - 2 - (K - 2), K - 2);
			} else if (xh == n) {
				// right anchored interval
				adp[idx] = my_choose(xl - 3 - (K - 2), K - 2);
			} else if (xl == 2 || xh == n - 1 || K == 2) {
				adp[idx] = 0;
			} else if (K == 3) {
				adp[idx] = 1;
			} else {
				// nondegenerate regular interval
				adp[idx] = 0;
				for (int i = 0; i <= K - 3; ++i) {
					adp[idx] += my_choose(xl - 3 - i, i) * my_choose(n - xh - 2 - (K - 3 - i), K - 3 - i);
				}
			}
		}
	}
#else
	
	// left anchored interval (xl == 1)
	for (int xh = 1; xh <= n; ++xh) {
		adp_l[xh - 1] = my_choose(n - xh - 1, K - 2); // (the xh == n case is irrelevant)
	}

	// right anchored interval (xh == n)
	for (int xl = 1; xl <= n; ++xl) {
		adp_r[xl - 1] = my_choose(xl - 2, K - 2); // (the xl == 1 case is irrelevant)
	}

	// nondegenerate regular interval
	for (int xd = 0; xd < n; ++xd) {
		adp[xd] = my_choose(n - xd - 3, K - 3); // (the xd > n - 3 case is irrelevant)
	}
#endif
}

void TestIO::compute_adp_independence_mk_single(int n, int K) {
	
	double log_denom = my_lchoose(n - 1, K - 1);
	
	// left anchored interval (xl == 1)
	for (int w = 1; w <= n; ++w) {
		if (n - w - 1 < 0 || K - 2 > n - w - 1 || K - 2 <0) {
			adp_l[w-1] = 0;
		}else{
		  adp_l[w-1] = exp(my_lchoose(n - w - 1, K - 2) - log_denom);
		}
	}

	// right anchored interval (xh == n)
	for (int w = 1; w <= n; ++w) {
		if (n - w - 1 < 0 || K - 2 > n - w - 1 || K - 2 <0) {
			adp_r[w-1] = 0;
		}else{
		  adp_r[w-1] = exp(my_lchoose(n - w - 1, K - 2) - log_denom);
		}
	}

	// nondegenerate regular interval
	for (int w = 1; w <= n; ++w) {
		if (n - w - 2 < 0 || K - 3 > n - w - 2 || K - 3 <0) {
		  adp[w-1] = 0;
		}else{
		  adp[w-1] = exp(my_lchoose(n - w - 2, K - 3) - log_denom);
		}
	}
}

void TestIO::compute_adp_k_sample(int n, int K) {
	
	// NOTE: in this test we would like to be able to use large samples, and the involved
	// binomial coefficients can get much larger than the DOUBLE_XMAX. We thus need to
	// take extra care to use logarithm binomial coefficients.
	double log_denom = my_lchoose(n - 1, K - 1);
	
	// edge-anchored interval (xi == 0)
	for (int w = 1; w < n; ++w) {
		if(n-w-1 < 0  || K-2 <0 || n- w - 1 < K-2){
		adp_l[w] = 0;
		}else{
		adp_l[w] = exp(my_lchoose(n - w - 1, K - 2) - log_denom);
		}
		
	}
	
	// mid interval
	for (int w = 1; w <= n - 2; ++w){
		if( n -w -2 <0 || K - 3 < 0 || n-w-2 <K-3){
			adp[w] = 0;
		}else{
			adp[w] = exp(my_lchoose(n - w - 2, K - 3) - log_denom);
		}

	}
	
}

void TestIO::compute_adp_k_sample_mk(int n, int K) { //compute adp choose constants for all k up to K
	
	for(int i=0; i<n*(K-1)+1;i++){
		adp_mk[i] = 0;
		adp_l_mk[i] = 0;
	}
	
	for (int i=0;i<K-1;i++){
		TestIO::compute_adp_k_sample(n, i+2); //our nr. partitions is from 2 to K
		for(int j=1;j<n;j++){
			adp_mk[i*n+j]=adp[j];
			adp_l_mk[i*n+j]=adp_l[j];
		}
		
		
	}
}


void TestIO::compute_adp_independence_mk(int n, int K) { //compute adp choose constants for all k up to K
	
	
	for(int i=0; i<n*(K-1)+1;i++){
		adp_mk[i] = 0;
		adp_l_mk[i] = 0;
		adp_r_mk[i] = 0;
	}
	
	for (int i=0;i<K-1;i++){
		TestIO::compute_adp_independence_mk_single(n, i+2 ); //our nr. partitions is from 2 to K
		for(int j=0;j<n-1;j++){
			adp_mk[i*n+j]=adp[j];
			adp_l_mk[i*n+j]=adp_l[j];
			adp_r_mk[i*n+j]=adp_r[j];
		}

	}
}

	void TestIO::declare_adp_independence(int n,int K){
		adp   = new double[n];
		adp_l = new double[n];
		adp_r = new double[n];
		
		for(int i=0;i<n;i++){
		  adp[i]=0;
		  adp_l[i]=0;
		  adp_r[i]=0;
		}
	}
	
	void TestIO::declare_adp_independence_mk(int n,int K){
		adp_mk = new double[n*(K-1)+1];
		adp_l_mk = new double[n*(K-1)+1];
		adp_r_mk = new double[n*(K-1)+1];
		
		
		for(int i=0;i<n*(K-1)+1;i++){
		  adp_mk[i]=0;
		  adp_l_mk[i]=0;
		  adp_r_mk[i]=0;
		}
		
	}
	
	void TestIO::declare_adp_k_sample(int n,int K){
		adp = new double[n];
		adp_l = new double[n];
		
		
		for(int i=0;i<n;i++){
		  adp[i]=0;
		  adp_l[i]=0;
		}
		
		
	}
	
	void TestIO::declare_adp_k_sample_mk(int n,int K){
			adp_mk = new double[n*(K-1)+1];
			adp_l_mk = new double[n*(K-1)+1];
			
			
			for(int i=0;i<n*(K-1)+1;i++){
			  adp_mk[i]=0;
			  adp_l_mk[i]=0;
			}
			
	}



double TestIO::my_choose(int n, int k) {
    if (n < 0) {
        return (0);
    }
	return (choose(n, k));
}


double TestIO::my_lchoose(int n, int k) {
    if (n < 0) {
        return (0);
    }
	if(k>n){
		return(0);
	}
	if(k<0){
		return(0);
	}
	return (lchoose(n, k));
}

// ResamplingTestParams
// ================================================================================================

ResamplingTestConfigurable::ResamplingTestConfigurable(ScoreType score_type, double w_sum, double w_max, double* score_params_r,
		int nr_perm, bool is_sequential, double alpha, double alpha0, double beta0, double eps, int base_seed, int nr_threads) :
		ScoreConfigurable(score_type, w_sum, w_max, score_params_r),
		nr_perm(nr_perm), is_sequential(is_sequential), alpha(alpha), alpha0(alpha0), beta0(beta0), eps(eps), base_seed(base_seed), nr_threads(nr_threads)
{
	if (this->nr_threads == 0) {
		this->nr_threads = get_available_nr_threads();
	}

	// NOTE: we round upward the number of permutations that the user specified,
	// to be a multiple of the number of threads. Package R code will later
	// discard any surplus permutations

	nr_perm_per_thread = ceil(double(this->nr_perm) / this->nr_threads);
	this->nr_perm = nr_perm_per_thread * this->nr_threads;
}

int ResamplingTestConfigurable::get_available_nr_threads(void) {
	// This could probably be improved, in particular made more general and accurate
#ifdef NO_THREADS
	return (1);
#else
#ifdef _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return (sysinfo.dwNumberOfProcessors);
#else
	return (sysconf(_SC_NPROCESSORS_ONLN));
#endif
#endif
}

// ScoreParams
// ================================================================================================

ScoreConfigurable::ScoreConfigurable(ScoreType score_type, double w_sum, double w_max, double* extra_params_r) :
		score_type(score_type), w_sum(w_sum), w_max(w_max)
{
	K = 0;
	correct_mi_bias = false;
	sig = 0;
	lambda = 0;
	Mk_Maxk=0;
	prior_length=0;
	prior=NULL;
	adp_mk_tables_nr=0;
	adp_mk_tables_m = NULL;
	adp_mk_tables_l = NULL;
	nnh = 0;
	nnh_lsb = 0;
	nnh_grid_cnt = 0;
	nnh_grid = NULL;
	uv_score_type = UV_GOF_WXN;
	nr_stats = 0;

	parse_params(score_type, extra_params_r);
}

ScoreConfigurable::~ScoreConfigurable() {
	// nothing for now
}

void ScoreConfigurable::parse_params(ScoreType st, double* extra_params_r) {
	nr_stats += BASE_NR_STATS;

	if (st == UV_IND_DDP || st == UV_IND_ADP) {
		K = extra_params_r[0];
		correct_mi_bias = extra_params_r[1];
	} else if(st == UV_IND_ADP_MK){
		correct_mi_bias = extra_params_r[0];
		equipartition_type = (int)extra_params_r[1];
		equipartition_nr_cells_m = (int)extra_params_r[2];
		adp_mk_tables_nr = (int) extra_params_r[3];
		adp_mk_tables_m = new int[adp_mk_tables_nr];
		adp_mk_tables_l = new int[adp_mk_tables_nr];
		int  pointer = 4;
		for(int i=0;i<adp_mk_tables_nr;i++){
			adp_mk_tables_m[i]=(int) extra_params_r[pointer];
			pointer ++;
			if(Mk_Maxk < adp_mk_tables_m[i]){
			Mk_Maxk = adp_mk_tables_m[i];
			}
		}
		for(int i=0;i<adp_mk_tables_nr;i++){
			adp_mk_tables_l[i]=(int) extra_params_r[pointer];
			pointer ++;
			if(Mk_Maxk < adp_mk_tables_l[i]){
			Mk_Maxk = adp_mk_tables_l[i];
			}
		}
	}else if (st == UV_KS_DS ||st == UV_KS_MDS ) {
		equipartition_type = (int)(extra_params_r[0]);
		equipartition_nr_cells_m = (int)(extra_params_r[1]);
		lambda = extra_params_r[2];
		Mk_Maxk=(int)( extra_params_r[3]);
		prior_length=(int)(extra_params_r[4]);
		prior= new double[prior_length];
		for(int i=0;i<prior_length;i++){
			prior[i]=extra_params_r[i+5];
		}
	} else if (st == CI_UVZ_NN || st == CI_MVZ_NN) {
		nnh = extra_params_r[0];
		nnh_lsb = extra_params_r[1];
	} else if (st == CI_UVZ_GAUSSIAN || st == CI_MVZ_GAUSSIAN) {
		sig = extra_params_r[0];
		nnh_lsb = extra_params_r[1];
	} else if (st == CI_UDF_ADP_MVZ_NN) {
		K = extra_params_r[0];
		correct_mi_bias = extra_params_r[1];
		nnh = extra_params_r[2];
		nnh_lsb = extra_params_r[3];
	} else if (st == CI_MVZ_NN_GRID_BW) {
		nnh_grid_cnt = extra_params_r[0];
		nnh_lsb = extra_params_r[1];
		nnh_grid = extra_params_r + 2; // NOTE: this means we assume that the parameters remain alive for the duration of the computation
		nr_stats += BASE_NR_STATS * nnh_grid_cnt;
	} else if (IS_UV_KS_XDP_TEST(st) || IS_UV_GOF_XDP_TEST(st)) {
		K = extra_params_r[0];
		equipartition_type = (int)extra_params_r[1];
		equipartition_nr_cells_m = (int)extra_params_r[2];
	} else if (st == MV_KS_HHG_EXTENDED || st == MV_IND_HHG_EXTENDED) {
		uv_score_type = (ScoreType)(extra_params_r[0]);
		parse_params(uv_score_type, extra_params_r + 1);
	} else if (st == UV_GOF_EXISTING) {
		K = 2;
	}

}
