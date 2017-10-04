/*
 * StatsComputer.cpp
 *
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <ctime>
#include <cstring>
#include <pthread.h>
#include <R.h>
#include <Rmath.h>

#include <Rinternals.h>
#ifdef WIN32
#include <stddef.h>
#endif

#include "StatsComputer.h"

using namespace std;

// This class does a very hackish version of polymorphism with respect to scores.
// I tried to redesign it with a proper class hierarchy, but it actually turns out to be
// less elegant and somewhat slower. So for the time being I'm leaving it the way it is.

StatsComputer::StatsComputer(TestIO& test_io, ScoreConfigurable& score_params,pthread_mutex_t *rng_mutex_param) :
		TestIO(test_io), ScoreConfigurable(score_params)
{
	//
	// Allocate and initialize temporary buffers
	//

	ScoreType st = score_type;
	rng_mutex = rng_mutex_param;
	if (IS_UVY_TEST(st)) {
		// A copy of y is created which will be permuted
		// NOTE: in all scores for which this is relevant, ys are univariate integers.
		y_perm = new int[xy_nrow];
		for (int i = 0; i < xy_nrow; ++i) {
			y_perm[i] = y[i];
		}
	} else {
		y_perm = NULL;
	}

	if (st == MV_TS_EXISTING){
		y0_idx = new int[y_counts[0]];
		y1_idx = new int[y_counts[1]];
	}else{
		y0_idx = y1_idx = NULL;
	}

	idx_1_to_n = idx_perm = idx_perm_inv = NULL;

	if (st == MV_IND_HHG_NO_TIES || st == MV_IND_HHG || st == MV_IND_HHG_EXTENDED || st == CI_UVZ_GAUSSIAN || st == CI_MVZ_GAUSSIAN || st == UV_IND_OPT_DDP2 || st == UV_IND_OPT_HOEFFDING) {
		idx_1_to_n = new int[xy_nrow];
		idx_perm = new int[xy_nrow];
		idx_perm_inv = new int[xy_nrow];

		for (int i = 0; i < xy_nrow; ++i) {
			idx_perm[i] = idx_perm_inv[i] = idx_1_to_n[i] = i;
		}
	} else if (st == CI_UVZ_NN || st == CI_MVZ_NN || st == CI_UDF_ADP_MVZ_NN || st == CI_MVZ_NN_GRID_BW) {
		idx_perm = new int[xy_nrow];
		idx_perm_inv = new int[xy_nrow];

		for (int i = 0; i < xy_nrow; ++i) {
			idx_perm[i] = idx_perm_inv[i] = i;
		}
	}

	if (st == UV_GOF_XDP2 || st == UV_GOF_XDP3 || st == UV_KS_XDP2 || st == UV_KS_XDP3 || st == UV_KS_XDP_MK || st == MV_KS_HHG_EXTENDED) {
		tbl_o = new double[nr_groups * K];
		tbl_e = new double[nr_groups * K];
	} else {
		tbl_o = tbl_e = NULL;
	}

	if (st == UV_KS_DS || (st == MV_KS_HHG_EXTENDED && uv_score_type == UV_KS_DS)) {
		ds_ctab = new int*[xy_nrow + 1];
		for (int k = 0; k < xy_nrow + 1; ++k) {
			ds_ctab[k] = new int[nr_groups];
		}

		ds_score = new double[xy_nrow + 1];
		ds_score_pearson = new double[xy_nrow + 1];
		ds_idx = new int[xy_nrow + 1];
		ds_counts = new double[nr_groups];
		
		mds_max_chi_by_k = NULL;
		mds_max_loglikelihood_by_k = NULL;
	} else if (st == UV_KS_MDS || (st == MV_KS_HHG_EXTENDED && uv_score_type == UV_KS_MDS)) {
		ds_ctab = new int*[xy_nrow + 1];
		for (int k = 0; k < xy_nrow + 1; ++k) {
			ds_ctab[k] = new int[nr_groups];
		}
		
		ds_ctab_bins = new int*[xy_nrow + 1];
		for (int k = 0; k < xy_nrow + 1; ++k) {
			ds_ctab_bins[k] = new int[nr_groups];
		}
		

		// FIXME this is wasteful (in space and thus time), we only need half of it, really
		ds_score = new double[(xy_nrow + 1) * (xy_nrow + 1)];
		ds_score_pearson = new double[(xy_nrow + 1) * (xy_nrow + 1)];

		ds_idx = NULL;
		ds_counts = NULL;
		
		mds_max_chi_by_k = new double[xy_nrow-1]; 
		mds_max_loglikelihood_by_k = new double[xy_nrow-1]; 
		
		for(int k=0;k<xy_nrow-1;++k){ // I prefer to get NA's for the ones that are not listed
			mds_max_chi_by_k[k]=NA_REAL; 
			mds_max_loglikelihood_by_k[k]=NA_REAL; 
		}

		
	} else {


		ds_ctab = NULL;
		ds_ctab_bins = NULL;
		ds_idx = NULL;
		ds_score = ds_counts = ds_score_pearson = NULL;
		mds_max_chi_by_k=NULL;
		mds_max_loglikelihood_by_k=NULL;
	}
	if(st==UV_KS_XDP_MK ){ //fixme: this bypasses the MV logic of the test (implemented by shachar).
		xdp_sc_mk = new double[K-1]; //DEBUG_LINUX
		xdp_sl_mk = new double[K-1];
		
		
		for(int k=0;k < K-1;++k){ // I prefer to get NA's for the ones that are not listed
			xdp_sc_mk[k]=NA_REAL; 
			xdp_sl_mk[k]=NA_REAL; 
		}
	}else{
		xdp_sc_mk=NULL;
		xdp_sl_mk=NULL;
	}
	
	if(st == UV_IND_ADP_MK){
		adp_ind_sc_mk = new double[adp_mk_tables_nr];
		adp_ind_sl_mk = new double[adp_mk_tables_nr];
		
		for(int k=0;k<adp_mk_tables_nr;k++){
			adp_ind_sc_mk[k]=NA_REAL;
			adp_ind_sl_mk[k]=NA_REAL;
		}
		
	}else{
		adp_ind_sc_mk = NULL;
		adp_ind_sl_mk = NULL;
	}
	
	
	if (st == MV_IND_HHG_NO_TIES || st == MV_IND_HHG || st == UV_IND_OPT_DDP2 || st == UV_IND_OPT_HOEFFDING) {
		hhg_gen_inversion_count = new int[xy_nrow];
		hhg_gen_source = new int[xy_nrow];
		hhg_gen_xy_perm = new int[xy_nrow];
		hhg_gen_xy_perm_temp = new int[xy_nrow];
		hhg_gen_y_rev = new int[xy_nrow];
		hhg_gen_left_buffer = new int[xy_nrow / 2];
		hhg_gen_right_buffer = new int[xy_nrow / 2 + (xy_nrow & 1)];
		hhg_gen_left_source_buffer = new int[xy_nrow / 2];
		hhg_gen_right_source_buffer = new int[xy_nrow / 2 + (xy_nrow & 1)];
	} else {
		hhg_gen_inversion_count = NULL;
		hhg_gen_source = NULL;
		hhg_gen_xy_perm = NULL;
		hhg_gen_xy_perm_temp = NULL;
		hhg_gen_y_rev = NULL;
		hhg_gen_left_buffer = NULL;
		hhg_gen_right_buffer = NULL;
		hhg_gen_left_source_buffer = NULL;
		hhg_gen_right_source_buffer = NULL;
	}

	if (st == MV_IND_HHG) {
		sorted_dx_gen.resize(xy_nrow);
		for (int k = 0; k < xy_nrow; ++k) {
			sorted_dx_gen[k].resize(xy_nrow);
		}
	}
    
	
	if (st == MV_IND_HHG_EXTENDED || st == MV_KS_HHG_EXTENDED) {
		uvs_n = xy_nrow - 1;
		uvs_x  = new double[uvs_n]; // x distances to point i, from all other points, sorted
		uvs_y  = new double[uvs_n]; // y distances to point i, from all other points, in the order of the sorted dx above
		uvs_xr = new double[uvs_n]; // ranks of the dx above
		uvs_yr = new int   [uvs_n]; // ranks of the dy above
	} else {
		uvs_n = 0;
		uvs_x = uvs_y = uvs_xr = NULL;
		uvs_yr = NULL;
	}

	if (st == MV_KS_HHG_EXTENDED) {
		uvs_yc = new int[nr_groups];
	} else {
		uvs_yc = NULL;
	}

	if (st == UV_IND_DDP || (st == MV_IND_HHG_EXTENDED && uv_score_type == UV_IND_DDP)) {
		x_ordered_by_y = new int[xy_nrow];
		y_ordered_by_x = new int[xy_nrow];
	} else {
		x_ordered_by_y = NULL;
		y_ordered_by_x = NULL;
	}

	dintegral_zero_based_idxs = 0; // only matters for IS_UV_DF_IND_TEST

	if (IS_UV_DF_KS_TEST(st) || st == UV_KS_AD || st == UV_KS_CVM_KS) {
		// OK, so in this case it is not a *double* integral; it is K single integrals
		// the last row will hold the overall x integral (ignoring y)
		dintegral_pn = xy_nrow + 1;
		double_integral = new int[(nr_groups + 1) * dintegral_pn];
	} else if (st == MV_KS_HHG_EXTENDED && (IS_UV_DF_KS_TEST(uv_score_type) || uv_score_type == UV_KS_AD || uv_score_type == UV_KS_CVM_KS)) {
		// as above
		dintegral_pn = xy_nrow;
		double_integral = new int[(nr_groups + 1) * dintegral_pn];
	} else if (IS_UV_DF_IND_TEST(st)) {
		dintegral_pn = xy_nrow + 2;
		dintegral_zero_based_idxs = !(st == UV_IND_DDP || st == UV_IND_ADP || st == UV_IND_ADP_MK );
		if(st == UV_IND_ADP_MK && equipartition_type == 1){
			double_integral = NULL; //since it is not used, and might be impossible to allocate
		}else{
			double_integral = new int[dintegral_pn * dintegral_pn];
		}
	} else if (st == MV_IND_HHG_EXTENDED && IS_UV_DF_IND_TEST(uv_score_type)) {
		dintegral_pn = xy_nrow + 1;
		double_integral = new int[dintegral_pn * dintegral_pn];
	} else if (st == CI_UDF_ADP_MVZ_NN) {
		dintegral_pn = nnh + 2;
		double_integral = new int[dintegral_pn * dintegral_pn];
	} else {
		dintegral_pn = 0;
		double_integral = NULL;
	}
	
	if(st == UV_IND_ADP_MK && equipartition_type == 1){
		dintegral_pn_eqp = equipartition_nr_cells_m + 1; 
		double_integral_eqp = new int[dintegral_pn_eqp * dintegral_pn_eqp];
	}else{
		dintegral_pn_eqp = 0; 
		double_integral_eqp = NULL;
	}

	if (st == UV_KS_KW || (st == MV_KS_HHG_EXTENDED && uv_score_type == UV_KS_KW)) {
		kw_rs = new double[nr_groups];
	} else {
		kw_rs = NULL;
	}
	
	if (st == CI_MVZ_NN_GRID_BW) {
		sum_chi_grid  = new double[nnh_grid_cnt];
		sum_like_grid = new double[nnh_grid_cnt];
		max_chi_grid  = new double[nnh_grid_cnt];
		max_like_grid = new double[nnh_grid_cnt];
	} else {
		sum_chi_grid  = sum_like_grid = max_chi_grid = max_like_grid = NULL;
	}

	if (st == CI_UDF_ADP_MVZ_NN) {
		nn_sorted_x.resize(nnh);
		nn_sorted_y.resize(nnh);
	}
	
	
	//
	// Assign resample and compute handlers
	//

	hhg_extended_uvs = NULL;

	switch (st) {
		case UV_GOF_WXN:
			compute_score = &StatsComputer::uv_gof_wxn;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_GOF_AD:
			compute_score = &StatsComputer::uv_gof_ad;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_GOF_CVM_KS:
			compute_score = &StatsComputer::uv_gof_cvm_ks;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_GOF_DCOV:
			compute_score = &StatsComputer::uv_gof_dcov;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_GOF_XDP2:
			compute_score = &StatsComputer::uv_gof_xdp2;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_GOF_XDP3:
			compute_score = &StatsComputer::uv_gof_xdp3;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_GOF_XDP:
			compute_score = &StatsComputer::uv_gof_xdp;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_KW:
			compute_score = &StatsComputer::uv_ks_kw;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_AD:
			compute_score = &StatsComputer::uv_ks_ad;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_CVM_KS:
			compute_score = &StatsComputer::uv_ks_cvm_ks;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_DCOV:
			compute_score = &StatsComputer::uv_ks_dcov;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_DS:
			compute_score = &StatsComputer::uv_ks_ds;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_MDS:
			compute_score = &StatsComputer::uv_ks_mds;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_XDP2:
			compute_score = &StatsComputer::uv_ks_xdp2;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_XDP3:
			compute_score = &StatsComputer::uv_ks_xdp3;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_KS_XDP:
			compute_score = &StatsComputer::uv_ks_xdp;
			resample = &StatsComputer::resample_univariate;
		break;
		
		case UV_KS_XDP_MK: //DEBUG_LINUX
			compute_score = &StatsComputer::uv_ks_xdp_mk;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_AD:
			compute_score = &StatsComputer::uv_ind_ad;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_CVM_KS:
			compute_score = &StatsComputer::uv_ind_cvm_ks;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_DCOV:
			compute_score = &StatsComputer::uv_ind_dcov;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_DDP2:
			compute_score = &StatsComputer::uv_ind_ddp2;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_DDP3_C:
			compute_score = &StatsComputer::uv_ind_ddp3_c;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_DDP3:
			compute_score = &StatsComputer::uv_ind_ddp3;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_DDP4:
			compute_score = &StatsComputer::uv_ind_ddp4;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_DDP:
			compute_score = &StatsComputer::uv_ind_ddp;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_ADP2:
			compute_score = &StatsComputer::uv_ind_adp2;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_ADP3_C:
			compute_score = &StatsComputer::uv_ind_adp3_c;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_ADP3:
			compute_score = &StatsComputer::uv_ind_adp3;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_ADP4:
			compute_score = &StatsComputer::uv_ind_adp4;
			resample = &StatsComputer::resample_univariate;
		break;

		case UV_IND_ADP:
			compute_score = &StatsComputer::uv_ind_adp;
			resample = &StatsComputer::resample_univariate;
		break;
		
		case UV_IND_ADP_MK:
			compute_score = &StatsComputer::uv_ind_adp_mk;
			resample = &StatsComputer::resample_univariate;
		break;

		case MV_TS_HHG:
			compute_score = &StatsComputer::mv_ts_hhg;
			resample = &StatsComputer::resample_univariate;
		break;

		case MV_KS_HHG:
			compute_score = &StatsComputer::mv_ks_hhg;
			resample = &StatsComputer::resample_univariate;
		break;

		case MV_KS_HHG_EXTENDED:
			compute_score = &StatsComputer::mv_ks_hhg_extended;
			resample = &StatsComputer::resample_univariate;

			switch (uv_score_type) {
				case UV_KS_KW 		: hhg_extended_uvs = &StatsComputer::uvs_ks_kw		; break;
				case UV_KS_AD		: hhg_extended_uvs = &StatsComputer::uvs_ks_ad		; break;
				case UV_KS_CVM_KS	: hhg_extended_uvs = &StatsComputer::uvs_ks_cvm_ks	; break;
				case UV_KS_DCOV		: hhg_extended_uvs = &StatsComputer::uvs_ks_dcov	; break;
				case UV_KS_DS		: hhg_extended_uvs = &StatsComputer::uvs_ks_ds		; break;
				case UV_KS_MDS		: hhg_extended_uvs = &StatsComputer::uvs_ks_mds		; break;
				case UV_KS_XDP2		: hhg_extended_uvs = &StatsComputer::uvs_ks_xdp2	; break;
				case UV_KS_XDP3		: hhg_extended_uvs = &StatsComputer::uvs_ks_xdp3	; break;
				case UV_KS_XDP		: hhg_extended_uvs = &StatsComputer::uvs_ks_xdp		; break;
				default: error("Unexpected univariate k-sample score specified"); hhg_extended_uvs = NULL; //note that uvs_ks_xdp_mk is not supported for this test!
			}
		break;

		case MV_IND_HHG_NO_TIES:
			compute_score = &StatsComputer::mv_ind_hhg_no_ties;
			resample = &StatsComputer::resample_multivariate; // actually y can be univariate here too and I could optimize for such a case
		break;

		case MV_IND_HHG:
			compute_score = &StatsComputer::mv_ind_hhg;
			resample = &StatsComputer::resample_multivariate; // actually y can be univariate here too and I could optimize for such a case
		break;

		case MV_IND_HHG_EXTENDED:
			compute_score = &StatsComputer::mv_ind_hhg_extended;
			resample = &StatsComputer::resample_multivariate;

			switch (uv_score_type) {
				case UV_IND_AD		: hhg_extended_uvs = &StatsComputer::uvs_ind_ad		; break;
				case UV_IND_CVM_KS	: hhg_extended_uvs = &StatsComputer::uvs_ind_cvm_ks	; break;
				case UV_IND_DCOV	: hhg_extended_uvs = &StatsComputer::uvs_ind_dcov	; break;
				case UV_IND_DDP2	: hhg_extended_uvs = &StatsComputer::uvs_ind_ddp2	; break;
				case UV_IND_DDP3_C	: hhg_extended_uvs = &StatsComputer::uvs_ind_ddp3_c	; break;
				case UV_IND_DDP3	: hhg_extended_uvs = &StatsComputer::uvs_ind_ddp3	; break;
				case UV_IND_DDP4	: hhg_extended_uvs = &StatsComputer::uvs_ind_ddp4	; break;
				case UV_IND_DDP		: hhg_extended_uvs = &StatsComputer::uvs_ind_ddp	; break;
				case UV_IND_ADP2	: hhg_extended_uvs = &StatsComputer::uvs_ind_adp2	; break;
				case UV_IND_ADP3_C	: hhg_extended_uvs = &StatsComputer::uvs_ind_adp3_c	; break;
				case UV_IND_ADP3	: hhg_extended_uvs = &StatsComputer::uvs_ind_adp3	; break;
				case UV_IND_ADP4	: hhg_extended_uvs = &StatsComputer::uvs_ind_adp4	; break;
				case UV_IND_ADP		: hhg_extended_uvs = &StatsComputer::uvs_ind_adp	; break;
				default: error("Unexpected univariate independence score specified"); hhg_extended_uvs = NULL;
			}
		break;

		case UV_GOF_EXISTING:
			compute_score = &StatsComputer::uv_gof_existing;
			resample = &StatsComputer::resample_univariate;
		break;

		case MV_TS_EXISTING:
			compute_score = &StatsComputer::mv_ts_existing;
			resample = &StatsComputer::resample_univariate;
		break;

		case MV_KS_EXISTING:
			compute_score = &StatsComputer::mv_ks_existing;
			resample = &StatsComputer::resample_univariate;
		break;

		case MV_IND_EXISTING:
			compute_score = &StatsComputer::mv_ind_existing;
			resample = &StatsComputer::resample_multivariate;
		break;

		case CI_UVZ_NN:
			compute_score = &StatsComputer::ci_uvz_nn;
			resample = &StatsComputer::resample_uvz_ci;
		break;

		case CI_UVZ_GAUSSIAN:
			compute_score = &StatsComputer::ci_uvz_gaussian;
			resample = &StatsComputer::resample_dummy;
		break;

		case CI_MVZ_NN:
			compute_score = &StatsComputer::ci_mvz_nn;
			resample = &StatsComputer::resample_mvz_ci;
		break;

		case CI_MVZ_GAUSSIAN:
			compute_score = &StatsComputer::ci_mvz_gaussian;
			resample = &StatsComputer::resample_dummy;
		break;

		case CI_UDF_ADP_MVZ_NN:
			compute_score = &StatsComputer::ci_udf_adp_mvz_nn;
			resample = &StatsComputer::resample_mvz_ci;
		break;

		case CI_MVZ_NN_GRID_BW:
			compute_score = &StatsComputer::ci_mvz_nn_grid;
			resample = &StatsComputer::resample_mvz_ci;
		break;

		case UV_IND_OPT_DDP2:
			compute_score = &StatsComputer::uv_ind_opt_ddp2;
			resample = &StatsComputer::resample_dummy;
		break;

		case UV_IND_OPT_HOEFFDING:
			compute_score = &StatsComputer::uv_ind_opt_hoeffding;
			resample = &StatsComputer::resample_dummy;
		break;

		
		default:
#ifdef DEBUG_CHECKS
			cerr << "Unexpected test type specified" << endl;
			exit(1);
#endif
			error("Unexpected score specified");
			compute_score = resample = NULL;
		break;
	}

	store_tables = false;
	should_randomize = false;

	min_w = min(w_sum, w_max);

    sum_chi = sum_like = max_chi = max_like = 0;
    max_sum_chi = max_sum_like = sum_max_chi = sum_max_like = 0;
	kahan_c_chi = kahan_c_like = 0;
	ng_chi = ng_like = 0;

	uvs_sc = uvs_mc = uvs_sl = uvs_ml = 0;
	uvs_y0 = 0;
}

StatsComputer::~StatsComputer() {
	delete[] y_perm;

	delete[] y0_idx;
	delete[] y1_idx;

	delete[] idx_1_to_n;
	delete[] idx_perm;
	delete[] idx_perm_inv;

	delete[] hhg_gen_inversion_count;
	delete[] hhg_gen_source;
	delete[] hhg_gen_xy_perm;
	delete[] hhg_gen_xy_perm_temp;
	delete[] hhg_gen_y_rev;
	delete[] hhg_gen_left_buffer;
	delete[] hhg_gen_right_buffer;
	delete[] hhg_gen_left_source_buffer;
	delete[] hhg_gen_right_source_buffer;

	delete[] uvs_x;
	delete[] uvs_y;
	delete[] uvs_xr;
	delete[] uvs_yr;
	delete[] uvs_yc;

	delete[] x_ordered_by_y;
	delete[] y_ordered_by_x;

	delete[] double_integral;
	delete[] double_integral_eqp;
  
	delete[] sum_chi_grid;
	delete[] sum_like_grid;
	delete[] max_chi_grid;
	delete[] max_like_grid;

	delete[] tbl_o;
	delete[] tbl_e;

	delete[] kw_rs;

	delete[] ds_score;
	delete[] ds_score_pearson;
	delete[] ds_idx;
	delete[] ds_counts;
	
	delete[] mds_max_chi_by_k; 
	delete[] mds_max_loglikelihood_by_k;
	delete[] xdp_sc_mk; //DEBUG_LINUX
	delete[] xdp_sl_mk;
	
	delete[] adp_ind_sc_mk;
	delete[] adp_ind_sl_mk;

	if (ds_ctab != NULL) {
		for (int k = 0; k < xy_nrow + 1; ++k) {
			delete[] ds_ctab[k];
		}
	}
	if (ds_ctab_bins != NULL) {
		for (int k = 0; k < xy_nrow + 1; ++k) {
			delete[] ds_ctab_bins[k];
		}
	}
	delete[] ds_ctab;
	delete[] ds_ctab_bins;
}

void StatsComputer::compute(void) {
	store_tables = tables_wanted;
	(this->*compute_score)();
	store_tables = false;
}

void StatsComputer::permute_and_compute(void) {
	// permute the y's (I'm doing it on top of the previous permutation, I don't see a problem with that)
	(this->*resample)();

	// compute statistic
	(this->*compute_score)();
}

// I could do away with this by always working on the stats themselves, but
// it is more convenient to work with local variables and then copy their
// contents to the correct buffer.

void StatsComputer::get_stats(double* stats) {
#if BASE_NR_STATS != 4
	#error ok, so need to update this...
#endif

	stats[0] = sum_chi;
	stats[1] = sum_like;
	stats[2] = max_chi;
	stats[3] = max_like;
	int offst = 4;
	// NOTE: this is safe to call even if not a grid test
	for (int i = 0; i < nnh_grid_cnt; ++i) {
		stats[offst + 0] = sum_chi_grid [i];
		stats[offst + 1] = sum_like_grid[i];
		stats[offst + 2] = max_chi_grid [i];
		stats[offst + 3] = max_like_grid[i];
		offst += 4;
	}
	
	//NOTE: this is writing into the k_stats section , and assumes no other permutations or tables are asked for.
	if (score_type == UV_KS_MDS || (score_type == MV_KS_HHG_EXTENDED && uv_score_type == UV_KS_MDS)) {
		//k_stats[0]=5.0;
		//Mk_Maxk
		//test_io
		//ScoreConfigurable::Mk_Maxk
		int current_max_k = Mk_Maxk;//ScoreConfigurable.
		for (int i = 0; i < current_max_k-1; ++i) {
			k_stats[i] = mds_max_loglikelihood_by_k [i];
			k_stats[i+(current_max_k-1)] = mds_max_chi_by_k [i];
		}
	}
	
	if(score_type == UV_KS_XDP_MK){ //this break the mv logic of the test used by shachar (since again, we pass the statistics through k_stats and not through the usual 4-pipes (SC,SL,MC,ML)
		for (int i = 0; i < (K-1); ++i) { //DEBUG_LINUX
			k_stats[i] = xdp_sc_mk [i];
			k_stats[i+(K-1)] = xdp_sl_mk [i];
		}
	}
	
	if(score_type == UV_IND_ADP_MK){
		for (int i = 0; i < adp_mk_tables_nr; ++i) { 
			k_stats[i] = adp_ind_sc_mk[i];
			k_stats[i+adp_mk_tables_nr] = adp_ind_sl_mk [i];
		}
	}

	if (score_type == MV_IND_HHG_EXTENDED || score_type == MV_KS_HHG_EXTENDED) {
		stats[offst + 0] = sum_max_chi;
		stats[offst + 1] = sum_max_like;
		stats[offst + 2] = max_sum_chi;
		stats[offst + 3] = max_sum_like;
		offst += 4;
	}
}





// Resampling functions
// ================================================================================================

// Generate a new random permutation of the y vector (Fisher-Yates)
// This can probably be optimized and may take a large percent of the CPU time
// Also, rand() is not really thread safe. On Windows it is kind of ok to use but
// on Linux this is a particularly bad idea since the state is process-wide.

void StatsComputer::resample_univariate(void) {
	R_rand_lock();
	for (int i = xy_nrow - 1; i > 0; --i) {
		int j = R_rand_wrapper_nolock() % (i + 1); // should use my_rand() implemented below //DEBUG_RAND

		int temp = y_perm[j];
		y_perm[j] = y_perm[i];
		y_perm[i] = temp;
	}
	R_rand_unlock();
}

void StatsComputer::resample_multivariate(void) {
	// In this case we don't really care about y values, but rather their
	// per (permuted) row sorted (permuted) indices.
	R_rand_lock();
	for (int i = 0; i < xy_nrow; ++i) {
		int j = R_rand_wrapper_nolock() % (i + 1); // should use my_rand() implemented below //DEBUG_RAND

		idx_perm[i] = idx_perm[j];
		idx_perm[j] = i;
	}

	for (int i = 0; i < xy_nrow; ++i) {
		idx_perm_inv[idx_perm[i]] = i;
	}
	R_rand_unlock();
}

void StatsComputer::resample_uvz_ci(void) {
	// A smoothed bootstrap. The smoothing is done over the values of z, so for each
	// sample i we sample uniformly (i.e. currently I am assuming a NN kernel for the
	// smoothing, and even assuming its width is the same as the HHG statistic's
	// kernel) from the z-neighborhood of i. We need to generate resampling indices
	// for x (stored in idx_perm_inv), and independent indices for y (stored in idx_perm).
	// The current implementation also assumes that the observations have been presorted
	// according to their z values (which is only meaningful for univariate z).

	// I guess this could be optimized...

	int bwh = (nnh_lsb >> 1);
	for (int i = 0; i < xy_nrow; ++i) {
		int i1 = max(0, i - bwh);
		int i2 = min(xy_nrow - 1, i + bwh);
		idx_perm    [i] = my_rand(i1, i2);
		idx_perm_inv[i] = my_rand(i1, i2);
	}
}

void StatsComputer::resample_mvz_ci(void) {
	// Same as the univariate version, but has to look in the sorted dz in order to
	// figure out the indices of neighbors

	for (int i = 0; i < xy_nrow; ++i) {
		int nn_x = my_rand(0, nnh_lsb - 1);
		int nn_y = my_rand(0, nnh_lsb - 1);
		idx_perm_inv[i] = (*sorted_dz)[i][nn_x].second;
		idx_perm    [i] = (*sorted_dz)[i][nn_y].second;
	}
}

void StatsComputer::resample_dummy(void) {
	should_randomize = true;
}

// Score compute functions
// ================================================================================================

void StatsComputer::uv_gof_wxn(void) {
	// TODO
}

void StatsComputer::uv_gof_ad(void) {
	// TODO
}

void StatsComputer::uv_gof_cvm_ks(void) {
	// TODO
}

void StatsComputer::uv_gof_dcov(void) {
	// TODO
}

void StatsComputer::uv_gof_xdp2(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_gof_xdp2();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_gof_xdp3(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_gof_xdp3();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_gof_xdp(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_gof_xdp();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ks_kw(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ks_kw();

	sum_chi  = uvs_sc;
	max_chi  = NA_REAL;
	sum_like = NA_REAL;
	max_like = NA_REAL;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ks_ad(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_ad();

	sum_chi  = uvs_sc;
	max_chi  = NA_REAL;
	sum_like = NA_REAL;
	max_like = NA_REAL;

	uvs_xr = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_cvm_ks(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_cvm_ks();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_dcov(void) {
	uvs_n  = xy_nrow;
	uvs_x  = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;
	uvs_y0 = 0; // FIXME it is not clear if this defines a coherent k-sample test, I need to think about it

	uvs_ks_dcov();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_x  = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_ds(void) {
	uvs_n  = xy_nrow;
	uvs_x  = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_ds();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_x  = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_mds(void) {
	uvs_n  = xy_nrow;
	uvs_x  = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_mds();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_x  = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_xdp2(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_xdp2();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_xdp3(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_xdp3();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

// K-sample DDP (also ADP) KxM test
void StatsComputer::uv_ks_xdp(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_xdp();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ks_xdp_mk(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;
	uvs_yc = y_counts;

	uvs_ks_xdp_mk();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
	uvs_yc = NULL;
}

void StatsComputer::uv_ind_ad(void) {
	// TODO
}

void StatsComputer::uv_ind_cvm_ks(void) {
	// TODO
}

void StatsComputer::uv_ind_dcov(void) {
	// TODO
}

void StatsComputer::uv_ind_ddp2(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_ddp2();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_ddp3_c(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_ddp3_c();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_ddp3(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_ddp3();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_ddp4(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_ddp4();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_ddp(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_ddp();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_adp2(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_adp2();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_adp3_c(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_adp3_c();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_adp3(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_adp3();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_adp4(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_adp4();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::uv_ind_adp(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_adp();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}


void StatsComputer::uv_ind_adp_mk(void) {
	uvs_n  = xy_nrow;
	uvs_xr = dx;
	uvs_yr = y_perm;

	uvs_ind_adp_mk();

	sum_chi  = uvs_sc;
	max_chi  = uvs_mc;
	sum_like = uvs_sl;
	max_like = uvs_ml;

	uvs_xr = NULL;
	uvs_yr = NULL;
}

void StatsComputer::mv_ts_hhg(void) {
	int n = xy_nrow;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);
    int i, j, k;
    int y_i, total_same_y, total_same_y_count_so_far, curr_dx_same_y_count;

	sum_chi  = 0;
	sum_like = 0;
	max_chi  = 0;
	max_like = 0;

    for (i = 0; i < n; i++) {
        k = 0;
        y_i = y_perm[i];
        total_same_y = y_counts[y_i]; // NOTE: this assumes that y itself is an index in [0, 1]
        total_same_y_count_so_far = 0;
        curr_dx_same_y_count = 0;

    	for (j = 0; j < n - 1; j++) {
        	k += ((*sorted_dx)[i][k].second == i); // exclude d(i,i)
    		curr_dx_same_y_count += (y_perm[(*sorted_dx)[i][k].second] == y_i);

        	if ((k == n - 1) || ((*sorted_dx)[i][k + 1 + ((*sorted_dx)[i][k+1].second == i)].first > (*sorted_dx)[i][k].first)) {
        		// found all duplicates at current dx

        		if (curr_dx_same_y_count > 0) {
        			// In the notation of the HHG test paper (describing the 2x2 table):
        			//
					// The value "A_1.": j+1 is the number of dx smaller or equal to
        			// the current dx, but we don't want to include the current point
        			// in the count so: j
        			//
					// The value "A_.1": in the two sample setup is
        			// necessarily total_same_y - 2
        			//
					// The value "A_11": total_same_y_count_so_far + curr_dx_same_y_count - 1
        			// where the minus one is for the current point again.
        			//
        			// All remaining table values can be computed from these.
        			//
					// When adding to the sum statistics, add the resulting chisq times-
					// curr_dx_same_y_count, since each such sample has the same table
					// (and all the others have degenerate tables that contribute nothing)

					a00 = total_same_y_count_so_far + curr_dx_same_y_count - 1;
					a01 = j - a00;
					a10 = total_same_y - 2 - a00;
					a11 = n - 2 - j - a10;

					// Note that it is expected this would only be necessary for the computing the
					// observed statistic (the permutation is identity).
					// It might be a better idea to create a separate copy of this function
					// and add this only in the copy.
					if (store_tables) {
						int row = i * n + (*sorted_dx)[i][k].second;
						obs_tbls[        row] = a00;
						obs_tbls[  n*n + row] = a01;
						obs_tbls[2*n*n + row] = a10;
						obs_tbls[3*n*n + row] = a11;
					}

#ifdef DEBUG_CHECKS
					if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
						cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
						exit(1);
					}
#endif

					accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, curr_dx_same_y_count);
        		}

        		// update/reset the counters in preparation for the next unique dx
        		total_same_y_count_so_far += curr_dx_same_y_count;
        		curr_dx_same_y_count = 0;
        	}

        	++k;
		}
    }
}

// NOTE: Actually now this is exactly the same as the 2-sample implementation (just the y_counts here would be of length nr_groups = "K")
void StatsComputer::mv_ks_hhg(void) {
	int n = xy_nrow;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);
    int i, j, k;
    int y_i, total_same_y, total_same_y_count_so_far, curr_dx_same_y_count;

	sum_chi  = 0;
	sum_like = 0;
	max_chi  = 0;
	max_like = 0;

    for (i = 0; i < n; i++) {
#ifdef DEBUG_PRINTS
    	int pi = i; // only for observed

		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "This point has y = " << y_perm[i] << ", which is the same in " << y_counts[(int)y_perm[i]] << " points." << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (k == i);
			cout << j << " (" << k << "): " << dx[k*n+i] << ", " << dy[k*n+i] << endl;
		}
		cout << "Marginally sorted distances dx (and src idx) for this row:" << endl;
		for (j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dx)[i][k].second == i);
			cout << j << ": " << (*sorted_dx)[i][k].first << " (" << (*sorted_dx)[i][k].second << ")" << endl;
		}
		cout << "Marginally sorted distances dy (and src idx) for this row:" << endl;
		for (j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == i);
			cout << j << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

        k = 0;
        y_i = y_perm[i];
        total_same_y = y_counts[y_i]; // NOTE: this assumes that y itself is an index in [0, nr_groups - 1]. Otherwise need to implement y_counts as an STL map or something
        total_same_y_count_so_far = 0;
        curr_dx_same_y_count = 0;

    	for (j = 0; j < n - 1; j++) {
        	k += ((*sorted_dx)[i][k].second == i); // exclude d(i,i)
    		curr_dx_same_y_count += (y_perm[(*sorted_dx)[i][k].second] == y_i);

        	if ((k == n - 1) || ((*sorted_dx)[i][k + 1 + ((*sorted_dx)[i][k+1].second == i)].first > (*sorted_dx)[i][k].first)) {
        		// found all duplicates at current dx

        		if (curr_dx_same_y_count > 0) {
        			// In the notation of the HHG test paper (describing the 2x2 table):
        			//
					// The value "A_1.": j+1 is the number of dx smaller or equal to
        			// the current dx, but we don't want to include the current point
        			// in the count so: j
        			//
					// The value "A_.1": in the two sample setup is
        			// necessarily total_same_y - 2
        			//
					// The value "A_11": total_same_y_count_so_far + curr_dx_same_y_count - 1
        			// where the minus one is for the current point again.
        			//
        			// All remaining table values can be computed from these.
        			//
					// When adding to the sum statistics, add the resulting chisq times-
					// curr_dx_same_y_count, since each such sample has the same table
					// (and all the others have degenerate tables that contribute nothing)

					a00 = total_same_y_count_so_far + curr_dx_same_y_count - 1;
					a01 = j - a00;
					a10 = total_same_y - 2 - a00;
					a11 = n - 2 - j - a10;

					// Note that it is expected this would only be necessary for the computing the
					// observed statistic (the permutation is identity).
					// It might be a better idea to create a separate copy of this function
					// and add this only in the copy.
					if (store_tables) {
						int row = i * n + (*sorted_dx)[i][k].second;
						obs_tbls[        row] = a00;
						obs_tbls[  n*n + row] = a01;
						obs_tbls[2*n*n + row] = a10;
						obs_tbls[3*n*n + row] = a11;
					}

#ifdef DEBUG_PRINTS
					cout << "with point at position " << j << " of the sorted dx: ";
					cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
					if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
						cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
						exit(1);
					}
#endif

					accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, curr_dx_same_y_count);
        		}

        		// update/reset the counters in preparation for the next unique dx
        		total_same_y_count_so_far += curr_dx_same_y_count;
        		curr_dx_same_y_count = 0;
        	}

        	++k;
		}

#ifdef DEBUG_PRINTS
		cout << "current stats: sum_like = " << sum_like << ", sum_chi = " << sum_chi << endl;
#endif
    }
}

void StatsComputer::mv_ks_hhg_extended(void) {
	int n = xy_nrow;
	int ranked_dx_sj, sj, ranked_dx_i;

	sum_chi      = 0;
	max_chi      = 0;
	sum_like     = 0;
	max_like     = 0;
	max_sum_chi  = 0;
	sum_max_chi  = 0;
	max_sum_like = 0;
	sum_max_like = 0;

	for (int i = 0; i < nr_groups; ++i) {
		uvs_yc[i] = y_counts[i];
	}

	for (int i = 0; i < n; ++i) {
		// This will probably be 1, but not necessarily due to the possibility of ties
		ranked_dx_i  = ranked_dx[i * n + i];

		// To make life more modular: first create simple arrays with the relevant n-1 distances in x and their matching y
		for (int j = 0, k = 0; j < n; ++j) {
			sj = (*sorted_dx)[i][j].second;
			if (sj != i) {
				ranked_dx_sj  = ranked_dx[sj * n + i];
				uvs_x[k] = (*sorted_dx)[i][j].first;
				uvs_xr[k] = ranked_dx_sj  - (ranked_dx_sj > ranked_dx_i);
				uvs_yr[k] = y_perm[sj];
				++k;
			}
		}

		// It is also necessary to update the uvs_yc so that they do not count the i'th sample
		// NOTE: I assume that all groups have at least two elements (otherwise I will for
		// some i get a group with zero elements)
		--uvs_yc[y_perm[i]];

		uvs_y0 = y_perm[i];

		// Then compute the univariate score and aggregate
		(this->*hhg_extended_uvs)();

		++uvs_yc[y_perm[i]];

		sum_chi      += uvs_sc;
		max_chi      = max(max_chi, uvs_mc);
		sum_like     += uvs_sl;
		max_like     = max(max_like, uvs_ml);
		sum_max_chi  += uvs_mc;
		sum_max_like +=	uvs_ml;
		max_sum_chi  = max(max_sum_chi, uvs_sc);
		max_sum_like = max(max_sum_like, uvs_sl);
	}

	sum_chi      /= n;
	sum_like     /= n;
	sum_max_chi  /= n;
	sum_max_like /=	n;
}

void StatsComputer::mv_ind_hhg_no_ties(void) {
	int n = xy_nrow, src;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	for (int i = 0; i < n; ++i) {
		int pi = idx_perm[i];

#ifdef DEBUG_PRINTS
		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (k == i);
			cout << j << " (" << k << "): " << dx[k*n+i] << ", " << dy[k*n+i] << endl;
		}
		cout << "Marginally sorted distances dx (and src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dx)[i][k].second == i);
			cout << j << ": " << (*sorted_dx)[i][k].first << " (" << (*sorted_dx)[i][k].second << ")" << endl;
		}
		cout << "Marginally sorted distances dy (and src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == i);
			cout << j << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

		// Use Yair's merge-sort-like implementation (assumes there are no ties)

		// NOTE: I am using the fact that the sorting of permuted y's is the same as
		// the sorting of original y's. I only have to go (a) to the permuted row instead
		// of the i'th row, and (b) pass the "sorted indices" through the permutation.

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == pi);
			src = idx_perm_inv[(*sorted_dy)[pi][k].second]; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_y_rev[src] = j;
		}

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dx)[i][k].second == i);
			src = (*sorted_dx)[i][k].second; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_xy_perm[j] = hhg_gen_y_rev[src];
			hhg_gen_source[j] = j;
			hhg_gen_inversion_count[j] = 0;
			hhg_gen_xy_perm_temp[j] = hhg_gen_xy_perm[j];
		}

		hhg_gen_inversions(hhg_gen_xy_perm_temp, hhg_gen_source, hhg_gen_inversion_count, n - 1);

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			a00 = j - hhg_gen_inversion_count[j];
			a01 = hhg_gen_inversion_count[j];
			a10 = hhg_gen_xy_perm[j] + hhg_gen_inversion_count[j] - j;
			a11 = n - hhg_gen_xy_perm[j] - hhg_gen_inversion_count[j] - 2;

			// Note that it is expected this would only be necessary for the computing the
			// observed statistic (the permutation is identity).
			// It might be a better idea to create a separate copy of this function
			// and add this only in the copy.
			if (store_tables) {
				k += ((*sorted_dx)[i][k].second == i);
				int row = i * n + (*sorted_dx)[i][k].second;
				obs_tbls[        row] = a00;
				obs_tbls[  n*n + row] = a01;
				obs_tbls[2*n*n + row] = a10;
				obs_tbls[3*n*n + row] = a11;
			}

#ifdef DEBUG_PRINTS
			cout << "with point at position " << j << " of the sorted dx: ";
			cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
			if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
				cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
				//exit(1);
			}
#endif

			accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, 1);
		}
	}
}

void StatsComputer::mv_ind_hhg(void) {
	int n = xy_nrow, src;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);

	// First, re-sort first according to x ascending, then according to y descending
	sort_xy_distances_per_row();

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	for (int i = 0; i < n; ++i) {
		int pi = idx_perm[i];

#ifdef DEBUG_PRINTS
		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (k == i);
			cout << j << " (" << k << "): " << dx[k*n+i] << ", " << dy[idx_perm[k]*n+pi] << endl;
		}
		cout << "Lexicographically [x asc, y desc] sorted distances dx (and dy, src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (sorted_dx_gen[i][k].i == i);
			cout << j << ": " << sorted_dx_gen[i][k].x << " (dy = " << sorted_dx_gen[i][k].y << ", src = " << sorted_dx_gen[i][k].i << ")" << endl;
		}
		cout << "Marginally sorted distances dy (and src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == i);
			cout << j << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

		// Use Yair's merge-sort-like implementation (assumes there are no ties)

		// NOTE: I am using the fact that the sorting of permuted y's is the same as
		// the sorting of original y's. I only have to go (a) to the permuted row instead
		// of the i'th row, and (b) pass the "sorted indices" through the permutation.

		double last_new_y = 0;
		int last_new_y_pos = -1;

		for (int j = n - 1, k = n - 1; j >= 1; --j, --k) {
			k -= ((*sorted_dy)[pi][k].second == pi);
			if (last_new_y_pos == -1 || (*sorted_dy)[pi][k].first != last_new_y) {
				last_new_y = (*sorted_dy)[pi][k].first;
				last_new_y_pos = j;
			}
			src = idx_perm_inv[(*sorted_dy)[pi][k].second]; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_y_rev[src] = last_new_y_pos;
		}

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (sorted_dx_gen[i][k].i == i);
			src = sorted_dx_gen[i][k].i; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_xy_perm[j] = hhg_gen_y_rev[src];
			hhg_gen_source[j] = j;
			hhg_gen_inversion_count[j] = 0;
			hhg_gen_xy_perm_temp[j] = hhg_gen_xy_perm[j];
		}

#ifdef DEBUG_PRINTS
		cout << "Contents of y_rev: (should be: rank of y in marginal sorting, listed in original order)" << endl;
		for (int j = 0; j < n - 1; ++j) {
			cout << j << ": " << hhg_gen_y_rev[j] << endl;
		}
		cout << "Contents of xy_perm: (should be: rank of y of points in order of x)" << endl;
		for (int j = 0; j < n - 1; ++j) {
			cout << j << ": " << hhg_gen_xy_perm[j] << endl;
		}
#endif

		hhg_gen_inversions(hhg_gen_xy_perm_temp, hhg_gen_source, hhg_gen_inversion_count, n - 1);

#ifdef DEBUG_PRINTS
		cout << "Contents of inversion_count: (should be: number of inversions in y from left)" << endl;
		for (int j = 0; j < n - 1; ++j) {
			cout << j << ": " << hhg_gen_inversion_count[j] << endl;
		}
#endif

		// (This was an alternative, that I did not pursue)
		// When ties are possible, what we want is, for every point j in the sorted order of x's,
		// (going with j from last to first): how many points have smaller or equal x (that would
		// be the last positions we encountered while walking with j that had a new value) and
		// also have a larger y (this is the number of inversions in y from all those points.
		// this can be decomposed as the inversions in y to the left of j, plus the inversions
		// in y to the right-and-same-x-value. The former we already have, the latter I think we
		// can compute in nlog(n) with a very similar algorithm to "inversions()").

		// (This is what I did)
		// An alternative to the above is not to sort the x's but instead sort according to
		// "x ascending, then y descending" (and so this really has to be redone for every
		// permutation...). If you do this, then the remainder of the algorithm remains unchanged.

		double last_new_x = 0;
		int last_new_x_pos = -1;

		for (int j = n - 2, k = n - 1; j >= 0; --j, --k) {
			k -= ((*sorted_dx)[i][k].second == i);
			if (last_new_x_pos == -1 || (*sorted_dx)[i][k].first != last_new_x) {
				last_new_x = (*sorted_dx)[i][k].first;
				last_new_x_pos = j;
			}

			a00 = last_new_x_pos - hhg_gen_inversion_count[j];
			a01 = hhg_gen_inversion_count[j];
			a10 = hhg_gen_xy_perm[j] - 1 + hhg_gen_inversion_count[j] - last_new_x_pos;
			a11 = n - hhg_gen_xy_perm[j] - hhg_gen_inversion_count[j] - 1;

			// Note that it is expected this would only be necessary for the computing the
			// observed statistic (the permutation is identity).
			// It might be a better idea to create a separate copy of this function
			// and add this only in the copy.
			if (store_tables) {
				int row = i * n + (*sorted_dx)[i][k].second;
				obs_tbls[        row] = a00;
				obs_tbls[  n*n + row] = a01;
				obs_tbls[2*n*n + row] = a10;
				obs_tbls[3*n*n + row] = a11;
			}

#ifdef DEBUG_PRINTS
			cout << "with point at position " << j << " of the sorted dx: ";
			cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
			if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
				cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
				//exit(1);
			}
#endif

			accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, 1);
		}

#ifdef DEBUG_PRINTS
		cout << "current stats: sum_like = " << sum_like << ", sum_chi = " << sum_chi << endl;
#endif
	}
}

void StatsComputer::mv_ind_hhg_extended(void) {
	int n = xy_nrow;
	int sj, psj;
	int ranked_dx_sj, ranked_dy_psj;

	sum_chi      = 0;
	max_chi      = 0;
	sum_like     = 0;
	max_like     = 0;
	max_sum_chi  = 0;
	sum_max_chi  = 0;
	max_sum_like = 0;
	sum_max_like = 0;

	for (int i = 0; i < n; ++i) {
		int pi = idx_perm[i];

		// These will probably be 1, but not necessarily due to the possibility of ties
		int ranked_dx_i  = ranked_dx[i  * n + i ];
		int ranked_dy_pi = ranked_dy[pi * n + pi];

		// To make life more modular: first create simple arrays with the relevant n-1 distances in x and in y
		for (int j = 0, k = 0; j < n; ++j) {
			sj = (*sorted_dx)[i][j].second;
			if (sj != i) {
				psj = idx_perm[sj]; // NOTE: psj will never be equal to pi
				ranked_dx_sj  = ranked_dx[sj  * n + i ];
				ranked_dy_psj = ranked_dy[psj * n + pi];
				uvs_x[k] = (*sorted_dx)[i][j].first;
				uvs_y[k] = dy[psj * n + pi];
				uvs_xr[k] = ranked_dx_sj  - (ranked_dx_sj  > ranked_dx_i );
				uvs_yr[k] = ranked_dy_psj - (ranked_dy_psj > ranked_dy_pi);
				++k;
			}
		}

		// Then compute the univariate score and aggregate
		(this->*hhg_extended_uvs)();

		sum_chi      += uvs_sc;
		max_chi      = max(max_chi, uvs_mc);
		sum_like     += uvs_sl;
		max_like     = max(max_like, uvs_ml);
		sum_max_chi  += uvs_mc;
		sum_max_like +=	uvs_ml;
		max_sum_chi  = max(max_sum_chi, uvs_sc);
		max_sum_like = max(max_sum_like, uvs_sl);
	}

	sum_chi      /= n;
	sum_like     /= n;
	sum_max_chi  /= n;
	sum_max_like /=	n;
}

// Existing GOF tests
// (used in simulations as a reference to compare the to performance of our tests)
void StatsComputer::uv_gof_existing(void) {
	// See general notes under hhg_gof_ddp() and uv_ks_cvm_ks().

	int n = xy_nrow;

	// will be abused for storing the statistics:
	sum_chi  = 0; // Cramer-von Mises (1928)
	max_chi  = 0; // Kolmogorov-Smirnov (1933)
	sum_like = 0; // LR-based Cramer-von Mises (Zhang 2006)
	max_like = 0; // LR-based Kolmogorov-Smirnov (Zhang 2006)

	// Could this be done faster?
	// Can collect individual tables if wanted

	double Eki, dt, d, lr;

	for (int i = 1; i < n; ++i) {
		Eki = n * null_dist[i];
		dt = i - Eki;
		d = dt * dt / n;
		lr = (i != 0 && i != n) ? (i * log(i / Eki) + (n - i) * log((n - i) / (n - Eki))) : 0;

		sum_chi += d;
		max_chi = max(max_chi, d); // NOTE: sometimes the K-S statistic is defined with an L1 norm and sometimes with L2
		sum_like += lr;
		max_like = max(max_like, lr);
	}
}

void StatsComputer::mv_ts_existing(void) {
	// This representation of the data simplifies the computations of edist and HT
	int j0 = 0, j1 = 0;
	for (int i = 0; i < xy_nrow; ++i) {
		if (y_perm[i] == 0) {
			y0_idx[j0++] = i;
		} else {
			y1_idx[j1++] = i;
		}
	}

	// Compute edist (a version of distance correlation for the 2-sample problem)
	sum_chi = compute_edist();

	// Compute Hall & Tajvidi (2002) test
	sum_like = compute_ht();
}

void StatsComputer::mv_ks_existing(void) {
	// TODO
}

void StatsComputer::mv_ind_existing(void) {
	// TODO
}

// Conditional independence test using the nearest neighbor kernel (univariate z)
void StatsComputer::ci_uvz_nn(void) {
	int n = xy_nrow;
	int nnhh = nnh / 2;
	int nnht = nnhh * 2 - 1;
	double nrmlz = 1.0 / nnht;
	int count[2][2];

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	double dx_ij, dx_ik, dy_ij, dy_ik;

	for (int i = nnhh; i < n - nnhh; ++i) { // FIXMECI what about the edges? completely ignored??
		// Can the following be computed faster? if nnhh is expected to be a small number
		// it doesn't seem like a good idea to compute in a fancy way like hhg_general(). Maybe
		// better to exploit the symmetry in the distance comparison matrices.

		int pi_x = idx_perm_inv[i];
		int pi_y = idx_perm    [i];

		for (int j = i - nnhh; j <= i + nnhh; ++j) {
			// A smoothed bootstrap that resamples both x and y
			if (j != i) {
				count[0][0] = 0;
				count[0][1] = 0;
				count[1][0] = 0;
				count[1][1] = 0;

				int pj_x = idx_perm_inv[j];
				int pj_y = idx_perm    [j];

				dx_ij = dx[pj_x * n + pi_x];
				dy_ij = dy[pj_y * n + pi_y];

				for (int k = i - nnhh; k <= i + nnhh; ++k) {
					if ((k != j) && (k != i)) { // can optimize: this condition can be avoided by unrolling to three loops
						int pk_x = idx_perm_inv[k];
						int pk_y = idx_perm    [k];

						dx_ik = dx[pk_x * n + pi_x];
						dy_ik = dy[pk_y * n + pi_y];

						++count[dx_ij > dx_ik][dy_ij > dy_ik];
					}
				}
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}
}

// Conditional independence test using a gaussian kernel (univariate z)
void StatsComputer::ci_uvz_gaussian(void) {
	int n = xy_nrow;
	double ebw = 3 * sig;
	double nrmlz = 1 / (2 * M_PI * sig * sig), enrmlz = -0.5 / (sig * sig), wijk;
	double count[2][2];

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	Rprintf("NOTE: THIS IS BROKEN\n");

	double dx_ij, dx_ik, dy_ij, dy_ik;
	int zli = 0, zhi = 0, ii, knn;

	for (int i = 0; i < n; ++i) {
		// Can do binary search, but if ebw is not too big that may be useless
		while ((zli < i) && (z[zli] < z[i] - ebw)) {
			++zli;
		}
		while ((zhi < n - 1) && (z[zhi] < z[i] + ebw)) {
			++zhi;
		}
		if (z[zhi] > z[i] + ebw) {
			--zhi;
		}
		knn = zhi - zli;

		for (ii = 0; ii < i - zli; ++ii) {
			idx_1_to_n[ii] = idx_perm[ii] = zli + ii;
		}
		for (; ii < knn; ++ii) {
			idx_1_to_n[ii] = idx_perm[ii] = zli + ii + 1;
		}
		if (should_randomize) {
			// FIXMECI smoothed bootstrap using a separate kernel
			for (ii = knn - 1; ii > 0; --ii) {
				int jj = my_R_rand_wrapper() % (ii + 1); // this can be heavily biased since RAND_MAX is 32K and taking the modulo will give higher probability to smaller numbers //DEBUG_RAND

				int temp = idx_perm[jj];
				idx_perm[jj] = idx_perm[ii];
				idx_perm[ii] = temp;
			}
		}

		for (int j = 0; j < knn; ++j) {
			count[0][0] = 0;
			count[0][1] = 0;
			count[1][0] = 0;
			count[1][1] = 0;

			int jj = idx_1_to_n[j];
			dx_ij = dx[jj * n + i];
			dy_ij = dy[idx_perm[j] * n + i];

			for (int k = 0; k < knn; ++k) {
				if (k != j) {
					int kk = idx_1_to_n[k];
					dx_ik = dx[kk * n + i];
					dy_ik = dy[idx_perm[k] * n + i];

					wijk = nrmlz * exp(enrmlz * ((z[jj] - z[i]) * (z[jj] - z[i]) + (z[kk] - z[i]) * (z[kk] - z[i])));
					count[dx_ij > dx_ik][dy_ij > dy_ik] += wijk;
				}
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}

	should_randomize = false;
}

// Conditional independence test using the nearest neighbor kernel (multivariate z)
void StatsComputer::ci_mvz_nn(void) {
	int n = xy_nrow;
	double nrmlz = 1.0 / (nnh - 1);
	double dx_ij, dx_ik, dy_ij, dy_ik;
	double count[2][2];

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	// NOTE: here I am using nnh nearest neighbors (not necessarily nnh/2 on each "side")
	// because it is not clear to me what is an equivalent notion of a "side" high dimension

	for (int i = 0; i < n; ++i) {
		int pi_x = idx_perm_inv[i];
		int pi_y = idx_perm    [i];

		for (int j = 1; j <= nnh; ++j) {
			count[0][0] = 0;
			count[0][1] = 0;
			count[1][0] = 0;
			count[1][1] = 0;

			// We want the j'th z-neighbor of the i'th sample, in the bootstrapped data.
			int jj = (*sorted_dz)[i][j].second; // the z's are not bootstrapped
			int pj_x = idx_perm_inv[jj]; // and this gives us the original sample that is the j'th neighbor in the bootstrapped data
			int pj_y = idx_perm    [jj];

			dx_ij = dx[pj_x * n + pi_x]; // FIXMECI though this is more readable, working on the transposed matrix may have better locality or reference
			dy_ij = dy[pj_y * n + pi_y];

			for (int k = 1; k < j; ++k) {
				int kk = (*sorted_dz)[i][k].second;
				int pk_x = idx_perm_inv[kk];
				int pk_y = idx_perm    [kk];

				dx_ik = dx[pk_x * n + pi_x];
				dy_ik = dy[pk_y * n + pi_y];

				++count[dx_ij > dx_ik][dy_ij > dy_ik];
			}
			for (int k = j + 1; k <= nnh; ++k) {
				int kk = (*sorted_dz)[i][k].second;
				int pk_x = idx_perm_inv[kk];
				int pk_y = idx_perm    [kk];

				dx_ik = dx[pk_x * n + pi_x];
				dy_ik = dy[pk_y * n + pi_y];

				++count[dx_ij > dx_ik][dy_ij > dy_ik];
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}

	// FIXMECI it may be necessary to normalize the sum statistics only by the number of tables
	// passing the w_sum/w_max conditions
	sum_chi /= nnh;
	sum_like /= nnh;
}

// Conditional independence test using a gaussian kernel (multivariate z)
void StatsComputer::ci_mvz_gaussian(void) {
	int n = xy_nrow;
	double count[2][2];
	double dx_ij, dx_ik, dy_ij, dy_ik;
	int knn;

	// FIXMECI 1: I need to compensate for the truncation I am doing (probably best to
	// precompute the normalizing constant in R and send it here as another kernel
	// parameter)

	// FIXMECI 2: allow specification of a general Gaussian kernel? (with a general
	// covariance matrix). I don't see justification for this.

	double ebw = 3 * sig;
	double nrmlz = 1 / (2 * M_PI * sig * sig), enrmlz = -0.5 / (sig * sig), wijk;

	Rprintf("NOTE: THIS IS BROKEN\n");

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	for (int i = 0; i < n; ++i) {
		// Can do binary search, but if ebw is not too big that may be useless
		for (knn = 0; (knn < n - 1) && ((*sorted_dz)[i][knn + 1].first < ebw); ++knn) {
			idx_1_to_n[knn] = idx_perm[knn] = (*sorted_dz)[i][knn + 1].second;
		}

		if (should_randomize) {
			// FIXMECI smoothed bootstrap using a separate kernel
			for (int ii = knn - 1; ii > 0; --ii) {
				int jj = my_R_rand_wrapper() % (ii + 1); // should use my_rand() implemented below //DEBUG_RAND

				int temp = idx_perm[jj];
				idx_perm[jj] = idx_perm[ii];
				idx_perm[ii] = temp;
			}
		}

		for (int j = 0; j < knn; ++j) {
			count[0][0] = 0;
			count[0][1] = 0;
			count[1][0] = 0;
			count[1][1] = 0;

			int jj = idx_1_to_n[j];
			dx_ij = dx[         jj * n + i];
			dy_ij = dy[idx_perm[j] * n + i];

			for (int k = 0; k < knn; ++k) {
				if (k != j) {
					int kk = idx_1_to_n[k];
					dx_ik = dx[         kk * n + i];
					dy_ik = dy[idx_perm[k] * n + i];

					wijk = nrmlz * exp(enrmlz * ((z[jj] - z[i]) * (z[jj] - z[i]) + (z[kk] - z[i]) * (z[kk] - z[i])));
					count[dx_ij > dx_ik][dy_ij > dy_ik] += wijk;
				}
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}

	should_randomize = false;
}

// Conditional independence test using the nearest neighbor kernel (multivariate z)
// In this version we go over a grid of bandwidth values and compute the regular
// test for each of them. The main set of statistics is then replaced with the
// maximum statistic across the grid.
void StatsComputer::ci_mvz_nn_grid(void) {
	double max_nnh_sum_chi  = 0;
	double max_nnh_max_chi  = 0;
	double max_nnh_sum_like = 0;
	double max_nnh_max_like = 0;

	for (int i = 0; i < nnh_grid_cnt; ++i) {
		nnh = nnh_grid[i];
		ci_mvz_nn();
		sum_chi_grid [i] = sum_chi;
		sum_like_grid[i] = sum_like;
		max_chi_grid [i] = max_chi;
		max_like_grid[i] = max_like;

		if (sum_chi > max_nnh_sum_chi) {
			max_nnh_sum_chi = sum_chi;
		}
		if (sum_like > max_nnh_sum_like) {
			max_nnh_sum_like = sum_like;
		}
		if (max_chi > max_nnh_max_chi) {
			max_nnh_max_chi = max_chi;
		}
		if (max_like > max_nnh_max_like) {
			max_nnh_max_like = max_like;
		}
	}

	sum_chi  = max_nnh_sum_chi;
	max_chi  = max_nnh_max_chi;
	sum_like = max_nnh_sum_like;
	max_like = max_nnh_max_like;
}

// A variant of the ADP statistic used for CI testing
// NOTE: this is currently NOT a distribution-free test!
void StatsComputer::ci_udf_adp_mvz_nn(void) {
	int n = xy_nrow;
	int rect_o, xh, yh;
	double cnt, rect_e, rect_c, rect_l;
	double edenom = 1.0 / nnh;

	sum_chi  = 0;
	sum_like = 0;
	max_chi  = 0; // NOTE:
	max_like = 0; // these will not be computed (the fast algorithm for ADP does not support max)

	for (int i = 0; i < n; ++i) {
		// 1. Start by computing the double integral of paired-sample indicators
		// NOTE: see further notes under compute_double_integral
		memset(double_integral, 0, sizeof(int) * dintegral_pn * dintegral_pn);

		// 1.1 Compute x- and y-ranks of the neighbors
		for (int j = 0; j < nnh; ++j) {
			// We want the j'th z-neighbor of the i'th sample, in the bootstrapped data.
			int jj = (*sorted_dz)[i][j + 1].second; // the z's are not bootstrapped
			int pj_x = idx_perm_inv[jj]; // and this gives us the original sample that is the j'th neighbor in the bootstrapped data
			int pj_y = idx_perm    [jj];
			nn_sorted_x[j].first = dx[pj_x];
			nn_sorted_y[j].first = dy[pj_y];
			nn_sorted_x[j].second = j;
			nn_sorted_y[j].second = j;
		}
		sort(nn_sorted_x.begin(), nn_sorted_x.end(), dbl_int_pair_comparator);
		sort(nn_sorted_y.begin(), nn_sorted_y.end(), dbl_int_pair_comparator);
		for (int j = 0; j < nnh; ++j) {
			nn_sorted_x[nn_sorted_x[j].second].first = j + 1;
			nn_sorted_y[nn_sorted_y[j].second].first = j + 1;
		}

		// 1.2. Populate the padded matrix with the indicator variables of whether there
		// is a point in the neighborhood
		for (int j = 0; j < nnh; ++j) {
			int rxj = nn_sorted_x[j].first;
			int ryj = nn_sorted_y[j].first;
			double_integral[ryj * dintegral_pn + rxj] = 1;
		}

		// 1.3. Then run linearly and compute the integral in one row-major pass
		int la = dintegral_pn;
		for (int j = 1; j < dintegral_pn; ++j) {
			int row_running_sum = 0;
			++la;
			for (int k = 1; k < dintegral_pn; ++k) {
				row_running_sum += double_integral[la];
				double_integral[la] = row_running_sum + double_integral[la - dintegral_pn];
				++la;
			}
		}

		// 2. Now use this to compute the ADP statistic for the neighborhood
		// NOTE: see further notes under hhg_udf_adp
		double nn_sum_chi = 0, nn_sum_like = 0;
		double nr_nonempty_cells = 0;
		kahan_c_chi = 0;
		kahan_c_like = 0;
		double kahan_t;

		for (int w = 1; w <= nnh; ++w) {
			for (int h = 1; h <= nnh; ++h) {
				for (int xl = 1; xl <= nnh - w + 1; ++xl) {
					for (int yl = 1; yl <= nnh - h + 1; ++yl) {
						xh = xl + w - 1;
						yh = yl + h - 1;
						cnt = count_adp_with_given_cell(xl, xh, yl, yh);
						if (cnt > 0) {
							rect_o = count_sample_points_in_rect(xl, xh, yl, yh);
							rect_e = w * h * edenom;
							rect_c = ((rect_o - rect_e) * (rect_o - rect_e)) / rect_e * cnt - kahan_c_chi;
							rect_l = ((rect_o > 0) ? (rect_o * log(rect_o / rect_e)) : 0) * cnt - kahan_c_like;

							kahan_t = nn_sum_chi + rect_c;
							kahan_c_chi = (kahan_t - nn_sum_chi) - rect_c;
							nn_sum_chi = kahan_t;

							kahan_t = nn_sum_like + rect_l;
							kahan_c_like = (kahan_t - nn_sum_like) - rect_l;
							nn_sum_like = kahan_t;

							if (rect_o > 0) {
								nr_nonempty_cells += cnt;
							}
						}
					}
				}
			}
		}

#ifdef XDP_NORMALIZE
		double nr_parts = choose(nnh - 1, K - 1); // this is actually only the sqrt of the nr parts
		nr_parts *= nr_parts;

		if (correct_mi_bias) {
			double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
			nn_sum_chi += mm_bias;
			nn_sum_like += mm_bias;
		}

		double normalizer = nr_parts * nnh; // gets us to the MI scale
		nn_sum_chi /= normalizer;
		nn_sum_like /= normalizer;
#endif

		sum_chi += nn_sum_chi;
		sum_like += nn_sum_like;
	}

#ifdef XDP_NORMALIZE
	sum_chi /= n;
	sum_like /= n;
#endif
}

// Univariate goodness-of-fit scores
// ================================================================================================

void StatsComputer::uvs_gof_wxn(void) {
	// TODO
}

void StatsComputer::uvs_gof_ad(void) {
	// TODO
}

void StatsComputer::uvs_gof_cvm_ks(void) {
	// TODO
}

void StatsComputer::uvs_gof_dcov(void) {
	// TODO
}

void StatsComputer::uvs_gof_xdp2(void) {
	// See general notes under hhg_gof_ddp()

	int n = uvs_n;

	uvs_sc  = 0;
	uvs_mc  = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Could this be done faster?
	// Can collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// K here is 2

	double chi, like, emin;

	for (int xi = 1; xi < n; ++xi) {
		tbl_o[0] = xi;
		tbl_o[1] = n - xi;
		tbl_e[0] = n * null_dist[xi];
		tbl_e[1] = n * (1 - null_dist[xi]);

		chi = (  (tbl_o[0] - tbl_e[0]) * (tbl_o[0] - tbl_e[0]) / tbl_e[0]
			   + (tbl_o[1] - tbl_e[1]) * (tbl_o[1] - tbl_e[1]) / tbl_e[1]);

		like = tbl_o[0] * log(tbl_o[0] / tbl_e[0]) + tbl_o[1] * log(tbl_o[1] / tbl_e[1]);

		emin = min(tbl_e[0], tbl_e[1]);

		accumulate_local_stats(chi, like, emin);
	}

	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
}

void StatsComputer::uvs_gof_xdp3(void) {
	// See general notes under hhg_gof_ddp()

		int n = uvs_n;

		uvs_sc  = 0;
		uvs_mc  = 0;
		uvs_sl = 0;
		uvs_ml = 0;

		ng_chi  = 0;
		ng_like = 0;

		// Could this be done in O(nlogn) time?
		// Can collect individual tables if wanted

		// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
		// K here is 3

		double chi, like, etmp, emin;

		for (int xi = 1; xi < n - 1; ++xi) {
			for (int xj = xi + 1; xj < n; ++xj) {
				tbl_o[0] = xi;
				tbl_o[1] = xj - xi;
				tbl_o[2] = n - xj;
				tbl_e[0] = n * null_dist[xi];
				tbl_e[1] = n * (null_dist[xj] - null_dist[xi]);
				tbl_e[2] = n * (1 - null_dist[xj]);

				chi = (  (tbl_o[0] - tbl_e[0]) * (tbl_o[0] - tbl_e[0]) / tbl_e[0]
					   + (tbl_o[1] - tbl_e[1]) * (tbl_o[1] - tbl_e[1]) / tbl_e[1]
					   + (tbl_o[2] - tbl_e[2]) * (tbl_o[2] - tbl_e[2]) / tbl_e[2]);

				like = 0;
				if (tbl_o[0] > 0) {
					like += tbl_o[0] * log(tbl_o[0] / tbl_e[0]);
				}
				if (tbl_o[1] > 0) {
					like += tbl_o[1] * log(tbl_o[1] / tbl_e[1]);
				}
				if (tbl_o[2] > 0) {
					like += tbl_o[2] * log(tbl_o[2] / tbl_e[2]);
				}

				etmp = min(tbl_e[0], tbl_e[1]);
				emin = min(etmp, tbl_e[2]);

				accumulate_local_stats(chi, like, emin);
			}
		}

		uvs_sc /= (double(n) * ng_chi);
		uvs_sl /= (double(n) * ng_like);
}

void StatsComputer::uvs_gof_xdp(void) {
	// It is assumed that y contains the null CDF, computed at the midpoints
	// between every two consecutive observations (so there are n-1 given
	// values between 0 and 1)

	int n = uvs_n;

	uvs_sc  = 0;
	uvs_sl = 0;

	// These can't be computed this way
	uvs_mc  = 0;
	uvs_ml = 0;

	// Can collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// For GOF testing, we could opt to partition at any point in the interval
	// between every two consecutive observations, and this choice actually matters.
	// The choice of the midpoints is arbitrary, but will do for now.

	// This algorithm works interval-wise rather than partition-wise.

	double normalized_cnt;
	double interval_o, interval_e, interval_c, interval_l;
	kahan_c_chi = 0;
	kahan_c_like = 0;
	double kahan_t;
	int xj, wmax;

	for (int xi = 0; xi < n; ++xi) {
		wmax = min(n - K - 1, n - xi);

		for (int w = 1; w <= wmax; ++w) {
			xj = xi + w;

			// could be slightly optimized by going over edge intervals in a separate loop from mid intervals
			if (xi == 0 || xj == n) {
				// edge interval
				normalized_cnt = adp_l[w];
			} else {
				// mid interval
				normalized_cnt = adp[w];
			}

			interval_o = xj - xi;
			interval_e = (((xj == n) ? 1.0 : null_dist[xj]) - null_dist[xi]) * n;

			interval_c = ((interval_o - interval_e) * (interval_o - interval_e)) / interval_e * normalized_cnt - kahan_c_chi;
			interval_l = ((interval_o > 0) ? (interval_o * log(interval_o / interval_e)) : 0) * normalized_cnt - kahan_c_like;

			kahan_t = uvs_sc + interval_c;
			kahan_c_chi = (kahan_t - uvs_sc) - interval_c;
			uvs_sc = kahan_t;

			kahan_t = uvs_sl + interval_l;
			kahan_c_like = (kahan_t - uvs_sl) - interval_l;
			uvs_sl = kahan_t;
		}
	}

	// the 1/n gets us to the MI scale
	uvs_sc /= n;
	uvs_sl /= n;
}

// Univariate two- and k-sample scores
// ================================================================================================

void StatsComputer::uvs_ks_kw(void) {
	int n = uvs_n, k, i;
	double xr_avg = 0.5 * (n + 1), tmp, numer = 0, denom = 0;

	for (k = 0; k < nr_groups; ++k) {
		kw_rs[k] = 0;
	}

	for (i = 0; i < n; ++i) {
		k = uvs_yr[i];
		kw_rs[k] += uvs_xr[i];
		tmp = uvs_xr[i] - xr_avg;
		denom += tmp * tmp;
	}

	for (k = 0; k < nr_groups; ++k) {
		tmp = kw_rs[k] / uvs_yc[k] - xr_avg;
		numer += uvs_yc[k] * tmp * tmp;
	}

	uvs_sc = numer / denom;
}

void StatsComputer::uvs_ks_ad(void) {
	// FIXME for the moment I am implementing the version that assumes there are no tied
	// observations (as appears in Scholz, 1986, Eq. 3). They have a more general
	// statistic in Eq. 4 and 5 to support ties.

	compute_single_integral(uvs_n, uvs_xr, uvs_yr);

	int n = uvs_n;
	double tmp, tmp_y;

	uvs_sc = 0; // will be abused for storing the statistic

	for (int yi = 0; yi < nr_groups; ++yi) {
		tmp_y = 0;
		for (int xi = 1; xi < n; ++xi) {
			tmp = n * double_integral[yi * dintegral_pn + xi] - xi * uvs_yc[yi];
			tmp_y += tmp * tmp / (xi * (n - xi));
		}
		uvs_sc += tmp_y / uvs_yc[yi];
	}

	uvs_sc /= n;
}

void StatsComputer::uvs_ks_cvm_ks(void) {
	// NOTE: The tests implemented here also exist in versions that use at-ranks grid for
	// inducing dichotomies. This creates the question of how to estimate the
	// empirical CDF at the partitioning points, and some people split each
	// rank and count it as being left on each side of the partition. I avoid
	// this by performing partitions at half-ranks as in ADP and our K-sample
	// test. The effect of this is expected to be negligible.
	// Also, I do not use w_sum/w_max. Maybe that's not fair but the intention
	// is to compare to the "plain vanilla" classical tests.

	// FIXME: need to think about ties: is it fine as it is?

	int n = uvs_n;
	compute_single_integral(n, uvs_xr, uvs_yr);

	// will be kind of abused for storing the statistics:
	uvs_sc = 0; // Cramer-von Mises (Kiefer, 1959)
	uvs_mc = 0; // Kolmogorov-Smirnov (Kiefer, 1959)
	uvs_sl = 0; // LR-based Cramer-von Mises (Zhang 2007)
	uvs_ml = 0; // LR-based Kolmogorov-Smirnov (Zhang 2007)

	// Could this be done faster?
	// May need to address roundoff issues
	// Can collect individual tables if wanted

	// NOTE: actually, this is very confusing, since Zhang doesn't simply switch from Pearson chi-squared
	// to the G statistic, rather he also changes the weights when he defines his KS, AD, and CvM

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// K here is 2

	double Nki, Nk, Eki, dt, d, lr, nrcp = 1.0 / n;

	for (int i = 1; i < n; ++i) {
		d = 0;
		lr = 0;
		for (int k = 0; k < nr_groups; ++k) {
			Nk = uvs_yc[k];
			Eki = Nk * nrcp * i;
			Nki = double_integral[k * dintegral_pn + i];
			dt = Nki - Eki;
			d += dt * dt / Nk;
			lr += (Nki != 0 && Nki != Nk) ? (Nki * log(Nki / Eki) + (Nk - Nki) * log((Nk - Nki) / (Nk - Eki))) : 0;
		}

		uvs_sc += d;
		uvs_mc = max(uvs_mc, d); // NOTE: sometimes the K-S statistic is defined with an L1 norm and sometimes with L2
		uvs_sl += lr;
		uvs_ml = max(uvs_ml, lr);
	}
}

void StatsComputer::uvs_ks_dcov(void) {
	int n = uvs_n;

	double d0 = 0, d1 = 0, ds = 0;
	int nr0 = 0, nr1 = 0;

	// NOTE: this is a naive implementation, which may not be very numerically stable.

	for (int i = 0; i < n; ++i) {
		if (uvs_yr[i] == uvs_y0) {
			d0 += uvs_x[i];
			++nr0;
		} else {
			d1 += uvs_x[i];
			++nr1;
		}

		ds += uvs_x[i] * uvs_x[i];
	}

	double d = d1 / nr1 - d0 / nr0;
	double mean_d = (d0 + d1) / n;
	double se_d = sqrt(ds / n - mean_d * mean_d);

	uvs_sc  = d; 		// the original dCov statistic
	uvs_mc  = fabs(d); 	// our "abs" variant
	uvs_sl  = d / se_d; // out "t test" variant
	uvs_ml  = 0; 		// not used yet
}

// NOTE: Dynamic Slicing code was originally written by Chao Ye and Bo Jiang. I've made some
// minimal modifications to fit it into my framework.
void StatsComputer::uvs_ks_ds(void) {
	int* x = uvs_yr;
	int len = uvs_n;
	int dim = nr_groups;

	double lpd = -lambda * log((double)(len));  //  penalty (log of prior odds) for each additional slice

	for (int k = 0; k < len + 1; ++k) {
		for (int j = 0; j < dim; ++j) {
			ds_ctab[k][j] = 0;
		}
	}

	int flagl = 1;
	int clpcount = 1;  //  clump count
	int clumpnum = 1;  //  clump number

	while (flagl < len) {
		if (x[flagl] - x[flagl - 1] != 0) {
			ds_ctab[clumpnum][x[flagl - 1]] = clpcount;
			clpcount = 1;
			clumpnum++;
		} else {
			clpcount++;
		}

		flagl++;
	}

	ds_ctab[clumpnum][x[len - 1]] = clpcount;

	for (int k = 1; k < clumpnum + 1; ++k) {
		for (int j = 0; j < dim; ++j) {
			ds_ctab[k][j] += ds_ctab[k - 1][j];
		}
	}

	//  dynamic programming
	for (int k = 0; k < clumpnum + 1; ++k) {
		ds_score[k] = 0;
		ds_idx[k] = -1;
	}

	int tc, cutpos;
	double tpcut, cutsc;

	for (int i = 1; i < clumpnum + 1; ++i) {
		//  j = 0
		cutsc = lpd + ds_score[0];
		tc = 0;

		for (int u = 0; u < dim; ++u) {
			tc += ds_ctab[i][u];
		}

		for (int u = 0; u < dim; ++u) {
			ds_counts[u] = (double)(ds_ctab[i][u]);
			if (ds_counts[u] > 1e-6) {
				cutsc += ds_counts[u] * log(ds_counts[u] / tc);
			}
		}

		cutpos = 0; //  j = 0 end
		for (int j = 1; j < i; ++j) {
			tpcut = lpd + ds_score[j];
			tc = 0;
			for (int u = 0; u < dim; ++u) {
				tc += ds_ctab[i][u] - ds_ctab[j][u];
			}
			for (int u = 0; u < dim; ++u) {
				ds_counts[u] = (double)(ds_ctab[i][u] - ds_ctab[j][u]);
				if (ds_counts[u] > 1e-6) {
					tpcut += ds_counts[u] * log(ds_counts[u] / tc);
				}
			}
			if (cutsc < tpcut) {
				cutsc = tpcut;
				cutpos = j;
			}
		}
		ds_score[i] = cutsc;
		ds_idx[i] = cutpos;
	}

	double mlik = ds_score[clumpnum] - lpd; // subtract null log-likelihood (assume one slice, i.e., no cut)
	for (int u = 0; u < dim; ++u) {
		if (ds_ctab[clumpnum][u] > 1e-6) {
			mlik -= ds_ctab[clumpnum][u] * log((double)(ds_ctab[clumpnum][u]) / len);
		}
	}

	uvs_sc  = mlik;	// maximized log-likelihood
	uvs_mc  = 0; 	// not used yet
	uvs_sl  = 0; 	// not used yet
	uvs_ml  = 0; 	// not used yet
}


void StatsComputer::uvs_ks_mds(void) {
	int* x = uvs_yr;
	int len = uvs_n;
	int dim = nr_groups;
	
	
	double *group_ratios = new double[dim];

	//
	// Figure out the clumps
	//

	// This is a table with rows representing the clumps of consecutive x
	// values (in the sorted sequence) that come from the same original
	// sample. First each row will have zeros in all places, except for one column,
	// whose number is the value of the original sample that all observations
	// in the clump came from. Then they compute a cumulative sum in columns.
	// So in the end the value at k,j is the number of samples of y==j that are
	// part of clump k or lower.
	
	int clpcount = 1;
	int clumpnum = 1;
	
		
	//used only in equipartition_type == 2
	int wanted_nr_bins = equipartition_nr_cells_m;//(int)(sqrt((double)(len)));
	
	if(equipartition_type == 0 || equipartition_type == 1 || equipartition_type == 2){
		//the regular clumping mechanism
		for (int k = 0; k < len + 1; ++k) {
			for (int j = 0; j < dim; ++j) {
				ds_ctab[k][j] = 0;
			}
		}
		
		
		for (int i = 0; i < len - 1; ++i) {
			if (x[i + 1] != x[i] || equipartition_type == 1 || equipartition_type == 2){
				ds_ctab[clumpnum][x[i]] = clpcount;
				clpcount = 1;
				++clumpnum;
			} else {
				++clpcount;
			}
		}

		ds_ctab[clumpnum][x[len - 1]] = clpcount;

		for (int k = 1; k < clumpnum + 1; ++k) {
			for (int j = 0; j < dim; ++j) {
				ds_ctab[k][j] += ds_ctab[k - 1][j];
			}
		}
	}
	if(equipartition_type == 2 ){
		int pointer = 0;
		int clump_counter=0;
		for(int i=1;i<=wanted_nr_bins;i++){
			if(i==wanted_nr_bins){
				pointer = len;
			}else{
				pointer = (int)(i * len /wanted_nr_bins);
			}
			for (int j = 0; j < dim; ++j){
				ds_ctab_bins[clump_counter + 1][j] = ds_ctab[pointer][j];
			}
			clump_counter ++;
		}
		//last row is the summation of the nr obs per group
		for (int j = 0; j < dim; ++j) { 
				ds_ctab_bins[clump_counter + 1][j] = ds_ctab[len][j];
		}
		
		
		clumpnum = clump_counter;
		for (int k = 1; k < clumpnum + 1; ++k){
			for (int j = 0; j < dim; ++j){
				ds_ctab[k][j] = ds_ctab_bins[k][j];
			}
		}
	}
	
	
	
	for(int j = 0; j < dim; ++j){
		group_ratios[j] = ((double)ds_ctab[clumpnum][j])/((double)len);
	}
	
	//
	// Dynamic programming
	//

	int cc = clumpnum + 1;
	for (int i = 0; i < cc; ++i) {
		for (int a = 0; a <= i; ++a) {
			ds_score[i * cc + a] = 0;
			ds_score_pearson[i * cc + a] = 0;
		}
	}

	int tc, tmpi;
	double tpcut, cutsc, tmpd;
	
	double tpcut_pearson,cutsc_pearson,temp_pearson;
	for (int i = 1; i < cc; ++i) {
		// a = 1 (only j = 0 is possible)
		tc = 0;
		for (int u = 0; u < dim; ++u) {
			tc += ds_ctab[i][u];
		}
		cutsc = 0;
		cutsc_pearson = 0;
		for (int u = 0; u < dim; ++u) {
			tmpd = ds_ctab[i][u];
			if (ds_ctab[i][u] != 0) {
				cutsc += tmpd * log(tmpd / tc);
			}
			temp_pearson = ((double)tc)*group_ratios[u]; //this is the expected score;
			cutsc_pearson += (tmpd - temp_pearson)*(tmpd -temp_pearson)/temp_pearson;
		}
		ds_score[i * cc + 1] = cutsc;
		ds_score_pearson[i * cc + 1] = cutsc_pearson; 

		for (int a = 2; a <= i; ++a) {
			cutsc = -DBL_MAX;
			cutsc_pearson= -DBL_MAX; 
			for (int j = a - 1; j < i; ++j) {
				tc = 0;
				for (int u = 0; u < dim; ++u) {
					tc += ds_ctab[i][u] - ds_ctab[j][u];
				}
				tpcut = ds_score[j * cc + a - 1];
				tpcut_pearson = ds_score_pearson[j * cc + a - 1];
				for (int u = 0; u < dim; ++u) {
					tmpi = ds_ctab[i][u] - ds_ctab[j][u];
					tmpd = tmpi;
					if (tmpi != 0) {
						tpcut += tmpd * log(tmpd / tc);
					}
					temp_pearson = ((double)tc)*group_ratios[u]; //this is the expected score; 
					tpcut_pearson += (tmpd -temp_pearson)*(tmpd -temp_pearson)/temp_pearson; 
				}

				if (cutsc < tpcut) {
					cutsc = tpcut;
				}
				if(cutsc_pearson <tpcut_pearson){ 
					cutsc_pearson = tpcut_pearson;
				}
			}

			ds_score[i * cc + a] = cutsc;
			ds_score_pearson[i * cc + a] = cutsc_pearson; 
		}
	}

	
	//computing constant to bring back to log-likelihood
	double ll_const=0;
	for (int u = 0; u < dim; ++u) {
		if (ds_ctab[clumpnum][u] != 0) {
			ll_const -= ds_ctab[clumpnum][u] * log((double) (ds_ctab[clumpnum][u]) / len);
		}
	}

	
	// Now choose the best score among constrained problems, when adding the penalty
#if 1 // with the suggested modified penalty
	double mlik = ds_score[clumpnum * cc + 2]-prior[0];   
	double mchi = ds_score_pearson[clumpnum * cc + 2]-prior[0]; 
	int maxa = min(cc-1,Mk_Maxk ); 
	int chosen_k=2;
	int chosen_k_pearson=2;
	if(prior_length<= maxa-2){
	mchi=-999.0;
	mlik=-999.0;
	}
	else{
	for (int a = 2; a <= maxa; ++a) {
		
		tmpd = ds_score[clumpnum * cc + a] - prior[a-2];
		
		if(mlik<tmpd){
		mlik=tmpd;
		chosen_k=a;
		}
		
		tmpd = ds_score_pearson[clumpnum * cc + a] - prior[a-2] ; 
		
		if(mchi<tmpd){
		mchi=tmpd;
		chosen_k_pearson=a;
		}
		
	}
	}
#else // sanity check: this should be identical to the original DS
	double mlik = ds_score[clumpnum * cc + 1];
	int chosen_k=1;
	double mchi = 0; //not supported in this section
	int chosen_k_pearson=1; //not supported in this section

	int maxa = min(cc-1, Mk_Maxk); 
	for (int a = 1; a <= maxa; ++a) {
		tmpd = ds_score[clumpnum * cc + a] - lambda*log(len) * (a - 1);
		if(mlik<tmpd){
		mlik=tmpd;
		chosen_k=a;
		}
	}
#endif

	//handle the k-specific maximum
	int a=2;
	for(;a<=maxa;a++){
		//the indexes are because mds_max_chi_by_k is indexed from zero, with k=2 in the first place. 
		mds_max_loglikelihood_by_k[a-2]=ds_score[clumpnum * cc + a]+ll_const;
		mds_max_chi_by_k[a-2]=ds_score_pearson[clumpnum * cc + a];
	}
	mlik += ll_const;
	
	uvs_sc  = chosen_k_pearson; 	// k- chosen for chi-square
	uvs_mc  = mchi; 	// maximized penalized chi-square
	uvs_sl  = chosen_k; 	// k chosen for maximum
	uvs_ml  = mlik; // maximized log-likelihood
	
	delete[] group_ratios;
	
}

void StatsComputer::uvs_ks_xdp2(void) {
	compute_single_integral(uvs_n, uvs_xr, uvs_yr);

	int n = uvs_n;
	double nrmlz = 1.0 / n;

	uvs_sc  = 0;
	uvs_mc  = 0;
	uvs_sl  = 0;
	uvs_ml  = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Could this be done faster?
	// Can collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// K here is 2

	double chi, like, etmp, emin;
	int yiK;

	for (int xi = 1; xi < n; ++xi) {
		chi = like = 0;
		emin = n;

		for (int yi = 0; yi < nr_groups; ++yi) {
			yiK = yi * K;

			tbl_o[yiK + 0] = double_integral[yi * dintegral_pn + xi];
			tbl_o[yiK + 1] = uvs_yc[yi] - double_integral[yi * dintegral_pn + xi];
			tbl_e[yiK + 0] = uvs_yc[yi] * double_integral[nr_groups * dintegral_pn + xi] * nrmlz;
			tbl_e[yiK + 1] = uvs_yc[yi] * (n - double_integral[nr_groups * dintegral_pn + xi]) * nrmlz;

			chi += (  (tbl_o[yiK    ] - tbl_e[yiK    ]) * (tbl_o[yiK    ] - tbl_e[yiK    ]) / tbl_e[yiK    ]
			        + (tbl_o[yiK + 1] - tbl_e[yiK + 1]) * (tbl_o[yiK + 1] - tbl_e[yiK + 1]) / tbl_e[yiK + 1]);

			if (tbl_o[yiK] > 0) {
				like += tbl_o[yiK] * log(tbl_o[yiK] / tbl_e[yiK]);
			}
			if (tbl_o[yiK + 1] > 0) {
				like += tbl_o[yiK + 1] * log(tbl_o[yiK + 1] / tbl_e[yiK + 1]);
			}

			etmp = min(tbl_e[yiK], tbl_e[yiK + 1]);
			emin = min(emin, etmp);
		}

		accumulate_local_stats(chi, like, emin);
	}

	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
}

void StatsComputer::uvs_ks_xdp3(void) {
	compute_single_integral(uvs_n, uvs_xr, uvs_yr);

	int n = uvs_n;
	double nrmlz = 1.0 / n;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Could this be done in O(nlogn) time?
	// Can collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// K here is 3

	double chi, like, etmp, emin;
	int yiK;

	for (int xi = 1; xi < n - 1; ++xi) {
		for (int xj = xi + 1; xj < n; ++xj) {
			chi = like = 0;
			emin = n;

			for (int yi = 0; yi < nr_groups; ++yi) {
				yiK = yi * K;

				tbl_o[yiK + 0] = double_integral[yi * dintegral_pn + xi];
				tbl_o[yiK + 1] = double_integral[yi * dintegral_pn + xj] - double_integral[yi * dintegral_pn + xi];
				tbl_o[yiK + 2] = uvs_yc[yi] - double_integral[yi * dintegral_pn + xj];
				tbl_e[yiK + 0] = uvs_yc[yi] * double_integral[nr_groups * dintegral_pn + xi] * nrmlz;
				tbl_e[yiK + 1] = uvs_yc[yi] * (double_integral[nr_groups * dintegral_pn + xj] - double_integral[nr_groups * dintegral_pn + xi]) * nrmlz;
				tbl_e[yiK + 2] = uvs_yc[yi] * (n - double_integral[nr_groups  * dintegral_pn + xj]) * nrmlz;

				chi += (  (tbl_o[yiK + 0] - tbl_e[yiK + 0]) * (tbl_o[yiK + 0] - tbl_e[yiK + 0]) / tbl_e[yiK + 0]
						+ (tbl_o[yiK + 1] - tbl_e[yiK + 1]) * (tbl_o[yiK + 1] - tbl_e[yiK + 1]) / tbl_e[yiK + 1]
						+ (tbl_o[yiK + 2] - tbl_e[yiK + 2]) * (tbl_o[yiK + 2] - tbl_e[yiK + 2]) / tbl_e[yiK + 2]);

				if (tbl_o[yiK + 0] > 0) {
					like += tbl_o[yiK + 0] * log(tbl_o[yiK + 0] / tbl_e[yiK + 0]);
				}
				if (tbl_o[yiK + 1] > 0) {
					like += tbl_o[yiK + 1] * log(tbl_o[yiK + 1] / tbl_e[yiK + 1]);
				}
				if (tbl_o[yiK + 2] > 0) {
					like += tbl_o[yiK + 2] * log(tbl_o[yiK + 2] / tbl_e[yiK + 2]);
				}

				etmp = min(tbl_e[yi * 2 + 0], tbl_e[yi * 2 + 1]);
				etmp = min(etmp, tbl_e[yi * 2 + 2]);
				emin = min(emin, etmp);
			}

			accumulate_local_stats(chi, like, emin);
		}
	}

	uvs_sc /= (double(n) * ng_chi);
	uvs_sl /= (double(n) * ng_like);
}

void StatsComputer::uvs_ks_xdp(void) {
	compute_single_integral(uvs_n, uvs_xr, uvs_yr);

	int n = uvs_n;
	double nrmlz = 1.0 / n;

	uvs_sc  = 0;
	uvs_sl = 0;

	// These can't be computed this way
	uvs_mc = 0;
	uvs_ml = 0;

	// Can collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// This algorithm works interval-wise rather than partition-wise.

	double normalized_cnt;
	double interval_o, interval_e, interval_c, interval_l;
	kahan_c_chi = 0;
	kahan_c_like = 0;
	double kahan_t;
	int xj, wmax;

	for (int xi = 0; xi < n; ++xi) {
		wmax = n-xi;

		for (int w = 1; w <= wmax; ++w) {
			xj = xi + w;

			// could be slightly optimized by going over edge intervals in a separate loop from mid intervals
			if (xi == 0 || xj == n) {
				// edge interval
				normalized_cnt = adp_l[w];
			} else {
				// mid interval
				normalized_cnt = adp[w];
			}

			for (int yi = 0; yi < nr_groups; ++yi){
				interval_o =               double_integral[yi        * dintegral_pn + xj] - double_integral[yi        * dintegral_pn + xi];
				interval_e = uvs_yc[yi] * (double_integral[nr_groups * dintegral_pn + xj] - double_integral[nr_groups * dintegral_pn + xi]) * nrmlz;

				interval_c = ((interval_o - interval_e) * (interval_o - interval_e)) / interval_e * normalized_cnt - kahan_c_chi;
				interval_l = ((interval_o > 0) ? (interval_o * log(interval_o / interval_e)) : 0) * normalized_cnt - kahan_c_like;

				kahan_t = uvs_sc + interval_c;
				kahan_c_chi = (kahan_t - uvs_sc) - interval_c;
				uvs_sc = kahan_t;

				kahan_t = uvs_sl + interval_l;
				kahan_c_like = (kahan_t - uvs_sl) - interval_l;
				uvs_sl = kahan_t;
			}
		}
	}

	// the 1/n gets us to the MI scale
	uvs_sc /= n;
	uvs_sl /= n;
}


void StatsComputer::uvs_ks_xdp_mk(void) {
	compute_single_integral(uvs_n, uvs_xr, uvs_yr);
	int perform_kahan = 0;
	int n = uvs_n;
	double nrmlz = 1.0 / n;

	uvs_sc  = 0;
	uvs_sl = 0;

	// These can't be computed this way
	uvs_mc = 0;
	uvs_ml = 0;

	// Can collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// This algorithm works interval-wise rather than partition-wise.

	
	double interval_o, interval_e, interval_c, interval_l,interval_xi_o;
	double *kahan_c_chi_mk, *kahan_c_like_mk,*kahan_c_chi_mk_edge,*kahan_c_like_mk_edge, kahan_t=0,temp=0;
	double *sum_chi_w,*sum_like_w,*sum_chi_w_edge,*sum_like_w_edge;
	int *ecdf_dictionary;
	int *cell_pointers;
	
	int wmax_limit;
	wmax_limit= n; //wmax can't reach more than this number
	
	if(equipartition_type == 0){
		wmax_limit= n; //wmax can't reach more than this number
	}else if(equipartition_type ==1){
		wmax_limit= equipartition_nr_cells_m;
	}
	kahan_c_chi_mk = new double[wmax_limit+1]; // note that the kahan summation in the first part reduces the error from O(n^2) to O(n), for the second part it reduces it from O(n) to O(sqrt n).
	kahan_c_like_mk = new double[wmax_limit+1];
	kahan_c_chi_mk_edge = new double[wmax_limit+1]; 
	kahan_c_like_mk_edge = new double[wmax_limit+1];
	sum_chi_w = new double[wmax_limit+1]; //for holding the sums, by w, of contribution of cells, arrays for edge and non edge cells, for chi square sums and likelihood
	sum_like_w = new double[wmax_limit+1];
	sum_chi_w_edge = new double[wmax_limit+1]; 
	sum_like_w_edge = new double[wmax_limit+1];
	cell_pointers = new int[wmax_limit+1];
	
	ecdf_dictionary = new int[nr_groups*(wmax_limit+1)];
	
	int xj, wmax;
	xj=0;
	for(int k=0;k<K-1;k++){ //the first entry is for k=2
		xdp_sc_mk[k] = 0;
		xdp_sl_mk[k] = 0;
	}
	for(int i=0;i<wmax_limit+1;i++){
		kahan_c_chi_mk[i] = 0;
		kahan_c_like_mk[i] = 0;
		kahan_c_chi_mk_edge[i] = 0;
		kahan_c_like_mk_edge[i] = 0;
		sum_chi_w[i] = 0;
		sum_like_w[i] = 0;
		sum_chi_w_edge[i] = 0;
		sum_like_w_edge[i] = 0;
	}
	
	if(equipartition_type == 0){
		for(int i=0;i<wmax_limit+1;i++){
			cell_pointers[i]=i;
		}
	}else if(equipartition_type ==1){
	
		for(int i=0;i<wmax_limit+1;i++){
			if(i != wmax_limit){
				cell_pointers[i]=(int)(((double) i*n)/((double)wmax_limit)); //since wmax_limit is sqrt(n)
			}else{
				cell_pointers[i]=n;
			}
		}
	}
	for(int yi=0;yi< nr_groups;++yi){
		for(int i=0;i<wmax_limit+1;i++){
			ecdf_dictionary[yi*(wmax_limit+1)+i] = double_integral[yi        * dintegral_pn + cell_pointers[i]] ;
		}
	}
	
	int pxi;
	int pxj;
	
	for (int yi = 0; yi < nr_groups; ++yi) { //for each group
		
		//we now do summation of interval_l, for cells of size w
		for (int xi = 0; xi < wmax_limit; ++xi) {
			pxi = cell_pointers[xi];
			wmax = wmax_limit - xi;
			
			
			interval_xi_o = ecdf_dictionary[yi        * (wmax_limit+1) + xi]; //the observed number of counts - up to points xi, for the current group

			for (int w = 1; w <= wmax; ++w) { // could be slightly optimized by going over edge intervals in a separate loop from mid intervals
					xj = xi + w;
					pxj = cell_pointers[xj];
					
					
					interval_o = ecdf_dictionary[yi        * (wmax_limit+1) + xj] - interval_xi_o;
					
					interval_e = ((double)uvs_yc[yi]) * ((double)(pxj-pxi)) * nrmlz;
					interval_c = ((interval_o - interval_e) * (interval_o - interval_e)) / interval_e; //note this is without the nrm.count , as shachar used to hold them (see other xdp functions).
					interval_l = ((interval_o > 0) ? (interval_o * log(interval_o / interval_e)) : 0); //same holds for this row.
										
					if (pxi == 0 || pxj == n){
						// edge interval
						if(perform_kahan == 1){ //perform kahan
							temp=interval_c-kahan_c_chi_mk_edge[w];
							kahan_t = sum_chi_w_edge[w] + temp;
							kahan_c_chi_mk_edge[w] = (kahan_t - sum_chi_w_edge[w]) - temp;
							sum_chi_w_edge[w]=kahan_t;
						
							temp=interval_l-kahan_c_like_mk_edge[w];
							kahan_t = sum_like_w_edge[w] + temp;
							kahan_c_like_mk_edge[w] = (kahan_t - sum_like_w_edge[w]) - temp;
							sum_like_w_edge[w]=kahan_t;
						}
						else{
							sum_chi_w_edge[w] += interval_c;
							sum_like_w_edge[w] += interval_l;
						}
					
						
										
					}else {
						//non edge interval
						if(perform_kahan == 1){ //perform kahan
							temp=interval_c-kahan_c_chi_mk[w];
							kahan_t = sum_chi_w[w] + temp;
							kahan_c_chi_mk[w] = (kahan_t - sum_chi_w[w]) - temp;
							sum_chi_w[w]=kahan_t;
						
							temp=interval_l-kahan_c_like_mk[w];
							kahan_t = sum_like_w[w] + temp;
							kahan_c_like_mk[w] = (kahan_t - sum_like_w[w]) - temp;
							sum_like_w[w]=kahan_t;
						}
						else{
							sum_chi_w[w] += interval_c;
							sum_like_w[w] += interval_l;
						}
					}
					
					
					
			}
			
		}
	
	}
	
	for(int i=0;i<wmax_limit+1;i++){
		kahan_c_chi_mk[i] = 0;
		kahan_c_like_mk[i] = 0;
	}
	
	wmax = wmax_limit;//-1;//wmax = n - K - 1;
	for(int k=0;k<K-1;k++){
		for (int w = 1; w <= wmax; ++w){
			if(perform_kahan ==1){
				//sum over chi edge
				temp= sum_chi_w_edge[w] * adp_l_mk[k*wmax+w] -kahan_c_chi_mk[k];
				kahan_t = xdp_sc_mk[k] + temp;
				kahan_c_chi_mk[k] = (kahan_t - xdp_sc_mk[k]) - temp;
				xdp_sc_mk[k]=kahan_t;
				
				//sum over chi no edge
				temp= sum_chi_w[w] * adp_mk[k*wmax+w] -kahan_c_chi_mk[k];
				kahan_t = xdp_sc_mk[k] + temp;
				kahan_c_chi_mk[k] = (kahan_t - xdp_sc_mk[k]) - temp;
				xdp_sc_mk[k]=kahan_t;
			
				//sum over like edge
				temp= sum_like_w_edge[w] * adp_l_mk[k*wmax+w] -kahan_c_like_mk[k];
				kahan_t = xdp_sl_mk[k] + temp;
				kahan_c_like_mk[k] = (kahan_t - xdp_sl_mk[k]) - temp;
				xdp_sl_mk[k]=kahan_t;
				
				//sum over like no edge
				temp= sum_like_w[w] * adp_mk[k*wmax+w] -kahan_c_like_mk[k];
				kahan_t = xdp_sl_mk[k] + temp;
				kahan_c_like_mk[k] = (kahan_t - xdp_sl_mk[k]) - temp;
				xdp_sl_mk[k]=kahan_t;
			}else{
			//code with no kahan summation:
				if(adp_l_mk[k*wmax+w] > 0){
					xdp_sc_mk[k] += sum_chi_w_edge[w] * adp_l_mk[k*wmax+w];
					xdp_sl_mk[k] += sum_like_w_edge[w] * adp_l_mk[k*wmax+w];
				}
				
				if(adp_mk[k*wmax+w] > 0){
					xdp_sc_mk[k] += sum_chi_w[w] * adp_mk[k*wmax+w];
					xdp_sl_mk[k] += sum_like_w[w] * adp_mk[k*wmax+w];
				}
			
			}
			
			
			
		}
	}
	
	
	
	delete[] kahan_c_chi_mk;
	delete[] kahan_c_like_mk;
	delete[] kahan_c_chi_mk_edge;
	delete[] kahan_c_like_mk_edge;
	delete[] sum_chi_w;
	delete[] sum_like_w;
	delete[] sum_chi_w_edge;
	delete[] sum_like_w_edge;
	delete[] cell_pointers;
	delete[] ecdf_dictionary;
	// the 1/n gets us to the MI scale
	for(int k=0;k<K-1;k++){
		xdp_sc_mk[k] = xdp_sc_mk[k]/n;
		xdp_sl_mk[k] = xdp_sl_mk[k]/n;
	}
}

// Univariate independence scores
// ================================================================================================

// NOTE: scores are always computed on (maybe a subset of) uvs_x, uvs_y, uvs_xr, uvs_yr, which are
// of size uvs_n. Other member variables may be used depending on the specific score.
// The following outputs are produced: uvs_sc, uvs_mc, uvs_sl, uvs_ml.

void StatsComputer::uvs_ind_ad(void) {
	// TODO
}

void StatsComputer::uvs_ind_cvm_ks(void) {
	// TODO
}

void StatsComputer::uvs_ind_dcov(void) {
	// TODO
}

// General notes about the univariate distribution-free test (UDF) scores:
// -----------------------------------------------------------------------------
//
// In this test type, the relevant inputs are rx (stored in <dx>) and ry (stored
// in <y> so that it can be permuted the same way as in the 2/k-sample tests).
//
// NOTE: in this C implementation, I expect ranks to be in 0,1,2,...,n-1 (except
// the high-k variants ddp() and adp() which work on 1-based ranks).
//
// The SPR test uses single points in the paired rank-rank sample space as loci
// for examining the local dependence. This is equivalent to Hoeffding's test.
// The PPR test can be seen as direct extension of this idea, which uses pairs
// of points in order to examine local dependence. In this case we can,
// hopefully, model better the distribution of points relative to one another.
// The TPR test goes on to look at triplets of points.
//
// Given the sample size <n>, the inputs <xr> and <yr> (the ranks of some <x>
// and <y> on their original scales) are merely two permutations of the numbers
// 1:n. We thus only need to focus on this <n> over <n> grid. Our statistic
// concerns the number of (x,y) points in the sample that fall in the rectangle
// whose diagonal connects two points from the set and has sides parallel to the
// axes. There are n(n - 1) such rectangles, so merely going over them is O(n^2),
// defining a lower bound for complexity. An algorithm that achieves this bound
// first computes the cumulative sum over the grid, of the indicator of whether
// the point (i,j) is in the set or not. Subsequently we can go over the n(n - 1)
// pairs of points and compute the number of points inside the rectangle they
// define in O(1) using the "double integral" we have computed.
//
// The '33' (3x3) variant is almost the same as the '22' (2,2) variant. The
// difference is that here, for each pair of points from the set, the 3x3
// contingency table for the 3x3 partition of the (x,y) plane is computed,
// whose center cell is defined by the rectangle used in the (2,2) variant.
//
// The "adp" variants are exactly the same as the "ddp" variants, only they
// run over all rank-rank pairs, not just those observed.
//
// -----------------------------------------------------------------------------

void StatsComputer::uvs_ind_ddp2(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all points
	int n = uvs_n;
	int nm1 = n - 1;
	double nm1d = n - 1;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

	int yi, xi;

	for (int i = 0; i < n; ++i) {
		yi = uvs_yr[i];
		xi = uvs_xr[i];

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
		// Can optimize: points with extreme ranks can be filtered beforehand...
		if (xi == 0 || xi == nm1 || yi == 0 || yi == nm1) {
			continue;
		}
#endif

		compute_spr_obs(xi, yi, n, dintegral_pn, nm1, nm1d);
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_ddp3_c(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	int nm1 = n - 1;
#endif
	int nm2 = n - 2;
	double nm2s = nm2 * nm2;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

	int yi, yj, xi, xj;
	int xr_lo, xr_hi, yr_lo, yr_hi;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			yi = uvs_yr[i]; xi = uvs_xr[i];
			yj = uvs_yr[j]; xj = uvs_xr[j];

		    if (xi < xj) {
		      xr_lo = xi;
		      xr_hi = xj;
		    } else {
		      xr_lo = xj;
		      xr_hi = xi;
		    }

		    if (yi < yj) {
		      yr_lo = yi;
		      yr_hi = yj;
		    } else {
		      yr_lo = yj;
		      yr_hi = yi;
		    }

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
			if (xr_lo == 0 || xr_hi == nm1 || yr_lo == 0 || yr_hi == nm1 ||
				xr_hi - xr_lo == 1 || yr_hi - yr_lo == 1)
			{
				continue;
			}
#endif

			compute_ppr_22(xr_lo, xr_hi, yr_lo, yr_hi, dintegral_pn, nm2, nm2s);
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_ddp3(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	int nm1 = n - 1;
#endif
	double nm2 = n - 2;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

	int yi, yj, xi, xj;
	int xr_lo, xr_hi, yr_lo, yr_hi;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			yi = uvs_yr[i]; xi = uvs_xr[i];
			yj = uvs_yr[j]; xj = uvs_xr[j];

		    if (xi < xj) {
		      xr_lo = xi;
		      xr_hi = xj;
		    } else {
		      xr_lo = xj;
		      xr_hi = xi;
		    }

		    if (yi < yj) {
		      yr_lo = yi;
		      yr_hi = yj;
		    } else {
		      yr_lo = yj;
		      yr_hi = yi;
		    }

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
			if (xr_lo == 0 || xr_hi == nm1 || yr_lo == 0 || yr_hi == nm1 ||
				xr_hi - xr_lo <= 1 || yr_hi - yr_lo <= 1)
			{
				continue;
			}
#endif

			compute_ppr_33(xr_lo, xr_hi, yr_lo, yr_hi, n, dintegral_pn, nm2);
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_ddp4(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	int nm1 = n - 1;
#endif
	double nm3 = n - 3;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

	int yi, xi, yj, xj, yk, xk, xl, xm, xh, yl, ym, yh;
	bool fij, fik, fjk;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			for (int k = j + 1; k < n; ++k) {
				yi = uvs_yr[i]; xi = uvs_xr[i];
				yj = uvs_yr[j]; xj = uvs_xr[j];
				yk = uvs_yr[k]; xk = uvs_xr[k];

				fij = (xi < xj); fik = (xi < xk); fjk = (xj < xk);

				if (fij && fjk) {
					xl = xi; xm = xj; xh = xk;
				} else if (fik && !fjk) {
					xl = xi; xm = xk; xh = xj;
				} else if (!fij && fik) {
					xl = xj; xm = xi; xh = xk;
				} else if (fjk && !fik) {
					xl = xj; xm = xk; xh = xi;
				} else if (!fik && fij) {
					xl = xk; xm = xi; xh = xj;
				} else {
					xl = xk; xm = xj; xh = xi;
				}

				fij = (yi < yj); fik = (yi < yk); fjk = (yj < yk);

				if (fij && fjk) {
					yl = yi; ym = yj; yh = yk;
				} else if (fik && !fjk) {
					yl = yi; ym = yk; yh = yj;
				} else if (!fij && fik) {
					yl = yj; ym = yi; yh = yk;
				} else if (fjk && !fik) {
					yl = yj; ym = yk; yh = yi;
				} else if (!fik && fij) {
					yl = yk; ym = yi; yh = yj;
				} else {
					yl = yk; ym = yj; yh = yi;
				}

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
				if (xl == 0 || xh == nm1 || yl == 0 || yh == nm1 ||
					xh - xm <= 1 || xm - xl <= 1 || yh - ym <= 1 || ym - yl <= 1) {
					// partition is degenerate (not sure if I should still allow it to contribute to the overall score
					continue;
				}
#endif

				compute_tpr(xl, xm, xh, yl, ym, yh, n, dintegral_pn, nm3);
			}
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_ddp(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Here we also need to order x by y and vice versa
	compute_ordered_ranks(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;

	uvs_sc = 0;
	uvs_sl = 0;

	// These can't be computed this way
	uvs_mc = 0;
	uvs_ml = 0;

	int rect_o;
	double rect_e, rect_c, rect_l, cnt;
	double edenom = 1.0 / (n - K + 1); // I wonder if the denominator here makes sense for the degenerate partitions
	kahan_c_chi = 0;
	kahan_c_like = 0;
	double	kahan_t;
	double nr_parts = 0;

	double nr_nonempty_cells = 0;

	for (int xl = 1; xl <= n; ++xl) {
		for (int xh = xl; xh <= n; ++xh) {
			for (int yl = 1; yl <= n; ++yl) {
				for (int yh = yl; yh <= n; ++yh) {
					cnt = count_ddp_with_given_cell(xl, xh, yl, yh);
					if (cnt > 0) {
						rect_o = count_sample_points_in_rect(xl, xh, yl, yh);
						rect_e = (xh - xl + 1) * (yh - yl + 1) * edenom;
						rect_c = ((rect_o - rect_e) * (rect_o - rect_e)) / rect_e * cnt - kahan_c_chi;
						rect_l = ((rect_o > 0) ? (rect_o * log(rect_o / rect_e)) : 0) * cnt - kahan_c_like;

						kahan_t = uvs_sc + rect_c;
						kahan_c_chi = (kahan_t - uvs_sc) - rect_c;
						uvs_sc = kahan_t;

						kahan_t = uvs_sl + rect_l;
						kahan_c_like = (kahan_t - uvs_sl) - rect_l;
						uvs_sl = kahan_t;

						nr_parts += cnt;

						if (rect_o > 0) {
							nr_nonempty_cells += cnt;
						}
					}
				}
			}
		}
	}

#ifdef XDP_NORMALIZE
  // NOTE: since there are degenerate partitions, nr_parts is not a multiple of
  // K^2... this is irrelevant for testing but may be relevant for estimating
  // mutual information (which DDP is not really intended for)

	nr_parts /= (K * K);

	if (correct_mi_bias) {
		double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
		uvs_sc += mm_bias;
		uvs_sl += mm_bias;
	}

	double normalizer = nr_parts * n; // the *extra* n factor get us to the MI scale.

	uvs_sc /= normalizer;
	uvs_sl /= normalizer;
#endif
}

void StatsComputer::uvs_ind_adp2(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all points
	int n = uvs_n;
	double nd = n;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	for (int xi = 1; xi < n; ++xi) {
		for (int yi = 1; yi < n; ++yi) {
#else
	int np1 = n + 1;
	for (int xi = 0; xi < np1; ++xi) {
		for (int yi = 0; yi < np1; ++yi) {
#endif
			compute_spr_all(xi, yi, n, dintegral_pn, nd);
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_adp3_c(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;
	int nm1 = n - 1;
	int nm2 = n - 2;
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	int nm3 = n - 3;
#endif
	double nm2s = nm2 * nm2;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	for (int xr_lo = 1; xr_lo < nm3; ++xr_lo) {
		for (int xr_hi = xr_lo + 2; xr_hi < nm1; ++xr_hi) {
			for (int yr_lo = 1; yr_lo < nm3; ++yr_lo) {
				for (int yr_hi = yr_lo + 2; yr_hi < nm1; ++yr_hi) {
#else
	for (int xr_lo = 0; xr_lo < nm1; ++xr_lo) {
		for (int xr_hi = xr_lo + 1; xr_hi < n; ++xr_hi) {
			for (int yr_lo = 0; yr_lo < nm1; ++yr_lo) {
				for (int yr_hi = yr_lo + 1; yr_hi < n; ++yr_hi) {
#endif
					// actually... this must be replaced with a special version that works on half ranks, see spr_all
					compute_ppr_22(xr_lo, xr_hi, yr_lo, yr_hi, dintegral_pn, nm2, nm2s);
				}
			}
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_adp3(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;
	int nm1 = n - 1;
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	int nm3 = n - 3;
#endif
	double nm2 = n - 2;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	for (int xr_lo = 1; xr_lo < nm3; ++xr_lo) {
		for (int xr_hi = xr_lo + 2; xr_hi < nm1; ++xr_hi) {
			for (int yr_lo = 1; yr_lo < nm3; ++yr_lo) {
				for (int yr_hi = yr_lo + 2; yr_hi < nm1; ++yr_hi) {
#else
	for (int xr_lo = 0; xr_lo < nm1; ++xr_lo) {
		for (int xr_hi = xr_lo + 1; xr_hi < n; ++xr_hi) {
			for (int yr_lo = 0; yr_lo < nm1; ++yr_lo) {
				for (int yr_hi = yr_lo + 1; yr_hi < n; ++yr_hi) {
#endif
					// actually... this must be replaced with a special version that works on half ranks, see spr_all
					compute_ppr_33(xr_lo, xr_hi, yr_lo, yr_hi, n, dintegral_pn, nm2);
				}
			}
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_adp4(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;
	int nm1 = n - 1;
	double nm3d = n - 3;

	uvs_sc = 0;
	uvs_mc = 0;
	uvs_sl = 0;
	uvs_ml = 0;

	ng_chi  = 0;
	ng_like = 0;

	// Can collect individual tables if wanted

#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
	int nm3 = n - 3;
	int nm5 = n - 5;

	for (int xl = 1; xl < nm5; ++xl) {
		for (int xm = xl + 2; xm < nm3; ++xm) {
			for (int xh = xm + 2; xh < nm1; ++xh) {
				for (int yl = 1; yl < nm5; ++yl) {
					for (int ym = yl + 2; ym < nm3; ++ym) {
						for (int yh = ym + 2; yh < nm1; ++yh) {
#else
	int nm2 = n - 2;

	for (int xl = 0; xl < nm2; ++xl) {
		for (int xm = xl + 1; xm < nm1; ++xm) {
			for (int xh = xm + 1; xh < n; ++xh) {
				for (int yl = 0; yl < nm2; ++yl) {
					for (int ym = yl + 1; ym < nm1; ++ym) {
						for (int yh = ym + 1; yh < n; ++yh) {
#endif
							// actually... this must be replaced with a special version that works on half ranks, see spr_all
							compute_tpr(xl, xm, xh, yl, ym, yh, n, dintegral_pn, nm3d);
						}
					}
				}
			}
		}
	}

#ifdef XDP_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	uvs_sc /= ng_chi;
	uvs_sl /= ng_like;
#endif
}

void StatsComputer::uvs_ind_adp(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(uvs_n, uvs_xr, uvs_yr);

	// Now compute the score using all rectangles
	int n = uvs_n;

	uvs_sc  = 0;
	uvs_sl = 0;

	// These can't be computed this way
	uvs_mc = 0;
	uvs_ml = 0;

	int rect_o, xh, yh;
	double cnt, rect_e, rect_c, rect_l;
	//double edenom = 1.0 / (n - K + 1); // I wonder if the denominator here makes sense for the degenerate partitions
	double edenom = 1.0 / n; // This makes much more sense (and is for partitions between ranks rather than at ranks)

	double nr_nonempty_cells = 0;

	// It is more stable numerically to accumulate the statistics in the order
	// of rectangle area, this is why I use the (xl, yl, w, h) enumeration.
	// While this helps stability, it may adversely affect memory locality...
	// Maybe with the Kahan summation, this is no longer necessary?

	kahan_c_chi = 0;
	kahan_c_like = 0;
	double kahan_t;

	for (int w = 1; w <= n; ++w) {
		for (int h = 1; h <= n; ++h) {
			for (int xl = 1; xl <= n - w + 1; ++xl) {
				for (int yl = 1; yl <= n - h + 1; ++yl) {
					xh = xl + w - 1;
					yh = yl + h - 1;
					cnt = count_adp_with_given_cell(xl, xh, yl, yh);
					if (cnt > 0) {
						// (in theory I could compute this for multiple K's at
						// once with relatively small overhead)
						rect_o = count_sample_points_in_rect(xl, xh, yl, yh);
						rect_e = w * h * edenom;
						rect_c = ((rect_o - rect_e) * (rect_o - rect_e)) / rect_e * cnt - kahan_c_chi;
						rect_l = ((rect_o > 0) ? (rect_o * log(rect_o / rect_e)) : 0) * cnt - kahan_c_like;

						kahan_t = uvs_sc + rect_c;
						kahan_c_chi = (kahan_t - uvs_sc) - rect_c;
						uvs_sc = kahan_t;

						kahan_t = uvs_sl + rect_l;
						kahan_c_like = (kahan_t - uvs_sl) - rect_l;
						uvs_sl = kahan_t;

						if (rect_o > 0) {
							nr_nonempty_cells += cnt;
						}
					}
				}
			}
		}
	}

#ifdef XDP_NORMALIZE
	double nr_parts = choose(n - 1, K - 1); // this is actually only the sqrt of the nr parts
	nr_parts *= nr_parts;

	if (correct_mi_bias) {
		double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
		uvs_sc += mm_bias;
		uvs_sl += mm_bias;
	}

	double normalizer = nr_parts * n; // gets us to the MI scale
	uvs_sc /= normalizer;
	uvs_sl /= normalizer;
#endif
}

void StatsComputer::uvs_ind_adp_mk(void) {
///*
	int perform_kahan = 0;
	long n = (long)uvs_n;
	
	
	uvs_sc  = 0;
	uvs_sl = 0;

	// These can't be computed this way
	uvs_mc = 0;
	uvs_ml = 0;

	long rect_o, xh, yh;
	double rect_e, rect_c, rect_l,rect_o_double;
	
	double edenom = 1.0 / n; // This makes much more sense (and is for partitions between ranks rather than at ranks)

	// First, compute the double integral (padded with a leading row and column of zeros)

	long wmax_limit = n;
	if(equipartition_type == 0){
		wmax_limit= n; //wmax can't reach more than this number
		compute_double_integral(uvs_n, uvs_xr, uvs_yr);
	}else if(equipartition_type ==1){
		wmax_limit= equipartition_nr_cells_m;//(int)(sqrt((double)(n))); //wmax can't reach more than this number
		//compute_double_integral(uvs_n, uvs_xr, uvs_yr);
		compute_double_integral_eqp(uvs_n, uvs_xr, uvs_yr,equipartition_nr_cells_m);
	}
	
	// Now compute the score using all rectangles

	
	//double nr_nonempty_cells = 0;
	long dim_size=wmax_limit+1; 
	double* cells_likelihood_array = new double[dim_size*dim_size*9];// for holding the counts by cell type and window sizes;
	double* cells_likelihood_kahan = new double[dim_size*dim_size*9];
	double* cells_chisq_array = new double[dim_size*dim_size*9];
	double* cells_chisq_kahan = new double[dim_size*dim_size*9];
	for(long i=0;i<dim_size*dim_size*9;i++){
		cells_likelihood_array[i]=0;
		cells_likelihood_kahan[i]=0;
		cells_chisq_array[i]=0;
		cells_chisq_kahan[i]=0;
	}
	
	
	long* cell_pointers = new long[wmax_limit+1];
	
	if(equipartition_type == 0){
		for(long i=0;i<wmax_limit+1;i++){
			cell_pointers[i]=i+1;
		}
		
	}else if(equipartition_type ==1){
	
		for(long i=0;i<wmax_limit+1;i++){
			if(i != wmax_limit){
				cell_pointers[i]=(long)(((double) (i*n))/((double)wmax_limit))+1; 
			}else{
				cell_pointers[i]=n+1;
			}
		}
	}
	
	

	double kahan_c_chi = 0, kahan_c_like = 0, kahan_t;
	long current_index=0;
	int current_cell_type=0;
	long pxl,pxh,pyl,pyh;
	long w_actual,h_actual;

	for (long w = 1; w <= wmax_limit; ++w){
		for (long h = 1; h <= wmax_limit; ++h){
			for (long xl = 1; xl <= wmax_limit - w + 1 ; ++xl){ 
				pxl = cell_pointers[xl - 1];
				for (long yl = 1; yl <= wmax_limit - h + 1 ; ++yl){ 
					pyl = cell_pointers[yl - 1];
										
					xh = xl + w; 
										
					yh = yl + h; 
										
					pxh = cell_pointers[xh - 1];
					pyh = cell_pointers[yh - 1];
					if (true) { //was if cnt>0
						current_cell_type = compute_adp_mk_cell_type(pxl, pxh -1 , pyl, pyh -1, uvs_n);
						current_index = current_cell_type *dim_size*dim_size + (w-1)*dim_size +h-1;
						if(equipartition_type == 0){
						rect_o = count_sample_points_in_rect(pxl , pxh -1, pyl , pyh - 1); 
						}else{ //(equipartition_type == 1) in this case
						rect_o =  count_sample_points_in_rect_eqp(xl , xh-1 , yl , yh -1 );
						}
						rect_o_double = (double)rect_o;
						w_actual = (pxh - pxl);
						h_actual = (pyh - pyl);
						rect_e = ((double)(w_actual)) * ((double)h_actual) * edenom;
						rect_c = ((rect_o_double - rect_e) * (rect_o_double - rect_e)) / rect_e - cells_chisq_kahan[current_index];
						rect_l = ((rect_o_double > 0) ? (rect_o_double * log(rect_o_double / rect_e)) : 0) - cells_likelihood_kahan[current_index];
						if(perform_kahan == 1){ //kahan summation
							kahan_t = cells_chisq_array[current_index] + rect_c;
							cells_chisq_kahan[current_index] = (kahan_t - cells_chisq_array[current_index]) - rect_c;
							cells_chisq_array[current_index] = kahan_t;

							kahan_t = cells_likelihood_array[current_index] + rect_l;
							cells_likelihood_kahan[current_index] = (kahan_t - cells_likelihood_array[current_index]) - rect_l;
							cells_likelihood_array[current_index] = kahan_t;
						}else{
							cells_chisq_array[current_index] +=  ((rect_o_double - rect_e) * (rect_o_double - rect_e)) / rect_e;
							cells_likelihood_array[current_index] += ((rect_o_double > 0) ? (rect_o_double * log(rect_o_double / rect_e)) : 0);
						}
						
					}
				}
			}
		}
	}
	
	double current_lrt;
	double current_chisq;
	double current_count;
	double current_add;
	double normalizer;
	normalizer = n; //nr_parts * n; // gets us to the MI scale

	for(int k=0;k<adp_mk_tables_nr;k++){
		current_lrt = 0;
		current_chisq = 0;
		kahan_c_chi = 0;
		kahan_c_like = 0;
		for(int cell_type=0;cell_type<9;cell_type++){
			for (int w = 1; w <= wmax_limit; w++) {
				for (int h = 1; h <= wmax_limit; h++) {
						current_index = cell_type *dim_size*dim_size + (w-1)*dim_size +h-1;
						current_count = count_adp_mk_cell_type(adp_mk_tables_m[k], adp_mk_tables_l[k], cell_type ,w ,h,wmax_limit);
						
						if((perform_kahan ==1) && current_count>0){ //kahan summation
						current_add = cells_chisq_array[current_index] * current_count;
						
						kahan_t = current_chisq + current_add - kahan_c_chi;
						kahan_c_chi = (kahan_t - current_add) - current_chisq;
						current_chisq = kahan_t;
											
						current_add = cells_likelihood_array[current_index] * current_count;
						
						kahan_t = current_lrt + current_add - kahan_c_like;
						kahan_c_like = (kahan_t - current_add) - current_lrt;
						current_lrt = kahan_t;
						
						}else if (current_count>0){
						current_add = cells_chisq_array[current_index] * current_count;
						current_chisq += current_add;
						
						current_add = cells_likelihood_array[current_index] * current_count;
						current_lrt += current_add;
						}
						
			}
			}
		}
		
		adp_ind_sc_mk[k] = current_chisq;
		adp_ind_sl_mk[k] = current_lrt;
		
		adp_ind_sc_mk[k] /= normalizer;
		adp_ind_sl_mk[k] /= normalizer;
		
	}
				
#ifdef XDP_NORMALIZE
	// MI bias correction with MM method is currently only available in single partition size function
	/* double nr_parts = choose(n - 1, K - 1); // this is actually only the sqrt of the nr parts
	nr_parts *= nr_parts;

	if (correct_mi_bias) {
		double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
		uvs_sc += mm_bias;
		uvs_sl += mm_bias;
	}*/

	//double normalizer = nr_parts * n; // gets us to the MI scale
	//uvs_sc /= normalizer;
	//uvs_sl /= normalizer;
#endif
//*/

uvs_sc = -1;
uvs_sl = -1;
uvs_mc = -1;
uvs_ml = -1;


delete[] cells_likelihood_array;
delete[] cells_likelihood_kahan;
delete[] cells_chisq_array;
delete[] cells_chisq_kahan;
delete[] cell_pointers;
}

// Misc. helper functions
// ================================================================================================

void StatsComputer::sort_xy_distances_per_row(void) {
	// Could be parallelized.

	for (int k = 0; k < xy_nrow; ++k) {
		for (int l = 0; l < xy_nrow; ++l) {
			sorted_dx_gen[k][l].x = dx[l * xy_nrow + k];
			sorted_dx_gen[k][l].y = dy[idx_perm[l] * xy_nrow + idx_perm[k]];
			sorted_dx_gen[k][l].i = l;
		}

		sort(sorted_dx_gen[k].begin(), sorted_dx_gen[k].end(), dbl_dbl_int_pair_comparator_xy);
	}
}



void StatsComputer::accumulate_2x2_contingency_table(double a00, double a01, double a10, double a11, double nrmlz, double reps) {
	double e00, e01, e10, e11, emin00_01, emin10_11, emin, current_chi, current_like;

	e00 = ((a00 + a01) * (a00 + a10)) * nrmlz;
	e01 = ((a00 + a01) * (a01 + a11)) * nrmlz;
	e10 = ((a10 + a11) * (a00 + a10)) * nrmlz;
	e11 = ((a10 + a11) * (a01 + a11)) * nrmlz;

	emin00_01 = min(e00, e01);
	emin10_11 = min(e10, e11);
	emin = min(emin00_01, emin10_11);

	if (emin > min_w) {
		current_chi = (a00 - e00) * (a00 - e00) / e00
					+ (a01 - e01) * (a01 - e01) / e01
					+ (a10 - e10) * (a10 - e10) / e10
					+ (a11 - e11) * (a11 - e11) / e11;
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi * reps;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((a00 > 0) ? (a00 * log(a00 / e00)) : 0)
				 + ((a01 > 0) ? (a01 * log(a01 / e01)) : 0)
				 + ((a10 > 0) ? (a10 * log(a10 / e10)) : 0)
				 + ((a11 > 0) ? (a11 * log(a11 / e11)) : 0);

	sum_like += current_like * reps;

	if (current_like > max_like) {
		max_like = current_like;
	}
}

void StatsComputer::accumulate_local_stats(double chi, double like, double emin) {
	double kahan_t, kahan_chi, kahan_like;

	if (emin > w_sum) {
		kahan_chi = chi - kahan_c_chi;
		kahan_t = uvs_sc + kahan_chi;
		kahan_c_chi = (kahan_t - uvs_sc) - kahan_chi;
		uvs_sc = kahan_t;
		++ng_chi;
	}

	if ((emin > w_max) && (chi > uvs_mc)) {
		uvs_mc = chi;
	}

	kahan_like = like - kahan_c_like;
	kahan_t = uvs_sl + kahan_like;
	kahan_c_like = (kahan_t - uvs_sl) - kahan_like;
	uvs_sl = kahan_t;
	++ng_like;

	if (like > uvs_ml) {
		uvs_ml = like;
	}
}

void StatsComputer::hhg_gen_inversions(int *permutation, int *source, int *inversion_count, int dim) {
    if (dim <= 1) {
        return;
    } else {
    	hhg_gen_inversions(permutation, source, inversion_count, dim / 2);
    	hhg_gen_inversions(&permutation[dim / 2], &source[dim / 2], inversion_count, dim - dim / 2);
    	hhg_gen_merge(permutation, source, inversion_count, dim);
    }
}

void StatsComputer::hhg_gen_merge(int *permutation, int *source, int *inversion_count, int dim) {
    int left_index = 0, right_index = 0;
    int i, half_dim = dim / 2, nleft = half_dim, nright = dim - half_dim;

    int* left = hhg_gen_left_buffer;
    int* right = hhg_gen_right_buffer;
    int* left_source = hhg_gen_left_source_buffer;
    int* right_source = hhg_gen_right_source_buffer;

    for (i = 0; i < half_dim; i++) {
        left[i] = permutation[i];
        left_source[i] = source[i];
        right[i] = permutation[i + half_dim];
        right_source[i] = source[i + half_dim];
    }

    if (nleft < nright) {
        right[i] = permutation[i + half_dim];
        right_source[i] = source[i + half_dim];
    }

    for (i = 0; i < dim; i++) {
        if ((left_index < half_dim) && (right_index < dim - half_dim)) {
             if (left[left_index] <= right[right_index]) { // I added "=" in order to support ties
                permutation[i] = left[left_index];
                source[i] = left_source[left_index];
                left_index++;
            } else {
                permutation[i] = right[right_index];
                source[i] = right_source[right_index];
                inversion_count[source[i]] += (half_dim - left_index);
                right_index++;
            }
        } else {
            if (left_index < half_dim) {
                permutation[i] = left[left_index];
                source[i] = left_source[left_index];
                left_index++;
            }

            if (right_index < dim - half_dim) {
                permutation[i] = right[right_index];
                source[i] = right_source[right_index];
                right_index++;
            }
        }
    }
}

void StatsComputer::compute_ordered_ranks(int n, double* xx, int* yy) {
	for (int i = 0; i < n; ++i) {
		y_ordered_by_x[(int)(xx[i]) - 1] = yy[i];
		x_ordered_by_y[      yy[i]  - 1] = xx[i];
	}
}

void StatsComputer::compute_single_integral(int n, double* xx, int* yy) {
	memset(double_integral, 0, sizeof(int) * (nr_groups + 1) * dintegral_pn);

	// Populate the padded matrix with the indicator variables of whether there
	// is a point in the set {(x_k, y_k)}_k=1^n in the grid slot (i, j) for i in 1,2,...,n and j in 1,2,...,nr_groups
	// NOTE: the last row holds the integral over all x, ignoring y
	for (int i = 0; i < n; ++i) {
		int yi = yy[i]; // assumed in 0..nr_groups-1, many ties, and won't be padded
		int xi = xx[i]; // assumed in 1..n, no ties, and will be padded
		double_integral[yi        * dintegral_pn + xi] = 1;
		double_integral[nr_groups * dintegral_pn + xi] = 1;
	}

	// Then run linearly and compute the integral row by row
	for (int k = 0; k <= nr_groups; ++k) {
		int row_running_sum = 0;
		for (int i = 1; i < dintegral_pn; ++i) {
			row_running_sum += double_integral[k * dintegral_pn + i];
			double_integral[k * dintegral_pn + i] = row_running_sum;
		}
	}
}

void StatsComputer::compute_double_integral(int n, double* xx, int* yy) {
	memset(double_integral, 0, sizeof(int) * dintegral_pn * dintegral_pn);

	// Populate the padded matrix with the indicator variables of whether there
	// is a point in the set {(x_k, y_k)}_k=1^n in the grid slot (i, j) for i,j in 1,2,...,n
	for (int i = 0; i < n; ++i) {
		int yi = yy[i] + dintegral_zero_based_idxs;
		int xi = xx[i] + dintegral_zero_based_idxs;
		double_integral[yi * dintegral_pn + xi] = 1;
	}

	// Then run linearly and compute the integral in one row-major pass
	int la = dintegral_pn;
	for (int i = 1; i < dintegral_pn; ++i) {
		int row_running_sum = 0;
		++la;
		for (int j = 1; j < dintegral_pn; ++j) {
			row_running_sum += double_integral[la];
			double_integral[la] = row_running_sum + double_integral[la - dintegral_pn];
			++la;
		}
	}
}

void StatsComputer::compute_double_integral_eqp(long n, double* xx, int* yy,long nr_atoms){
	dintegral_pn_eqp = nr_atoms+1;
	memset(double_integral_eqp, 0, sizeof(int) * dintegral_pn_eqp * dintegral_pn_eqp);

	// Populate the padded matrix with the indicator variables of whether there
	// is a point in the set {(x_k, y_k)}_k=1^n in the grid slot (i, j) for i,j in 1,2,...,n
	long yi;
	long xi;
	for (long i = 0; i < n; ++i) {
		yi = (long)ceil(((double)(yy[i]*nr_atoms))/((double)(n))) + dintegral_zero_based_idxs;
		xi = (long)ceil(((double)(xx[i]*nr_atoms))/((double)(n))) + dintegral_zero_based_idxs;
		//REprintf("i: %d xi: %d yi: %d \n",i,xi,yi); //HHHHH check that inverse of atoms find is good
		double_integral_eqp[yi * dintegral_pn_eqp + xi] += 1;
	}

	// Then run linearly and compute the integral in one row-major pass
	long la = dintegral_pn_eqp;
	for (long i = 1; i < dintegral_pn_eqp; ++i) {
		long row_running_sum = 0;
		++la;
		for (long j = 1; j < dintegral_pn_eqp; ++j) {
			row_running_sum += double_integral_eqp[la];
			double_integral_eqp[la] = row_running_sum + double_integral_eqp[la - dintegral_pn_eqp];
			//REprintf("i: %d j: %d tot: %d \n",i,j,double_integral_eqp[la]);
			++la;
		}
	}
}

int StatsComputer::count_sample_points_in_rect(int xl, int xh, int yl, int yh) {
	return (  double_integral[ yh      * dintegral_pn + xh    ]
	        - double_integral[ yh      * dintegral_pn + xl - 1]
	        - double_integral[(yl - 1) * dintegral_pn + xh    ]
	        + double_integral[(yl - 1) * dintegral_pn + xl - 1]);
}


long StatsComputer::count_sample_points_in_rect_eqp(long xl, long xh, long yl, long yh) {
	return (  double_integral_eqp[ yh      * dintegral_pn_eqp + xh    ]
	        - double_integral_eqp[ yh      * dintegral_pn_eqp + xl - 1]
	        - double_integral_eqp[(yl - 1) * dintegral_pn_eqp + xh    ]
	        + double_integral_eqp[(yl - 1) * dintegral_pn_eqp + xl - 1]);
}



// This counts the number of data driven partitions of the data, that contain the given rectangle as a cell.
double StatsComputer::count_ddp_with_given_cell(int xl, int xh, int yl, int yh) {
	int n = uvs_n, m;

	// This may not be necessary in practical usage scenarios
	if ((xl == 1 && xh == n) || (yl == 1 && yh == n)) {
		return (0);
	}

	if (xl == 1 && yl == 1) {
		// ll corner cell
		if ((y_ordered_by_x[xh] > yh) && (x_ordered_by_y[yh] > xh)) {
			m = count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if (x_ordered_by_y[yh] == xh + 1) {
				// The hh point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xh + 1, and the point
				// with y = yh + 1, which do not split the rectangle, and the remaining
				// K-3 points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (xl == 1 && yh == n) {
		// lh corner cell
		if ((y_ordered_by_x[xh] < yl) && (x_ordered_by_y[yl - 2] > xh)) {
			m = count_sample_points_in_rect(xh + 2, n, 1, yl - 2);

			if (x_ordered_by_y[yl - 2] == xh + 1) {
				// The hl point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xh + 1, and the point
				// with y = yl - 1, which do not split the rectangle, and the remaining K-3
				// points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (xh == n && yl == 1) {
		// hl corner cell
		if ((y_ordered_by_x[xl - 2] > yh) && (x_ordered_by_y[yh] < xl)) {
			m = count_sample_points_in_rect(1, xl - 2, yh + 2, n);

			if (x_ordered_by_y[yh] == xl - 1) {
				// The lh point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xl - 1, and the point
				// with y = yh + 1, which do not split the rectangle, and the remaining K-3
				// points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (xh == n && yh == n) {
		// hh corner cell
		if ((y_ordered_by_x[xl - 2] < yl) && (x_ordered_by_y[yl - 2] < xl)) {
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2);

			if (x_ordered_by_y[yl - 2] == xl - 1) {
				// The ll point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xl - 1, and the point
				// with y = yl - 1, which do not split the rectangle, and the remaining K-3
				// points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (yl == 1) {
		// bottom edge cell
		if (  (x_ordered_by_y[yh    ] < xl || x_ordered_by_y[yh] > xh)
			&& y_ordered_by_x[xl - 2] > yh && y_ordered_by_x[xh] > yh)
		{
			m = count_sample_points_in_rect(1, xl - 2, yh + 2, n)
			  + count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if ((x_ordered_by_y[yh] == xl - 1) || (x_ordered_by_y[yh] == xh + 1)) {
				// One of the top corners is in the sample, so this point has to be
				// chosen, so does the point with the remaining x boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with y = yh + 1, the point with x = xl - 1, and
				// the point with x = xh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else if (yh == n) {
		// top edge cell
		if (  (x_ordered_by_y[yl - 2] < xl || x_ordered_by_y[yl - 2] > xh)
			&& y_ordered_by_x[xl - 2] < yl && y_ordered_by_x[xh    ] < yl) {
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2)
			  + count_sample_points_in_rect(xh + 2, n, 1, yl - 2);

			if ((x_ordered_by_y[yl - 2] == xl - 1) || (x_ordered_by_y[yl - 2] == xh + 1)) {
				// One of the bottom corners are in the sample, so this point has to be
				// chosen, so does the point with the remaining x boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with y = yh + 1, the point with x = xl - 1, and
				// the point with x = xh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else if (xl == 1) {
		// left edge cell
		if (  (y_ordered_by_x[xh    ] < yl || y_ordered_by_x[xh] > yh)
			&& x_ordered_by_y[yl - 2] > xh && x_ordered_by_y[yh] > xh) {
			m = count_sample_points_in_rect(xh + 2, n, 1, yl - 2)
			  + count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if ((y_ordered_by_x[xh] == yl - 1) || (y_ordered_by_x[xh] == yh + 1)) {
				// One of the right corners are in the sample, so this point has to be
				// chosen, so does the point with the remaining y boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with x = xh + 1, the point with y = yl - 1, and
				// the point with y = yh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else if (xh == n) {
		// right edge cell
		if (  (y_ordered_by_x[xl - 2] < yl || y_ordered_by_x[xl - 2] > yh)
			&& x_ordered_by_y[yl - 2] < xl && x_ordered_by_y[yh    ] < xl)
		{
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2)
			  + count_sample_points_in_rect(1, xl - 2, yh + 2, n);

			if ((y_ordered_by_x[xl - 2] == yl - 1) || (y_ordered_by_x[xl - 2] == yh + 1)) {
				// One of the left corners are in the sample, so this point has to be
				// chosen, so does the point with the remaining y boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with x = xl - 1, the point with y = yl - 1, and
				// the point with y = yh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else {
		// inner cell
		if (   (y_ordered_by_x[xl - 2] < yl || y_ordered_by_x[xl - 2] > yh)
			&& (y_ordered_by_x[xh    ] < yl || y_ordered_by_x[xh    ] > yh)
			&& (x_ordered_by_y[yl - 2] < xl || x_ordered_by_y[yl - 2] > xh)
			&& (x_ordered_by_y[yh    ] < xl || x_ordered_by_y[yh    ] > xh))
		{
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2)
			  + count_sample_points_in_rect(1, xl - 2, yh + 2, n)
			  + count_sample_points_in_rect(xh + 2, n, 1, yl - 2)
			  + count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if (   ((x_ordered_by_y[yl - 2] == xl - 1) && (x_ordered_by_y[yh] == xh + 1))
				|| ((x_ordered_by_y[yl - 2] == xh + 1) && (x_ordered_by_y[yh] == xl - 1))) {
				// Two opposing corners are sample points. These must be chosen, and the
				// rest K-3 points are chosen from m.
				return (choose(m, K - 3));
			} else if (   (x_ordered_by_y[yl - 2] == xl - 1) || (x_ordered_by_y[yl - 2] == xh + 1)
					   || (x_ordered_by_y[yh    ] == xl - 1) || (x_ordered_by_y[yh    ] == xh + 1)) {
				// Exactly one corner is a sample point and has to be chosen. We also have
				// to choose the two points with x and y coordinates of the opposing
				// corner. The remaining K-4 points are chosen from m.
				return (choose(m, K - 4));
			} else {
				// No corner is a sample point, so we must choose the four points with
				// the bounding x's and y's (which we can). The remaining K-5 points are
				// chosen from m.
				return (choose(m, K - 5));
			}
		}
	}

	// If we reached here, then the specified rectangle can never be a cell in any
	// data driven partition of the given dataset
	return (0);
}

// This counts the number of unrestricted partitions of the data.
double StatsComputer::count_adp_with_given_cell(int xl, int xh, int yl, int yh) {
#if 0
	// This is quite silly, and I can probably speed this up considerably
	// by blocking in a way that takes cache hierarchy into account.
	return (adp[(xl - 1) * uvs_n + xh - 1] * adp[(yl - 1) * uvs_n + yh - 1]);
#else
	double cx, cy;

	if (xl == 1) {
		// left anchored interval
		cx = adp_l[xh - 1];
	} else if (xh == uvs_n) {
		// right anchored interval
		cx = adp_r[xl-1];
	} else {
		cx = adp[xh - xl];
	}

	if (yl == 1) {
		// left anchored interval
		cy = adp_l[yh - 1];
	} else if (yh == uvs_n) {
		// right anchored interval
		cy = adp_r[yl-1];
	} else {
		cy = adp[yh - yl];
	}

	return (cx * cy);
#endif
}


// This counts the number of unrestricted partitions of the data.
int StatsComputer::compute_adp_mk_cell_type(long xl, long xh, long yl, long yh,long nr_atoms) {
	int cx, cy;
	int rc=-1;
	if (xl == 1) {
		// left anchored interval
		cx =1;
	} else if (xh == nr_atoms) {
		// right anchored interval
		cx = 3;
	} else {
		cx = 2;
	}

	if (yl == 1) {
		// left anchored interval
		cy = 1;
	} else if (yh == nr_atoms) {
		// right anchored interval
		cy = 3;
	} else {
		cy = 2;
	}
	if(cy == 1 && cx == 1){
		rc=0;
	}else if(cy == 1 && cx == 2){
		rc=1;
	}else if(cy == 1 && cx == 3){
		rc=2;
	}else if(cy == 2 && cx == 1){
		rc=3;
	}else if(cy == 2 && cx == 2){
		rc=4;
	}else if(cy == 2 && cx == 3){
		rc=5;
	}else if(cy == 3 && cx == 1){
		rc=6;
	}else if(cy == 3 && cx == 2){
		rc=7;
	}else if(cy == 3 && cx == 3){
		rc=8;	
	}
	
	return (rc);

}

double StatsComputer::count_adp_mk_cell_type(int M,int L,int type, long w, long h,long nr_atoms) {
	double cx,cy;
	
	if (type == 0 || type == 3 || type ==6) {
		// left anchored interval
		cx = adp_l_mk[(M-2)*nr_atoms + w-1];//[xh - 1];
	} else if (type == 2 || type == 5 || type ==8) {
		// right anchored interval
		cx = adp_r_mk[(M-2)*nr_atoms + w-1];
	} else {
		cx = adp_mk[(M-2)*nr_atoms + w-1];
	}

	if (type == 0 || type ==1 || type ==2) {
		// left anchored interval
		cy = adp_l_mk[(L-2)*nr_atoms + h-1];
	} else if (type ==6 || type ==7 || type ==8) {
		// right anchored interval
		cy = adp_r_mk[(L-2)*nr_atoms + h-1];
	} else {
		cy = adp_mk[(L-2)*nr_atoms + h-1];
	}
	
	return((double)cx*cy);
}

double StatsComputer::compute_ht(void) {
	int n = xy_nrow;
    int x_c_i, total_other_y_count_so_far;
	double expected_total_other_y_count_so_far, nm1i = 1.0 / (n - 1);
	double ht0 = 0, ht1 = 0, tmp;
    int i, j;
	int n0 = y_counts[0];
	int n1 = y_counts[1];

	// NOTE: I'll compute the HT test in the version they refer to as "T", with
	// "gamma" being 2, and with "w" being all ones. They also have an "S"
	// version which is comparable maybe with the HHG "max" variants.
	// Also, for the moment, ties are not handled.

	// We already have the y1_idx and y0_idx from above, which define the
	// permuted partitioning of "Z" to "X" and "Y". We also already have sorted
	// dx's (these are the distances between the pooled sample "Z").

    for (i = 0; i < n0; i++) {
        x_c_i = y0_idx[i];
        total_other_y_count_so_far = 0; // this is "Mij" in the HT paper

    	for (j = 0; j < n1; j++) {
    		total_other_y_count_so_far += (y_perm[(*sorted_dx)[x_c_i][j].second] == 1);
    		expected_total_other_y_count_so_far = n1 * j * nm1i;
    		tmp = (total_other_y_count_so_far - expected_total_other_y_count_so_far);
    		ht0 += tmp * tmp;
		}
    }

    for (i = 0; i < n1; i++) {
        x_c_i = y1_idx[i];
        total_other_y_count_so_far = 0; // this is "Nij" in the HT paper

    	for (j = 0; j < n0; j++) {
    		total_other_y_count_so_far += (y_perm[(*sorted_dx)[x_c_i][j].second] == 0);
    		expected_total_other_y_count_so_far = n0 * j * nm1i;
    		tmp = (total_other_y_count_so_far - expected_total_other_y_count_so_far);
    		ht1 += tmp * tmp;
		}
    }

    return (ht0 / n0 + ht1 / n1);
}

double StatsComputer::compute_edist(void) {
	int n0 = y_counts[0];
	int n1 = y_counts[1];

	double sum01 = 0;
	double sum00 = 0;
	double sum11 = 0;

	for (int i = 0; i < n0; ++i) {
		for (int j = 0; j < n1; ++j) {
			sum01 += dx[y1_idx[j] * xy_nrow + y0_idx[i]];
		}
	}

	for (int i = 0; i < n0; ++i) {
		for (int j = 0; j < n0; ++j) {
			sum00 += dx[y0_idx[j] * xy_nrow + y0_idx[i]];
		}
	}

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n1; ++j) {
			sum11 += dx[y1_idx[j] * xy_nrow + y1_idx[i]];
		}
	}

	double ret = (2.0 / (n0 * n1)) * sum01 - (1.0 / (n0 * n0)) * sum00 - (1.0 / (n1 * n1)) * sum11;
	ret *= ((double)(n0 * n1)) / (n0 + n1);
	return (ret);
}

void StatsComputer::compute_spr_obs(int xi, int yi, int n, int pn, int nm1, double nm1d) {
	int a11, a12, a21, a22;
	double e11, e12, e21, e22, current_chi, current_like;
	double emin11_12, emin21_22, emin;

	a11 = double_integral[n  * pn + xi] - double_integral[(yi + 1) * pn + xi    ];
	a12 = double_integral[n  * pn + n ] - double_integral[ n       * pn + xi + 1] - double_integral[(yi + 1) * pn + n] + double_integral[(yi + 1) * pn + xi + 1];
	a21 = double_integral[yi * pn + xi];
	a22 = double_integral[yi * pn + n ] - double_integral[yi * pn + xi + 1];

	e11 =        xi  * (nm1 - yi) / nm1d;
	e12 = (nm1 - xi) * (nm1 - yi) / nm1d;
	e21 =        xi  *        yi  / nm1d;
	e22 = (nm1 - xi) *        yi  / nm1d;

	emin11_12 = min(e11, e12);
	emin21_22 = min(e21, e22);
	emin = min(emin11_12, emin21_22);

	if (emin > min_w) {
#ifdef XDP_ALLOW_DEGENERATE_PARTITIONS
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0);
#else
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22;
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		uvs_sc += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > uvs_mc)) {
		uvs_mc = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
	             + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
	             + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
	             + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0);

	uvs_sl += current_like;
	++ng_like;
	if (current_like > uvs_ml) {
		uvs_ml = current_like;
	}
}

void StatsComputer::compute_spr_all(int xi, int yi, int n, int pn, double nd) {
	int a11, a12, a21, a22;
	double e11, e12, e21, e22, current_chi, current_like;
	double emin11_12, emin21_22, emin;

	a11 = double_integral[n  * pn + xi] - double_integral[yi * pn + xi];
	a12 = double_integral[n  * pn + n ] - double_integral[ n * pn + xi] - double_integral[yi * pn + n] + double_integral[yi * pn + xi];
	a21 = double_integral[yi * pn + xi];
	a22 = double_integral[yi * pn + n ] - double_integral[yi * pn + xi];

	e11 =      xi  * (n - yi) / nd;
	e12 = (n - xi) * (n - yi) / nd;
	e21 =      xi  *      yi  / nd;
	e22 = (n - xi) *      yi  / nd;

	emin11_12 = min(e11, e12);
	emin21_22 = min(e21, e22);
	emin = min(emin11_12, emin21_22);

	if (emin > min_w) {
#ifdef XDP_ALLOW_DEGENERATE_PARTITIONS
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0);
#else
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22;
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		uvs_sc += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > uvs_mc)) {
		uvs_mc = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
	             + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
	             + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
	             + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0);

	uvs_sl += current_like;
	++ng_like;
	if (current_like > uvs_ml) {
		//Rprintf("xi: %d, yi:%d  ml:%lf \n\r",xi,yi,current_like);
		uvs_ml = current_like;
	}
}

void StatsComputer::compute_ppr_22(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int pn, int nm2, double nm2s) {
	int pij_num, Aij, Bij;
	double pij, qij, current_chi, current_like;
	double emin;

    // in the notation of the grant document:
    pij_num = (yr_hi - yr_lo - 1) * (xr_hi - xr_lo - 1);

    Aij = double_integral[yr_hi * pn + xr_hi    ] + double_integral[(yr_lo + 1) * pn + xr_lo + 1]
        - double_integral[yr_hi * pn + xr_lo + 1] - double_integral[(yr_lo + 1) * pn + xr_hi    ];
    Bij = nm2 - Aij;

    pij = pij_num / nm2s;
    qij = 1 - pij;

	emin = min(pij, qij) * nm2s;

	if (emin > min_w) {
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
		current_chi = ((Aij - nm2 * pij) * (Aij - nm2 * pij) / (nm2 * pij * (1 - pij)));
#else
		current_chi = (pij * (1 - pij) > 0) ? ((Aij - nm2 * pij) * (Aij - nm2 * pij) / (nm2 * pij * (1 - pij))) : 0;
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		uvs_sc += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > uvs_mc)) {
		uvs_mc = current_chi;
	}

	current_like = ((Aij > 0) ? (Aij * log(Aij / (pij * nm2))) : 0)
		         + ((Bij > 0) ? (Bij * log(Bij / (qij * nm2))) : 0);

	uvs_sl += current_like;
	++ng_like;
	if (current_like > uvs_ml) {
		uvs_ml = current_like;
	}
}

void StatsComputer::compute_ppr_33(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int n, int pn, double nm2) {
	int a11, a21, a31, a12, a22, a32, a13, a23, a33;
	double e11, e21, e31, e12, e22, e32, e13, e23, e33, current_chi, current_like;
	double emin1, emin2, emin3, emin4, emin;

	a11 = double_integral[(n    ) * pn + xr_lo] - double_integral[(yr_hi + 1) * pn + xr_lo    ];
	a21 = double_integral[(yr_hi) * pn + xr_lo] - double_integral[(yr_lo + 1) * pn + xr_lo    ];
	a31 = double_integral[(yr_lo) * pn + xr_lo];
	a12 = double_integral[(n    ) * pn + xr_hi] + double_integral[(yr_hi + 1) * pn + xr_lo + 1] - double_integral[(n    ) * pn + xr_lo + 1] - double_integral[(yr_hi + 1) * pn + xr_hi];
	a22 = double_integral[(yr_hi) * pn + xr_hi] + double_integral[(yr_lo + 1) * pn + xr_lo + 1] - double_integral[(yr_hi) * pn + xr_lo + 1] - double_integral[(yr_lo + 1) * pn + xr_hi];
	a32 = double_integral[(yr_lo) * pn + xr_hi] - double_integral[(yr_lo    ) * pn + xr_lo + 1];
	a13 = double_integral[(n    ) * pn + n    ] + double_integral[(yr_hi + 1) * pn + xr_hi + 1] - double_integral[(n    ) * pn + xr_hi + 1] - double_integral[(yr_hi + 1) * pn + n    ];
	a23 = double_integral[(yr_hi) * pn + n    ] + double_integral[(yr_lo + 1) * pn + xr_hi + 1] - double_integral[(yr_hi) * pn + xr_hi + 1] - double_integral[(yr_lo + 1) * pn + n    ];
	a33 = double_integral[(yr_lo) * pn + n    ] - double_integral[(yr_lo    ) * pn + xr_hi + 1];

#if 0
	int rs1, rs2, rs3, cs1, cs2, cs3;

	rs1 = a11 + a12 + a13;
	rs2 = a21 + a22 + a23;
	rs3 = a31 + a32 + a33;
	cs1 = a11 + a21 + a31;
	cs2 = a12 + a22 + a32;
	cs3 = a13 + a23 + a33;

	e11 = rs1 * cs1 / nm2;
	e21 = rs2 * cs1 / nm2;
	e31 = rs3 * cs1 / nm2;
	e12 = rs1 * cs2 / nm2;
	e22 = rs2 * cs2 / nm2;
	e32 = rs3 * cs2 / nm2;
	e13 = rs1 * cs3 / nm2;
	e23 = rs2 * cs3 / nm2;
	e33 = rs3 * cs3 / nm2;
#else
	e11 = (xr_lo            ) * (n     - 1 - yr_hi) / nm2;
	e21 = (xr_lo            ) * (yr_hi - 1 - yr_lo) / nm2;
	e31 = (xr_lo            ) * (yr_lo            ) / nm2;
	e12 = (xr_hi - 1 - xr_lo) * (n     - 1 - yr_hi) / nm2;
	e22 = (xr_hi - 1 - xr_lo) * (yr_hi - 1 - yr_lo) / nm2;
	e32 = (xr_hi - 1 - xr_lo) * (yr_lo            ) / nm2;
	e13 = (n     - 1 - xr_hi) * (n     - 1 - yr_hi) / nm2;
	e23 = (n     - 1 - xr_hi) * (yr_hi - 1 - yr_lo) / nm2;
	e33 = (n     - 1 - xr_hi) * (yr_lo            ) / nm2;
#endif

	emin1 = min(e11, e21);
	emin2 = min(e31, e12);
	emin3 = min(e22, e32);
	emin4 = min(e13, e23);
	emin1 = min(emin1, emin2);
	emin2 = min(emin3, emin4);
	emin1 = min(emin1, emin2);
	emin = min(emin1, e33);

	if (emin > min_w) {
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a31 - e31) * (a31 - e31) / e31
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22
					+ (a32 - e32) * (a32 - e32) / e32
					+ (a13 - e13) * (a13 - e13) / e13
					+ (a23 - e23) * (a23 - e23) / e23
					+ (a33 - e33) * (a33 - e33) / e33;
#else
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e31 > 0) ? ((a31 - e31) * (a31 - e31) / e31) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0)
					+ ((e32 > 0) ? ((a32 - e32) * (a32 - e32) / e32) : 0)
					+ ((e13 > 0) ? ((a13 - e13) * (a13 - e13) / e13) : 0)
					+ ((e23 > 0) ? ((a23 - e23) * (a23 - e23) / e23) : 0)
					+ ((e33 > 0) ? ((a33 - e33) * (a33 - e33) / e33) : 0);
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		uvs_sc += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > uvs_mc)) {
		uvs_mc = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
				 + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
				 + ((a31 > 0) ? (a31 * log(a31 / e31)) : 0)
				 + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
				 + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0)
				 + ((a32 > 0) ? (a32 * log(a32 / e32)) : 0)
				 + ((a13 > 0) ? (a13 * log(a13 / e13)) : 0)
				 + ((a23 > 0) ? (a23 * log(a23 / e23)) : 0)
				 + ((a33 > 0) ? (a33 * log(a33 / e33)) : 0);

	uvs_sl += current_like;
	++ng_like;
	if (current_like > uvs_ml) {
		uvs_ml = current_like;
	}
}

void StatsComputer::compute_tpr(int xl, int xm, int xh, int yl, int ym, int yh, int n, int pn, double nm3) {
	int a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44;
	double e11, e12, e13, e14, e21, e22, e23, e24, e31, e32, e33, e34, e41, e42, e43, e44;
	double current_chi, current_like;
	double emin1, emin2, emin3, emin4, emin5, emin6, emin7, emin8, emin;

    a11 = double_integral[(n ) * pn + xl] - double_integral[(yh + 1) * pn + xl    ];
    a21 = double_integral[(yh) * pn + xl] - double_integral[(ym + 1) * pn + xl    ];
    a31 = double_integral[(ym) * pn + xl] - double_integral[(yl + 1) * pn + xl    ];
    a41 = double_integral[(yl) * pn + xl];

    a12 = double_integral[(n ) * pn + xm] - double_integral[(n     ) * pn + xl + 1] - double_integral[(yh + 1) * pn + xm] + double_integral[(yh + 1) * pn + xl + 1];
    a22 = double_integral[(yh) * pn + xm] - double_integral[(yh    ) * pn + xl + 1] - double_integral[(ym + 1) * pn + xm] + double_integral[(ym + 1) * pn + xl + 1];
    a32 = double_integral[(ym) * pn + xm] - double_integral[(ym    ) * pn + xl + 1] - double_integral[(yl + 1) * pn + xm] + double_integral[(yl + 1) * pn + xl + 1];
    a42 = double_integral[(yl) * pn + xm] - double_integral[(yl    ) * pn + xl + 1];

    a13 = double_integral[(n ) * pn + xh] - double_integral[(n     ) * pn + xm + 1] - double_integral[(yh + 1) * pn + xh] + double_integral[(yh + 1) * pn + xm + 1];
    a23 = double_integral[(yh) * pn + xh] - double_integral[(yh    ) * pn + xm + 1] - double_integral[(ym + 1) * pn + xh] + double_integral[(ym + 1) * pn + xm + 1];
    a33 = double_integral[(ym) * pn + xh] - double_integral[(ym    ) * pn + xm + 1] - double_integral[(yl + 1) * pn + xh] + double_integral[(yl + 1) * pn + xm + 1];
    a43 = double_integral[(yl) * pn + xh] - double_integral[(yl    ) * pn + xm + 1];

    a14 = double_integral[(n ) * pn + n ] - double_integral[(n     ) * pn + xh + 1] - double_integral[(yh + 1) * pn + n ] + double_integral[(yh + 1) * pn + xh + 1];
    a24 = double_integral[(yh) * pn + n ] - double_integral[(yh    ) * pn + xh + 1] - double_integral[(ym + 1) * pn + n ] + double_integral[(ym + 1) * pn + xh + 1];
    a34 = double_integral[(ym) * pn + n ] - double_integral[(ym    ) * pn + xh + 1] - double_integral[(yl + 1) * pn + n ] + double_integral[(yl + 1) * pn + xh + 1];
    a44 = double_integral[(yl) * pn + n ] - double_integral[(yl    ) * pn + xh + 1];

#if 0
	int rs1, rs2, rs3, rs4, cs1, cs2, cs3, cs4;

	rs1 = a11 + a12 + a13 + a14;
	rs2 = a21 + a22 + a23 + a24;
	rs3 = a31 + a32 + a33 + a34;
	rs4 = a41 + a42 + a43 + a44;
	cs1 = a11 + a21 + a31 + a41;
	cs2 = a12 + a22 + a32 + a42;
	cs3 = a13 + a23 + a33 + a43;
	cs4 = a14 + a24 + a34 + a44;

	e11 = rs1 * cs1 / nm3;
	e21 = rs2 * cs1 / nm3;
	e31 = rs3 * cs1 / nm3;
	e41 = rs4 * cs1 / nm3;
	e12 = rs1 * cs2 / nm3;
	e22 = rs2 * cs2 / nm3;
	e32 = rs3 * cs2 / nm3;
	e42 = rs4 * cs2 / nm3;
	e13 = rs1 * cs3 / nm3;
	e23 = rs2 * cs3 / nm3;
	e33 = rs3 * cs3 / nm3;
	e43 = rs4 * cs3 / nm3;
	e14 = rs1 * cs4 / nm3;
	e24 = rs2 * cs4 / nm3;
	e34 = rs3 * cs4 / nm3;
	e44 = rs4 * cs4 / nm3;
#else
	e11 = (xl         ) * (n  - 1 - yh) / nm3;
	e21 = (xl         ) * (yh - 1 - ym) / nm3;
	e31 = (xl         ) * (ym - 1 - yl) / nm3;
	e41 = (xl         ) * (yl         ) / nm3;
	e12 = (xm - 1 - xl) * (n  - 1 - yh) / nm3;
	e22 = (xm - 1 - xl) * (yh - 1 - ym) / nm3;
	e32 = (xm - 1 - xl) * (ym - 1 - yl) / nm3;
	e42 = (xm - 1 - xl) * (yl         ) / nm3;
	e13 = (xh - 1 - xm) * (n  - 1 - yh) / nm3;
	e23 = (xh - 1 - xm) * (yh - 1 - ym) / nm3;
	e33 = (xh - 1 - xm) * (ym - 1 - yl) / nm3;
	e43 = (xh - 1 - xm) * (yl         ) / nm3;
	e14 = (n  - 1 - xh) * (n  - 1 - yh) / nm3;
	e24 = (n  - 1 - xh) * (yh - 1 - ym) / nm3;
	e34 = (n  - 1 - xh) * (ym - 1 - yl) / nm3;
	e44 = (n  - 1 - xh) * (yl         ) / nm3;
#endif

	emin1 = min(e11, e21);
	emin2 = min(e31, e41);
	emin3 = min(e12, e22);
	emin4 = min(e32, e42);
	emin5 = min(e13, e23);
	emin6 = min(e33, e43);
	emin7 = min(e14, e24);
	emin8 = min(e34, e44);
	emin1 = min(emin1, emin2);
	emin2 = min(emin3, emin4);
	emin3 = min(emin5, emin6);
	emin4 = min(emin7, emin8);
	emin1 = min(emin1, emin2);
	emin2 = min(emin3, emin4);
	emin = min(emin1, emin2);

	if (emin > min_w) {
#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a31 - e31) * (a31 - e31) / e31
					+ (a41 - e41) * (a41 - e41) / e41
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22
					+ (a32 - e32) * (a32 - e32) / e32
					+ (a42 - e42) * (a42 - e42) / e42
					+ (a13 - e13) * (a13 - e13) / e13
					+ (a23 - e23) * (a23 - e23) / e23
					+ (a33 - e33) * (a33 - e33) / e33
					+ (a43 - e43) * (a43 - e43) / e43
					+ (a14 - e14) * (a14 - e14) / e14
					+ (a24 - e24) * (a24 - e24) / e24
					+ (a34 - e34) * (a34 - e34) / e34
					+ (a44 - e44) * (a44 - e44) / e44;
#else
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e31 > 0) ? ((a31 - e31) * (a31 - e31) / e31) : 0)
					+ ((e41 > 0) ? ((a41 - e41) * (a41 - e41) / e41) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0)
					+ ((e32 > 0) ? ((a32 - e32) * (a32 - e32) / e32) : 0)
					+ ((e42 > 0) ? ((a42 - e42) * (a42 - e42) / e42) : 0)
					+ ((e13 > 0) ? ((a13 - e13) * (a13 - e13) / e13) : 0)
					+ ((e23 > 0) ? ((a23 - e23) * (a23 - e23) / e23) : 0)
					+ ((e33 > 0) ? ((a33 - e33) * (a33 - e33) / e33) : 0)
					+ ((e43 > 0) ? ((a43 - e43) * (a43 - e43) / e43) : 0)
					+ ((e14 > 0) ? ((a14 - e14) * (a14 - e14) / e14) : 0)
					+ ((e24 > 0) ? ((a24 - e24) * (a24 - e24) / e24) : 0)
					+ ((e34 > 0) ? ((a34 - e34) * (a34 - e34) / e34) : 0)
					+ ((e44 > 0) ? ((a44 - e44) * (a44 - e44) / e44) : 0);
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		uvs_sc += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > uvs_mc)) {
		uvs_mc = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
				 + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
				 + ((a31 > 0) ? (a31 * log(a31 / e31)) : 0)
				 + ((a41 > 0) ? (a41 * log(a41 / e41)) : 0)
				 + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
				 + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0)
				 + ((a32 > 0) ? (a32 * log(a32 / e32)) : 0)
				 + ((a42 > 0) ? (a42 * log(a42 / e42)) : 0)
				 + ((a13 > 0) ? (a13 * log(a13 / e13)) : 0)
				 + ((a23 > 0) ? (a23 * log(a23 / e23)) : 0)
				 + ((a33 > 0) ? (a33 * log(a33 / e33)) : 0)
				 + ((a43 > 0) ? (a43 * log(a43 / e43)) : 0)
				 + ((a14 > 0) ? (a14 * log(a14 / e14)) : 0)
				 + ((a24 > 0) ? (a24 * log(a24 / e24)) : 0)
				 + ((a34 > 0) ? (a34 * log(a34 / e34)) : 0)
				 + ((a44 > 0) ? (a44 * log(a44 / e44)) : 0);

	uvs_sl += current_like;
	++ng_like;
	if (current_like > uvs_ml) {
		uvs_ml = current_like;
	}
}

void StatsComputer::uv_ind_opt_ddp2(void) {
	uvs_ind_opt_ddp2();
}

void StatsComputer::uvs_ind_opt_ddp2(void) {
	long n = xy_nrow, src;
    long a00, a01, a10, a11;
	long b11, b12, b21, b22;
	long  nm1 = n-1;
	double nm1d = (double) (n-1);
	long xi,yi;
	double e11, e12, e21, e22, current_chi, current_like;
	double emin11_12, emin21_22, emin;
	
	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;
	long ng_chi=0;ng_like=0;
	long i=0;
	long pi = 0; // idx_perm[i];
//#define DEBUG_PRINTS
#ifdef DEBUG_PRINTS
		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (int  k = 0; k < n ; ++k) {
			cout << k << " (" << k << "): " << dx[k] << ", " << dy[k] << endl;
		}
		for (int k = 0; k < n ; ++k) {
			
			cout << k << ": " << (*sorted_dx)[i][k].first << " (" << (*sorted_dx)[i][k].second << ")" << endl;
		}
		for (int k = 0; k < n ; ++k) {
			cout << k << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

		// Use Yair's merge-sort-like implementation (assumes there are no ties)

		// NOTE: I am using the fact that the sorting of permuted y's is the same as
		// the sorting of original y's. I only have to go (a) to the permuted row instead
		// of the i'th row, and (b) pass the "sorted indices" through the permutation.

		for (long j = 0; j < n ; ++j ) {
			src = idx_perm_inv[(*sorted_dy)[pi][j].second]; // NOTE: k may be different than from the line above
			hhg_gen_y_rev[src] = j;
		}

		for (long j = 0; j < n ; ++j){
			src = (*sorted_dx)[i][j].second; // NOTE: k may be different than from the line above
			hhg_gen_xy_perm[j] = hhg_gen_y_rev[src];
			hhg_gen_source[j] = j;
			hhg_gen_inversion_count[j] = 0;
			hhg_gen_xy_perm_temp[j] = hhg_gen_xy_perm[j];
		}

		hhg_gen_inversions(hhg_gen_xy_perm_temp, hhg_gen_source, hhg_gen_inversion_count, n);

		for (long j = 0, k = 0; j < n ; ++j, ++k) {
			a00 = j - hhg_gen_inversion_count[j] ;
			a01 = hhg_gen_inversion_count[j] ;
			a10 = hhg_gen_xy_perm[j] + hhg_gen_inversion_count[j] - j ;
			a11 = n - hhg_gen_xy_perm[j] - hhg_gen_inversion_count[j] - 1 ;

			// Note that it is expected this would only be necessary for the computing the
			// observed statistic (the permutation is identity).
			// It might be a better idea to create a separate copy of this function
			// and add this only in the copy.
			
			if(true){

			    b11 =a01;
				b12 =a11;
				b21 =a00;
				b22 =a10;
				
				xi = j;
				yi = a00+a10;
				#ifdef DEBUG_PRINTS
					cout << "xi " << xi << " yi "<< yi<< endl;
				#endif
				
				#ifndef XDP_ALLOW_DEGENERATE_PARTITIONS
				// Can optimize: points with extreme ranks can be filtered beforehand...
				if (xi == 0 || xi == nm1 || yi == 0 || yi == nm1){
					continue;
				}
				#endif

				e11 =        xi  * (nm1 - yi) / nm1d;
				e12 = (nm1 - xi) * (nm1 - yi) / nm1d;
				e21 =        xi  *        yi  / nm1d;
				e22 = (nm1 - xi) *        yi  / nm1d;
				//cout<< " e11 "<<e11 << 
				emin11_12 = min(e11, e12);
				emin21_22 = min(e21, e22);
				emin = min(emin11_12, emin21_22);

				if (emin > min_w) {
				#ifdef XDP_ALLOW_DEGENERATE_PARTITIONS
					current_chi = ((e11 > 0) ? ((b11 - e11) * (b11 - e11) / e11) : 0)
								+ ((e21 > 0) ? ((b21 - e21) * (b21 - e21) / e21) : 0)
								+ ((e12 > 0) ? ((b12 - e12) * (b12 - e12) / e12) : 0)
								+ ((e22 > 0) ? ((b22 - e22) * (b22 - e22) / e22) : 0);
				#else
					current_chi = (b11 - e11) * (b11 - e11) / e11
								+ (b21 - e21) * (b21 - e21) / e21
								+ (b12 - e12) * (b12 - e12) / e12
								+ (b22 - e22) * (b22 - e22) / e22;
				#endif
				} else {
					current_chi = 0;
				}

				if (emin > w_sum){
					sum_chi += current_chi;
					++ng_chi;
				}

				if ((emin > w_max) && (current_chi > max_chi)){
					max_chi = current_chi;
				}

				current_like = ((b11 > 0) ? ((double)(b11) * log((double)b11 / e11)) : 0)
							 + ((b21 > 0) ? ((double)(b21) * log((double)b21 / e21)) : 0)
							 + ((b12 > 0) ? ((double)(b12) * log((double)b12 / e12)) : 0)
							 + ((b22 > 0) ? ((double)(b22) * log((double)b22 / e22)) : 0);

				sum_like += current_like;
				++ng_like;
				if (current_like > max_like) {
					max_like = current_like;
				}
			
			}
			
			if (store_tables) {
				//k += ((*sorted_dx)[i][k].second == i);
				int row = i * n + (*sorted_dx)[i][k].second;
				obs_tbls[        row] = a00;
				obs_tbls[  n*n + row] = a01;
				obs_tbls[2*n*n + row] = a10;
				obs_tbls[3*n*n + row] = a11;
			}

#ifdef DEBUG_PRINTS
			cout << "with point at position " << j << " of the sorted dx: ";
			cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
			if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
				cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
				//exit(1);
			}
#endif
			//accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, 1);
		}
		ng_like *=n;
		ng_chi *=n;
		sum_like /= ng_like;
		sum_chi /= ng_chi;
}

void StatsComputer::uv_ind_opt_hoeffding(void) {
	uvs_ind_opt_hoeffding();
}

void StatsComputer::uvs_ind_opt_hoeffding(void){
	long n = xy_nrow, src;
    long a00 ,  a10;
	//long a01, a11;
	//long b11, b12, b22;
	long b21;
    
	//long  nm1 = n-1;
	//double nm1d = (double) (n-1);
	long xi,yi;
//	double e11, e12, e21, e22;
	
	double hoeffding_delta;
	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;
	
	
	long i=0;
	long pi = 0; // idx_perm[i];
//#define DEBUG_PRINTS
#ifdef DEBUG_PRINTS
		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (int  k = 0; k < n ; ++k) {
			cout << k << " (" << k << "): " << dx[k] << ", " << dy[k] << endl;
		}
		for (int k = 0; k < n ; ++k) {
			
			cout << k << ": " << (*sorted_dx)[i][k].first << " (" << (*sorted_dx)[i][k].second << ")" << endl;
		}
		for (int k = 0; k < n ; ++k) {
			cout << k << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

		// Use Yair's merge-sort-like implementation (assumes there are no ties)

		// NOTE: I am using the fact that the sorting of permuted y's is the same as
		// the sorting of original y's. I only have to go (a) to the permuted row instead
		// of the i'th row, and (b) pass the "sorted indices" through the permutation.

		for (long j = 0; j < n ; ++j ) {
			src = idx_perm_inv[(*sorted_dy)[pi][j].second]; // NOTE: k may be different than from the line above
			hhg_gen_y_rev[src] = j;
		}

		for (long j = 0; j < n ; ++j){
			src = (*sorted_dx)[i][j].second; // NOTE: k may be different than from the line above
			hhg_gen_xy_perm[j] = hhg_gen_y_rev[src];
			hhg_gen_source[j] = j;
			hhg_gen_inversion_count[j] = 0;
			hhg_gen_xy_perm_temp[j] = hhg_gen_xy_perm[j];
		}

		hhg_gen_inversions(hhg_gen_xy_perm_temp, hhg_gen_source, hhg_gen_inversion_count, n);

		for (long j = 0, k = 0; j < n ; ++j, ++k) {
			a00 = j - hhg_gen_inversion_count[j] ;
			//a01 = hhg_gen_inversion_count[j] ;
			a10 = hhg_gen_xy_perm[j] + hhg_gen_inversion_count[j] - j ;
			//a11 = n - hhg_gen_xy_perm[j] - hhg_gen_inversion_count[j] - 1 ;

			// Note that it is expected this would only be necessary for the computing the
			// observed statistic (the permutation is identity).
			// It might be a better idea to create a separate copy of this function
			// and add this only in the copy.
			
			//b11 =a01;
			//b12 =a11;
			b21 =a00;
			//b22 =a10;
				
			xi = j;
			yi = a00 + a10;
			#ifdef DEBUG_PRINTS
				cout << "xi " << xi << " yi "<< yi<< endl;
			#endif
			
			//e11 =        xi  * (nm1 - yi) / nm1d;
			//e12 = (nm1 - xi) * (nm1 - yi) / nm1d;
			//e21 =        xi  *        yi  / nm1d;
			//e22 = (nm1 - xi) *        yi  / nm1d;
				
			hoeffding_delta = ((double)b21 + 1.0)/n - ((double)xi + 1.0)* ((double)yi + 1.0)/ ((double)n * (double)n);
			sum_like += hoeffding_delta * hoeffding_delta;
			

#ifdef DEBUG_PRINTS
			cout << "with point at position " << j << " of the sorted dx: ";
			cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
			if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
				cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
				//exit(1);
			}
#endif
			
		}
		//sum_like *= (double)n;
}


inline int StatsComputer::my_rand(int lo, int hi) {
	// return a random number between lo and hi inclusive.

#if 0
	// perhaps more accurate, but slow and undeterministic.
	// rand() is not a high quality PRNG anyway

    int divisor = RAND_MAX / (hi + 1);
    int retval;

    do {
        retval = my_R_rand_wrapper() / divisor; //DEBUG_RAND
    } while (retval > hi);

    return (retval + lo);
#else
    // slightly skewed but faster
    return (my_R_rand_wrapper() % (hi - lo + 1) + lo); //DEBUG_RAND
#endif
}

int StatsComputer::my_R_rand_wrapper(){
	double result;
	pthread_mutex_lock(rng_mutex);
	GetRNGstate();
	result = (int)(RAND_MAX  * unif_rand());
	PutRNGstate();
	pthread_mutex_unlock(rng_mutex);
	return result;
}

//function assumes it has access - see the other lock function
int StatsComputer::R_rand_wrapper_nolock(){
	double result;
	GetRNGstate();
	result = (int)(RAND_MAX  * unif_rand());
	PutRNGstate();
	return result;
}

void StatsComputer::R_rand_lock(){
	pthread_mutex_lock(rng_mutex);
}
void StatsComputer::R_rand_unlock(){
	pthread_mutex_unlock(rng_mutex);
}

double logfactorial(int n)
{ 
	double i=0,fact=1; 
	if(n<=1) { 
		return(0); 
	} 
	else{ 
		for(i=1;i<=n;i++){ 
			fact=fact+log(i); 
		} 
	return(fact); 
	} 
}
