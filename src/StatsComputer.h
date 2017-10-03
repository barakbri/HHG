/*
 * StatsComputer.h
 *
 */

#ifndef STATSCOMPUTER_H_
#define STATSCOMPUTER_H_

#include "HHG.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

class StatsComputer : public TestIO, public ScoreConfigurable {
public:
	StatsComputer(TestIO& test_input, ScoreConfigurable& score_params, pthread_mutex_t *rng_mutex_param);
	virtual ~StatsComputer();

	void compute(void);
	void permute_and_compute(void);
	void get_stats(double* stats);

protected:
	void resample_univariate(void);
	void resample_multivariate(void);
	void resample_uvz_ci(void);
	void resample_mvz_ci(void);
	void resample_dummy(void);

	void uv_gof_wxn(void);
	void uv_gof_ad(void);
	void uv_gof_cvm_ks(void);
	void uv_gof_dcov(void);
	void uv_gof_xdp2(void);
	void uv_gof_xdp3(void);
	void uv_gof_xdp(void);

	void uv_ks_kw(void);
	void uv_ks_ad(void);
	void uv_ks_cvm_ks(void);
	void uv_ks_dcov(void);
	void uv_ks_ds(void);
	void uv_ks_mds(void);
	void uv_ks_xdp2(void);
	void uv_ks_xdp3(void);
	void uv_ks_xdp(void);
	void uv_ks_xdp_mk(void);
	

	void uv_ind_ad(void);
	void uv_ind_cvm_ks(void);
	void uv_ind_dcov(void);
	void uv_ind_ddp2(void);
	void uv_ind_ddp3_c(void);
	void uv_ind_ddp3(void);
	void uv_ind_ddp4(void);
	void uv_ind_ddp(void);
	void uv_ind_adp2(void);
	void uv_ind_adp3_c(void);
	void uv_ind_adp3(void);
	void uv_ind_adp4(void);
	void uv_ind_adp(void);
	void uv_ind_adp_mk(void);
	
	void uv_ind_opt_ddp2(void);
	void uv_ind_opt_hoeffding(void);

	void mv_ts_hhg(void);
	void mv_ks_hhg(void);
	void mv_ks_hhg_extended(void);

	void mv_ind_hhg_no_ties(void);
	void mv_ind_hhg(void);
	void mv_ind_hhg_extended(void);

	void uv_gof_existing(void);
	void mv_ts_existing(void);
	void mv_ks_existing(void);
	void mv_ind_existing(void);

	void ci_uvz_nn(void);
	void ci_uvz_gaussian(void);
	void ci_mvz_nn(void);
	void ci_mvz_gaussian(void);
	void ci_udf_adp_mvz_nn(void);
	void ci_mvz_nn_grid(void);

	double compute_ht(void);
	double compute_edist(void);

	void sort_xy_distances_per_row(void);

	void accumulate_2x2_contingency_table(double a00, double a01, double a10, double a11, double nrmlz, double reps);
	void accumulate_local_stats(double chi, double like, double emin);
	void hhg_gen_inversions(int *permutation, int *source, int *inversion_count, int dim);
	void hhg_gen_merge(int *permutation, int *source, int *inversion_count, int dim);
	void compute_ordered_ranks(int n, double* xx, int* yy);
	void compute_single_integral(int n, double* xx, int* yy);
	void compute_double_integral(int n, double* xx, int* yy);
	void compute_double_integral_eqp(long n, double* xx, int* yy,long nr_atoms);
	int count_sample_points_in_rect(int xl, int xh, int yl, int yh);
	long count_sample_points_in_rect_eqp(long xl, long xh, long yl, long yh);
	double count_ddp_with_given_cell(int xl, int xh, int yl, int yh);
	double count_adp_with_given_cell(int xl, int xh, int yl, int yh);
	int compute_adp_mk_cell_type(long xl, long xh, long yl, long yh,long nr_atoms);
	double count_adp_mk_cell_type(int M,int L,int type, long w, long h,long nr_atoms);

	int my_rand(int lo, int hi);
	int my_R_rand_wrapper();
	int R_rand_wrapper_nolock();
	void R_rand_lock();
	void R_rand_unlock();
	
	void compute_spr_obs(int xi, int yi, int n, int pn, int nm1, double nm1d);
	void compute_spr_all(int xi, int yi, int n, int pn, double nd);
	void compute_ppr_22(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int pn, int nm2, double nm2s);
	void compute_ppr_33(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int n, int pn, double nm2);
	void compute_tpr(int xl, int xm, int xh, int yl, int ym, int yh, int n, int pn, double nm3);

	void uvs_gof_wxn(void);
	void uvs_gof_ad(void);
	void uvs_gof_cvm_ks(void);
	void uvs_gof_dcov(void);
	void uvs_gof_xdp2(void);
	void uvs_gof_xdp3(void);
	void uvs_gof_xdp(void);

	void uvs_ks_kw(void);
	void uvs_ks_ad(void);
	void uvs_ks_cvm_ks(void);
	void uvs_ks_dcov(void);
	void uvs_ks_ds(void);
	void uvs_ks_mds(void);
	void uvs_ks_xdp2(void);
	void uvs_ks_xdp3(void);
	void uvs_ks_xdp(void);
	void uvs_ks_xdp_mk(void); //DEBUG_LINUX

	void uvs_ind_ad(void);
	void uvs_ind_cvm_ks(void);
	void uvs_ind_dcov(void);
	void uvs_ind_ddp2(void);
	void uvs_ind_ddp3_c(void);
	void uvs_ind_ddp3(void);
	void uvs_ind_ddp4(void);
	void uvs_ind_ddp(void);
	void uvs_ind_adp2(void);
	void uvs_ind_adp3_c(void);
	void uvs_ind_adp3(void);
	void uvs_ind_adp4(void);
	void uvs_ind_adp(void);
	void uvs_ind_adp_mk(void);
	
	void uvs_ind_opt_ddp2(void);
	void uvs_ind_opt_hoeffding(void);
	
	void (StatsComputer::*compute_score)(void);
	void (StatsComputer::*resample)(void);
	void (StatsComputer::*hhg_extended_uvs)(void);

protected:
	double min_w;
	bool should_randomize;
	bool store_tables;

    double sum_chi, sum_like, max_chi, max_like;
    double *sum_chi_grid, *sum_like_grid, *max_chi_grid, *max_like_grid;
	double *mds_max_chi_by_k; 
	double *mds_max_loglikelihood_by_k; 
	double *xdp_sc_mk, *xdp_sl_mk; //DEBUG_LINUX
	double *adp_ind_sc_mk, *adp_ind_sl_mk; //DEBUG_LINUX
    double max_sum_chi, max_sum_like, sum_max_chi, sum_max_like;

    int *y_perm;
	int *y0_idx, *y1_idx; // indices of samples with y_i == 0 and 1 respectively
	int *idx_1_to_n;
	int *idx_perm, *idx_perm_inv;
	
	pthread_mutex_t *rng_mutex;
	
protected:
	// FIXME one day the univariate score, and scores in general, will be objects
	// For now the fact that some fields are shared and some are supposed to be
	// exclusive to univariate scores (rather than the external multivariate score)
	// is not clear. So I'm trying to gather all these fields in the following block:

	int uvs_n;
	double *uvs_x, *uvs_y;
	double *uvs_xr; // this is also actually integer, but for backward compatibility with higher levels I keep it double
	int *uvs_yr;
	double uvs_sc, uvs_mc, uvs_sl, uvs_ml;
	int* uvs_yc;
	double uvs_y0;

	double kahan_c_chi, kahan_c_like;
	int ng_chi, ng_like;
	int *x_ordered_by_y, *y_ordered_by_x;
	double *tbl_o, *tbl_e;
	int* double_integral;
	int* double_integral_eqp;
	int dintegral_zero_based_idxs;
	int dintegral_pn;
	int dintegral_pn_eqp;
	double *kw_rs;
	int **ds_ctab, *ds_idx,**ds_ctab_bins;
	double *ds_score, *ds_counts , *ds_score_pearson; //originally ds_score was the only score for the likelihood as implemented by the dynamic slicing paper, we later added a pearson variant.

protected:
	int *hhg_gen_inversion_count, *hhg_gen_source, *hhg_gen_xy_perm, *hhg_gen_xy_perm_temp, *hhg_gen_y_rev;
	int *hhg_gen_left_buffer, *hhg_gen_right_buffer, *hhg_gen_left_source_buffer, *hhg_gen_right_source_buffer;

	struct dbl_dbl_int {
		double x;
		double y;
		int i;
	};

	typedef std::vector< std::vector<dbl_dbl_int> > dbl_dbl_int_matrix;

	static inline bool dbl_int_pair_comparator(const dbl_int_pair& l, const dbl_int_pair& r) {
		return l.first < r.first;
	}

	static inline bool dbl_dbl_int_pair_comparator_xy(const dbl_dbl_int& l, const dbl_dbl_int& r) {
		return ((l.x < r.x) || ((l.x == r.x) && (l.y > r.y)));
	}

	dbl_dbl_int_matrix sorted_dx_gen;

	dbl_int_pair_vector nn_sorted_x, nn_sorted_y;
};

double logfactorial(int n);


#endif /* STATSCOMPUTER_H_ */
