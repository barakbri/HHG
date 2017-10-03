#ifndef HHG_H_
#define HHG_H_

#include <vector>

//#define NO_THREADS
//#define DEBUG_THREADS
//#define DEBUG_CHECKS
//#define DEBUG_PRINTS // NOTE: don't use this with R_INTERFACE, and you probably should define NO_THREADS (I don't print to the R console, and I don't synchronize printing)
//#define ST_DEBUG_PRINTS
//#define DATAIN_DEBUG_PRINTS
//#define DEBUG_LIMIT_COHORT_SIZE
//#define DEBUG_LIMIT_NR_RANGES

#undef ERROR
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

typedef enum {
	//
	// Univariate tests
	//

	// Goodness of fit tests
	UV_GOF_WXN 			= 0, 	// Wilcoxon
	UV_GOF_AD			= 1,	// Anderson-Darling
	UV_GOF_CVM_KS		= 2,	// Cramer-von Mises and Kolmogorov-Smirnov
	UV_GOF_DCOV			= 3,	// Distance covariance
	UV_GOF_XDP2			= 4,	// ADP (also DDP) with K = 2
	UV_GOF_XDP3			= 5,	// ADP (also DDP) with K = 3
	UV_GOF_XDP			= 6,	// ADP (also DDP) for arbitrary K (but no max variant, and faster for K >= 4)

	// K-sample tests
	UV_KS_KW			= 7,	// Kruskal-Wallis
	UV_KS_AD			= 8,	// Anderson-Darling
	UV_KS_CVM_KS		= 9,	// Cramer-von Mises and Kolmogorov-Smirnov
	UV_KS_DCOV			= 10,	// Distance covariance
	UV_KS_DS			= 11,	// Dynamic slicing
	UV_KS_XDP2			= 12,	// ADP (also DDP) with K = 2
	UV_KS_XDP3			= 13,	// ADP (also DDP) with K = 3
	UV_KS_XDP			= 14,	// ADP (also DDP) with K = 2

	// Independence tests
	UV_IND_AD			= 15,	// Anderson-Darling
	UV_IND_CVM_KS		= 16,	// Cramer-von Mises and Kolmogorov-Smirnov
	UV_IND_DCOV			= 17,	// Distance covariance
	UV_IND_DDP2			= 18,	// DDP with K = 2
	UV_IND_DDP3_C		= 19,	// DDP with K = 3 (using center cell only)
	UV_IND_DDP3			= 20,	// DDP with K = 3
	UV_IND_DDP4			= 21,	// DDP with K = 4
	UV_IND_DDP			= 22,	// DDP with arbitrary K (but no max variant, and faster for K >= 4)
	UV_IND_ADP2			= 23,	// ADP with K = 2
	UV_IND_ADP3_C		= 24,	// ADP with K = 3 (using center cell only)
	UV_IND_ADP3			= 25,	// ADP with K = 3
	UV_IND_ADP4			= 26,	// ADP with K = 4
	UV_IND_ADP			= 27,	// ADP with arbitrary K (but no max variant, and faster for K >= 3)

	//
	// Multivariate tests
	//

	// TODO MV_GOF tests

	// Two and K-sample
	MV_TS_HHG			= 28,	// Two-sample HHG
	MV_KS_HHG			= 29,	// K-sample HHG
	MV_KS_HHG_EXTENDED	= 30,	// K-sample HHG extended

	// Independence
	MV_IND_HHG_NO_TIES	= 31,	// Independence HHG without ties
	MV_IND_HHG			= 32,	// Independence HHG
	MV_IND_HHG_EXTENDED = 33,	// Independence HHG extended
	
	//
	// Extras and experiments
	//

	UV_GOF_EXISTING		= 34,
	MV_TS_EXISTING		= 35,
	MV_KS_EXISTING		= 36,
	MV_IND_EXISTING		= 37,

	CI_UVZ_NN			= 38,
	CI_UVZ_GAUSSIAN		= 39,
	CI_MVZ_NN			= 40,
	CI_MVZ_GAUSSIAN		= 41,
	CI_UDF_ADP_MVZ_NN	= 42,
	CI_MVZ_NN_GRID_BW	= 43,

	//Multi Partition Versions
	UV_KS_MDS			= 44,
	UV_KS_XDP_MK        = 45,
	UV_IND_ADP_MK       = 46,
	
	//Optimized Versions of DDP2 and hoeffding
	UV_IND_OPT_DDP2		= 47,
	UV_IND_OPT_HOEFFDING= 48,
} ScoreType;

// These are univariate *distribution-free* independence tests
#define IS_UV_DF_IND_TEST(tt) ( \
		(tt) == UV_IND_DDP2       || (tt) == UV_IND_ADP2     || (tt) == UV_IND_DDP3_C   || \
		(tt) == UV_IND_ADP3_C     || (tt) == UV_IND_DDP3     || (tt) == UV_IND_ADP3     || \
		(tt) == UV_IND_DDP4       || (tt) == UV_IND_ADP4     || (tt) == UV_IND_DDP      || \
		(tt) == UV_IND_ADP        || (tt) == UV_IND_ADP_MK   )

// These are univariate *distribution-free* K-sample tests
#define IS_UV_DF_KS_TEST(tt) ((tt) == UV_KS_XDP2 || (tt) == UV_KS_XDP3 || (tt) == UV_KS_XDP || (tt) == UV_KS_XDP_MK)

// These are univariate *distribution-free* goodness-of-fit tests
#define IS_UV_DF_GOF_TEST(tt) ((tt) == UV_GOF_XDP2 || (tt) == UV_GOF_XDP3 || (tt) == UV_GOF_XDP)

// These are CI tests with MV Z
#define IS_CI_MVZ_TEST(tt) ((tt) == CI_MVZ_NN         || (tt) == CI_MVZ_GAUSSIAN || \
		                    (tt) == CI_UDF_ADP_MVZ_NN || (tt) == CI_MVZ_NN_GRID_BW)

#define IS_UV_GOF_XDP_TEST IS_UV_DF_GOF_TEST
#define IS_UV_KS_XDP_TEST IS_UV_DF_KS_TEST

// These are all the 2- and k-sample tests
#define IS_KS_TEST(tt) ( \
		(tt) == MV_TS_HHG || (tt) == MV_KS_HHG || (tt) == MV_KS_HHG_EXTENDED || IS_UV_DF_KS_TEST(tt) || \
        (tt) == UV_KS_KW  || (tt) == UV_KS_AD  || (tt) == UV_KS_CVM_KS       || (tt) == UV_KS_DCOV   || \
        (tt) == UV_KS_DS  || (tt) == UV_KS_MDS || (tt) == MV_TS_EXISTING || (tt) == UV_KS_XDP_MK)

#define IS_UVY_TEST(tt) (IS_KS_TEST(tt) || IS_UV_DF_IND_TEST(tt) || (tt) == CI_UDF_ADP_MVZ_NN)

#define BASE_NR_STATS 4

//#define XDP_ALLOW_DEGENERATE_PARTITIONS
#define XDP_NORMALIZE

typedef std::vector< std::vector<double> > matrix;
typedef std::pair<double, int> dbl_int_pair;
typedef std::vector<dbl_int_pair> dbl_int_pair_vector;
typedef std::vector< std::vector<dbl_int_pair> > dbl_int_pair_matrix;

// The following is actually a conglomeration or flattening of what should have been a
// Configurable class hierarchy... for now I'm settling for this.

struct ScoreConfigurable {
	ScoreConfigurable(ScoreType score_type, double w_sum, double w_max, double* score_params_r);
	~ScoreConfigurable();

	ScoreType score_type;
	double w_sum;
	double w_max;
	int K; // used for HHGCI NN kernel width, and for order of ADP/DDP partitions
	bool correct_mi_bias;
	double sig;
	double lambda; //lambda and Mk_Maxk are three parameters used for the ds, mds functions.
	int equipartition_type; // type of equipartition to perform (for MDS, 0 is none ,1 is no clumps, 2 is placemints in m bins) (for XDP KS, 0 for none , 1 for  m bins)
	int equipartition_nr_cells_m;
	int Mk_Maxk;
	int adp_mk_tables_nr;
	int* adp_mk_tables_m;
	int* adp_mk_tables_l;
	int prior_length; // prior to be used for MDS 
	double* prior; 
	int nnh; // NN kernel width used for our statistic in the CI test
	int nnh_lsb; // NN kernel width used for locally smoothed bootstrap (for computing p-values in the CI test)
	int nnh_grid_cnt;
	double* nnh_grid;
	ScoreType uv_score_type; // a univariate score to use in the extended HHG test
	int nr_stats;

protected:
	void parse_params(ScoreType st, double* extra_params_r);
};

struct ResamplingTestConfigurable : public ScoreConfigurable {
	ResamplingTestConfigurable(ScoreType score_type, double w_sum, double w_max, double* score_params_r,
			int nr_perm, bool is_sequential, double alpha, double alpha0, double beta0, double eps, int base_seed, int nr_threads);

	int nr_perm;
	bool is_sequential;
	double alpha, alpha0, beta0, eps;
	int base_seed;
	int nr_threads;

	int nr_perm_per_thread;

protected:
	int get_available_nr_threads(void);
};

struct TestIO {
	TestIO(int xy_nrow, int y_ncol, double* dx, double* dy, double* y, bool tables_wanted, bool perm_stats_wanted, ResamplingTestConfigurable& score_params);
	~TestIO();

	void release(void);

	//
	// Inputs
	//

	int xy_nrow, y_ncol;

	// These are raw inputs
	double *dx, *dy, *y;
	double *z, *dz, *null_dist;

	// These inputs are created by preprocessing
	dbl_int_pair_matrix *sorted_dx, *sorted_dy, *sorted_dz;
	int *ranked_dx, *ranked_dy;
	int nr_groups, *y_counts;
	double *adp, *adp_l, *adp_r;
	double *adp_mk, *adp_l_mk, *adp_r_mk; 
	
	
	int is_equipartition;
	int equipartition_m_nr_bins;
	int choose_nr_atoms;
	//
	// Outputs
	//

	bool tables_wanted;
	bool perm_stats_wanted;
	bool k_stats_wanted; //for MDS, max by k data , or for ks_xdp_mk SC and SL by k nr. partitions
	bool debug_vec_wanted;

	SEXP R_output;
	double *obs_stats, *obs_tbls, *pvals, *perm_stats,*k_stats,*debug_vec;
	const static int DEBUG_VEC_SIZE=10000; //could be changed later;
	protected:
	void allocate_outputs(ResamplingTestConfigurable& score_params);
	void preprocess(ResamplingTestConfigurable& score_params);
	void count_unique_y(void);
	void compute_distances(void);
	void sort_x_distances_per_row(void);
	void sort_y_distances_per_row(void);
	void sort_z_distances_per_row(void);
	void rank_x_distances_per_row(void);
	void rank_y_distances_per_row(void);
	
	void sort_x_distances_opt(void);
	void sort_y_distances_opt(void);
	
	void declare_adp_independence(int n,int K);
	void declare_adp_independence_mk(int n,int K);
	void declare_adp_k_sample(int n,int K);
	void declare_adp_k_sample_mk(int n,int K);
	
	void compute_adp_independence(int n, int K);
	void compute_adp_independence_mk_single(int n, int K);
	void compute_adp_k_sample(int n, int M);
	void compute_adp_k_sample_mk(int n, int M);
	void compute_adp_independence_mk(int n, int M);

	double my_choose(int n, int k);
	double my_lchoose(int n, int k);

	static inline bool dbl_int_pair_comparator(const dbl_int_pair& l, const dbl_int_pair& r) { return l.first < r.first; }
	static inline bool int_comparator(const int& l, const int& r) { return l < r; }
};

#endif /* HHG_H_ */
