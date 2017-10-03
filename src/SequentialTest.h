/*
 * SequentialTest.h
 *
 * A sequential testing procedure.
 *
 * Currently using Wald's sequential procedure, on the a specialized HHG test,
 * but maybe I'll derive and override in the future with other procedures.
 *
 */
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#ifndef SEQUENTIALTEST_H_
#define SEQUENTIALTEST_H_

#include "StatsComputer.h"

struct Compute_permutations_thread_arg; // see below

class SequentialTest : public TestIO, public ResamplingTestConfigurable {
public:
	SequentialTest(TestIO& test_io, ResamplingTestConfigurable& resampling_test_params);
	virtual ~SequentialTest();

	void run();

protected:
	void reset(void);
	bool update_sequential(int statistic_idx, bool is_null_more_extreme);
	bool update_sequential_all(double* perm_stats);

	double lA, lB;
	double exp1, exp2;

	double* llr;
	int* pvalc;
	bool* stopped_high;
	bool* stopped_low;
	int* perm_counter;

	double* local_perm_stats;

	StatsComputer** scs;
	volatile bool stop_all_flag;
	pthread_mutex_t mutex;
	pthread_mutex_t rng_mutex;
public: // this shouldn't be public but is made so due to pthread/c++ limitations
	void compute_one_permutation(int t);
	void compute_permutations(Compute_permutations_thread_arg* carg);
};

// Helper object for use with threading
struct Compute_permutations_thread_arg {
	SequentialTest* seq;
	int t; // thread number
	bool done_flag;

	Compute_permutations_thread_arg(SequentialTest* seq, int t) : seq(seq), t(t), done_flag(false) {}
};

#endif /* SEQUENTIALTEST_H_ */
