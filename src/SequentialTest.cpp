/*
 * SequentialTest.cpp
 *
 */

#ifdef _WIN32
#include <windows.h>
#include <stddef.h>
#else
#include <unistd.h>
#endif

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

#include "SequentialTest.h"

#ifdef GWAS_INTERFACE
static char* get_strtime(void) { time_t now = time(0); char* strtime = asctime(localtime(&now)); strtime[strlen(strtime) - 1] = '\0'; return strtime; }
#endif

using namespace std;

// Parallelization helper functions
// ================================================================================================

#ifndef NO_THREADS

extern "C" {

static void* compute_permutations_thread(void* arg) {
	Compute_permutations_thread_arg* carg = reinterpret_cast<Compute_permutations_thread_arg*>(arg);

	(carg->seq)->compute_permutations(carg); // runs permutations and updates the sequential test state through a mutex, eventually stopping early or after all permutations are completed

	carg->done_flag = true;
	return NULL; // implicitly calls pthread_exit()
}

} // extern "C"

#endif

// ================================================================================================

SequentialTest::SequentialTest(TestIO& test_io, ResamplingTestConfigurable& resampling_test_params) :
		TestIO(test_io), ResamplingTestConfigurable(resampling_test_params)
{
#ifndef _WIN32
	
#endif

	double p0 = alpha * (1 + eps);
	double p1 = alpha * (1 - eps);

	exp1 = log((p1 / (1 - p1)) / (p0 / (1 - p0)));
	exp2 = log((1 - p1) / (1 - p0));

	lA = log((1 - beta0) / alpha0);
	lB = log(beta0 / (1 - alpha0));

	llr          = new double[nr_stats];
	pvalc        = new int   [nr_stats];
	stopped_high = new bool  [nr_stats];
	stopped_low  = new bool  [nr_stats];
	perm_counter = new int   [nr_stats];

	reset();

	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&rng_mutex, NULL);

	scs = new StatsComputer*[nr_threads];
	for (int t = 0; t < nr_threads; ++t) {
		scs[t] = new StatsComputer(test_io, *this, &rng_mutex);
	}

	if (perm_stats_wanted) {
		local_perm_stats = perm_stats; // we will move this along the output vector
	} else {
		local_perm_stats = new double[nr_stats];
	}
}

SequentialTest::~SequentialTest() {
	delete[] llr;
	delete[] pvalc;
	delete[] stopped_high;
	delete[] stopped_low;
	delete[] perm_counter;

	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&rng_mutex);

	for (int t = 0; t < nr_threads; ++t) {
		delete scs[t];
	}
	delete[] scs;

	if (!perm_stats_wanted) {
		delete[] local_perm_stats;
	}
}

void SequentialTest::reset(void) {
	for (int i = 0; i < nr_stats; ++i) {
		llr[i] = 0;
		pvalc[i] = 1;
		stopped_high[i] = false;
		stopped_low[i] = false;
		perm_counter[i] = 0;
	}
}

void SequentialTest::run(void) {
	reset();

	// Compute observed statistics
	scs[0]->compute();
	scs[0]->get_stats(obs_stats);

#ifdef ST_DEBUG_PRINTS
	Rprintf("Observed stats: %g, %g, %g, %g\n", obs_stats[0], obs_stats[1], obs_stats[2], obs_stats[3]);
#endif

	if (nr_perm > 0) {
		// Compute many null statistics and estimate the p-value
		stop_all_flag = false;

#ifdef NO_THREADS
		for (int k = 0; k < nr_perm; ++k) {
			scs[0]->permute_and_compute();
			scs[0]->get_stats(local_perm_stats);

#ifdef ST_DEBUG_PRINTS
			cout << "Perm " << k << " stats:";
			for (int i = 0; i < nr_stats; ++i) {
				cout << " " << local_perm_stats[i];
			}
			cout << endl;
#endif

			// count more extreme values of the test statistics than that observed
			bool stp_all = update_sequential_all(local_perm_stats);

			if (perm_stats_wanted) {
				local_perm_stats += nr_stats;
			}

			if (stp_all) {
				break;
			}

#ifdef GWAS_INTERFACE
			int k_progress = nr_perm / 100 + 1;
			if ((k % k_progress) == 0) {
				cout << get_strtime() << ": permutation test " << round((100.0 * k) / (nr_perm + 1)) << "% complete" << endl;
			}
#endif
		}

#else // => use threads

		pthread_t* worker_threads = new pthread_t[nr_threads];

		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED); // due to an extremely annoying bug, I am not using joinable threads
		pthread_attr_setstacksize(&attr, 4 * 1024 * 1024); // if we have a maximal feasible xy_nrow = 10000, then this should leave plenty of room for recursive calls {1,2,...,log2(xy_nrow)}

		Compute_permutations_thread_arg** worker_args = new Compute_permutations_thread_arg*[nr_threads];

		for (int t = 0; t < nr_threads; ++t) {
			worker_args[t] = new Compute_permutations_thread_arg(this, t);
			pthread_create(&(worker_threads[t]), &attr, compute_permutations_thread, (void*)(worker_args[t]));
		}

		// Manually join the worker threads (might be unnecessarily slow for small datasets)
		bool all_done = false;
		while (!all_done) {
			all_done = true;
			for (int t = 0; t < nr_threads; ++t) {
				all_done &= worker_args[t]->done_flag;
			}

#ifdef _WIN32
			Sleep(100);
#else
			usleep(100);
#endif
		}

		for (int t = 0; t < nr_threads; ++t) {
			delete worker_args[t];
		}

		delete[] worker_args;
		pthread_attr_destroy(&attr);
		delete[] worker_threads;

#endif // NO_THREADS

#ifdef GWAS_INTERFACE
		if (stop_all_flag) {
			cout << get_strtime() << ": Sequential testing says stop on all statistics. Moving on." << endl;
		} else {
			cout << get_strtime() << ": Computed all permutations!" << endl;
		}
#endif

		// Compute p-values
		for (int i = 0; i < nr_stats; ++i) {
			pvals[i] = ((double)(pvalc[i])) / (nr_perm + 1);
		}
	} else {
		for (int i = 0; i < nr_stats; ++i) {
			pvals[i] = NA_REAL;
		}
	}
}

// NOTE: This function runs in the context of a worker thread
void SequentialTest::compute_permutations(Compute_permutations_thread_arg* carg) {
	int t = carg->t;

	// rand() is not thread safe, and its behavior varies markedly between Linux (process state, usually) and Windows (thread state, usually)
	// The following is not a proper solution, but it'll have to do for now.
#ifdef _WIN32
	
#else
	// Should use drand48_r() on Linux. For now it's not critical since I don't really use multiple
	// concurrent threads when I'm working on the cluster/cloud.
#endif

	for (int k = 0; k < nr_perm_per_thread; ++k) {
		scs[t]->permute_and_compute();

		pthread_mutex_lock(&mutex);

		scs[t]->get_stats(local_perm_stats);
		

#ifdef ST_DEBUG_PRINTS
		if (t == 0) {
			cout << "Perm " << k << " stats:";
			for (int i = 0; i < nr_stats; ++i) {
				cout << " " << local_perm_stats[i];
			}
			cout << endl;
		}
#endif

		// count more extreme values of the test statistics than that observed
		bool stp_all = update_sequential_all(local_perm_stats);

		if (perm_stats_wanted) {
			local_perm_stats += nr_stats;
		}

		if (stp_all) {
			pthread_mutex_unlock(&mutex);
			break;
		}

#ifdef GWAS_INTERFACE
		int k_progress = nr_perm_per_thread / 100 + 1;
		if ((t == 0) && (k % k_progress == 0)) {
			cout << get_strtime() << ": permutation test " << round((100.0 * k) / nr_perm_per_thread) << "% complete" << endl;
		}
#endif

		pthread_mutex_unlock(&mutex);
	}
}

bool SequentialTest::update_sequential(int statistic_idx, bool is_null_more_extreme) {
	if (!is_sequential) {
		pvalc[statistic_idx] += is_null_more_extreme;
		++perm_counter[statistic_idx];
		return false;
	}

	if (stopped_high[statistic_idx]) {
		return true;
	}

	pvalc[statistic_idx] += is_null_more_extreme;
	llr[statistic_idx] += is_null_more_extreme * exp1 + exp2;
	++perm_counter[statistic_idx];

	if ((!stopped_low[statistic_idx]) && (llr[statistic_idx] <= lB)) {
		// mark as useless and stop permuting
		pvalc[statistic_idx] = nr_perm + 1;
		stopped_high[statistic_idx] = true;
		return true;
	} else if (llr[statistic_idx] >= lA) {
		// we can probably stop permuting and conclude x is truly
		// associated with y. But we want the full p-value
		// in this case so that we can follow up with a data driven
		// multiplicity adjustment such as BH. So we have to just make
		// sure that we wont stop in the other direction (unlikely but
		// possible).
		stopped_low[statistic_idx] = true;
	}

	return false;
}

bool SequentialTest::update_sequential_all(double* perm_stats) {
	bool stp_all = true;

	for (int i = 0; i < nr_stats; ++i) {
		stp_all &= update_sequential(i, perm_stats[i] >= obs_stats[i]);
	}

	stop_all_flag = stp_all;

	return (stp_all);
}
