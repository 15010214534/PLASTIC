/* PLASTIC v0.1
/* Header file for the integrate and fire network program
 * Copyright Guy Billings 2007
 * Contact: g.billings@ucl.ac.uk
 */
/*This file is part of PLASTIC.

    PLASTIC is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PLASTIC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PLASTIC; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
 
#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define __GSL_RANGE_CHECK_OFF__

extern gsl_rng *rgen;

struct double_list {
	double element;
	struct double_list *next_element;
	};
typedef struct double_list LIST;
typedef LIST* NEXT_LIST_UNIT;
 
void assign_value(int, int, int, double *, double);
double retrieve_value(int, int, int, double*);
NEXT_LIST_UNIT double_push_back(double, NEXT_LIST_UNIT);
void list_to_array(NEXT_LIST_UNIT,int,double*);
void error_out(int);
void deallocate_list(NEXT_LIST_UNIT,int);
void copy_gsl_mat(gsl_matrix *, double *, int , int);
void import_gsl_matrix(gsl_matrix *out, double *in, int k, int l);
void copy_double_arr(double *, double *, int);
double * build_double_array(int);
void input_activity_random(int,double *);
double exp_dist(double);
double custom_rand(double, gsl_rng *);
void scal_mult_arr(double,int, double *);
void scal_add_arr(double,int, double *);
double gaussrand(double);
void  update_V(int, double, double, double*, double*, double*, double*, double, double, double);
gsl_matrix* build_gsl_matrix(int, int, double, int);
void update_weights(int,int,int, gsl_matrix*,gsl_matrix*, gsl_matrix*, gsl_matrix*, double, double*,double,int,int,gsl_matrix*,double,gsl_matrix*,double,double,double,double,gsl_matrix*,gsl_matrix*,double,double); 
void update_weights_dep(int,int,int, gsl_matrix*,gsl_matrix*, gsl_matrix*, gsl_matrix*, double, double*,double,double,int,int,gsl_matrix*,gsl_matrix*,double,double,double,double,double,double,double,gsl_matrix*,gsl_matrix*);
void syn_input(int,int,gsl_matrix*, gsl_matrix *,gsl_matrix*, gsl_matrix *, double*, double*, double*,double,double,double);
void total_input(int, int, double*, gsl_matrix*, gsl_matrix*);
void mixed_rate(int inputs, double *activity, double *rates, double cor_rate, double *flag, double base, double var_rate);
void inhm_rate(int inputs, double *activity, double *rates, double *target_rates, double *rate_vars);
void drive_noise(int, double, double, double, double *);
void localised_stimulus(int, double *, double, double, double, int, double);
void pre_fire(int, double, double*, double*,double *);
void post_fire(int, double, double, double*, double*);
int the_terminator(int, int,int,int,int);
double overlap_vec(int,int ,gsl_matrix*, gsl_matrix*, gsl_matrix*);
double sum_arr(int, double*);
void any_of(int, double*, int*);
void decay(int, int, gsl_matrix*, gsl_matrix*, gsl_matrix*, 
gsl_matrix*,double,double,double,double,gsl_matrix *,gsl_matrix*,int,double,double,double,int,double,int,gsl_matrix*,int);
void record_spikes(int*,int,int,double*,NEXT_LIST_UNIT*,NEXT_LIST_UNIT*,NEXT_LIST_UNIT*,int*,double, int);
void average_weights(int, int, double*, gsl_matrix *, gsl_matrix *, NEXT_LIST_UNIT *,double,NEXT_LIST_UNIT *);
void accum_double_array(double*, double*, int);
void find_frequency(int,double*,double*,NEXT_LIST_UNIT*,double);
void find_connections(int, int, gsl_matrix*, double *);
int test_stimloc(int,int);
void euclid(int,int,gsl_matrix*,gsl_matrix*,gsl_matrix*,NEXT_LIST_UNIT*);
void pearson(int,int,gsl_matrix*,gsl_matrix*,gsl_matrix*,NEXT_LIST_UNIT*);
void mult_matrix(gsl_matrix *, gsl_matrix *, gsl_matrix *);
int find_recurrents(int,gsl_matrix*);
void track_weight(int,gsl_matrix*,gsl_matrix*,NEXT_LIST_UNIT*);
void total_current(int, int, gsl_matrix*, gsl_matrix*, double *);
void out_current(double*, double, double*, NEXT_LIST_UNIT*, int timesteps, int done);
void snapshot(gsl_matrix*,double*,int,int,int,int);
