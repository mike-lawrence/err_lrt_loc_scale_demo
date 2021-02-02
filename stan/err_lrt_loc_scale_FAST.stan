data{

	// num_lrt: number of trials
	int<lower=1> num_lrt ;

	// lrt: lrt on each trial
	vector[num_lrt] lrt ;

	// num_err: number of trials
	int<lower=1> num_err ;

	// err: err on each trial
	int<lower=0,upper=1> err[num_err] ;

	// num_err_cells: number of error summary cells
	int<lower=1> num_err_cells ;

	// err_num_obs: number of error observations per cell of the summary
	int err_num_obs[num_err_cells] ;

	// err_num_err: number of errors per cell of the summary
	int err_num_err[num_err_cells] ;

	// num_subj: number of subj
	int<lower=1> num_subj ;

	// num_rows_uW: num rows in uW
	int<lower=1> num_rows_uW ;

	// num_cols_uW: num cols in uW
	int<lower=1> num_cols_uW ;

	// uW: unique entries in the within predictor matrix
	matrix[num_rows_uW,num_cols_uW] uW ;

	// lrt_index: index of each lrt in flattened subject-by-condition value matrix
	int lrt_index[num_lrt] ;

	// err_index: index of each lrt in flattened subject-by-condition value matrix
	int err_index[num_err] ;

	// err_cell_index: index of each err summary in flattened subject-by-condition value matrix
	int err_cell_index[num_err_cells] ;

}
transformed data{

	// lrt_mean: mean lrt value
	real lrt_mean = mean(lrt) ;

	// lrt_sd: sd of lrts
	real lrt_sd = sd(lrt) ;

	// lrt_: lrtervations scaled to have zero mean and unit variance
	vector[num_lrt] lrt_ = (lrt-lrt_mean)/lrt_sd ;

	// compute observed intercept for the error data
	real obs_err_intercept = logit(mean(to_vector(err_num_err)./to_vector(err_num_obs))) ;

	// num_coef: lrt_loc, lrt_scale, err_logodds
	int num_coef = num_cols_uW*3 ;

}
parameters{

	//for parameters below, trailing underscore denotes that they need to be un-scaled in generated quantities

	// coef_mean_: mean (across subj) for each coefficient
	row_vector[num_coef] mean_coef_ ;

	// coef_sd_: sd (across subj) for each coefficient
	vector<lower=0>[num_coef] sd_coef_ ;

	// chol_corr: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[num_coef] chol_corr ;

	// multi_normal_helper: a helper variable for implementing non-centered parameterization
	matrix[num_coef,num_subj] multi_normal_helper ;

}
model{

	////
	// Priors
	////

	// normal(0,1) priors on all means
	mean_coef_ ~ std_normal() ;

	// normal(0,1) priors on all sds
	sd_coef_ ~ std_normal() ;

	// relatively-flat prior on correlations
	chol_corr ~ lkj_corr_cholesky(1) ;

	// multi_normal_helper must have normal(0,1) prior for non-centered parameterization
	to_vector(multi_normal_helper) ~ std_normal() ;


	////
	// Compute latent/hierarchical quantities
	////

	// compute coefficients for each subject/condition
	matrix[num_subj,num_coef] subj_coef_ = (
		rep_matrix(mean_coef_,num_subj)
		+ transpose(
			diag_pre_multiply(sd_coef_,chol_corr)
			* multi_normal_helper
		)
	) ;

	// Loop over subj and conditions to compute unique entries in design matrix
	matrix[num_rows_uW,num_subj] lrt_loc_for_subj_cond ;
	matrix[num_rows_uW,num_subj] lrt_scale_for_subj_cond ;
	matrix[num_rows_uW,num_subj] err_logodds_for_subj_cond ;
	for(this_subj in 1:num_subj){
		for(this_condition in 1:num_rows_uW){
			lrt_loc_for_subj_cond[this_condition,this_subj] = (
				dot_product(
					subj_coef_[this_subj,1:num_cols_uW]
					, uW[this_condition]
				)
			) ;
			lrt_scale_for_subj_cond[this_condition,this_subj] = (
				exp(
					dot_product(
						subj_coef_[this_subj,(num_cols_uW+1):(num_cols_uW*2)]
						, uW[this_condition]
					)
				)
			) ;
			err_logodds_for_subj_cond[this_condition,this_subj] = (
				obs_err_intercept
				+ dot_product(
					subj_coef_[this_subj,(num_cols_uW*2+1):(num_cols_uW*3)]
					, uW[this_condition]
				)
			) ;
		}
		// // slightly less explicit than the above loop but equally fast:
		// value_for_subj_cond[,this_subj] = rows_dot_product(
		// 	rep_matrix(
		// 		subj_coef_[this_subj]
		// 		, num_rows_uW
		// 	)
		// 	, W
		// ) ;
	}


    ////
    // Compute likelihood
    ////

	// lrt likelihood
	lrt_ ~ normal(
		to_vector(lrt_loc_for_subj_cond)[lrt_index]
		, to_vector(lrt_scale_for_subj_cond)[lrt_index]
	) ;

	//// err likelihood (SLOW version):
	// err ~ bernoulli_logit(
	// 	to_vector(err_logodds_for_subj_cond)[err_index]
	// ) ;

	// err likelihood (FAST version):
	err_num_err ~ binomial(
		err_num_obs
		, inv_logit(to_vector(err_logodds_for_subj_cond))[err_cell_index]
	) ;

}
