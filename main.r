#preamble (options, installs, imports & custom functions) ----

options(warn=1) #really should be default in R
`%!in%` = Negate(`%in%`) #should be in base R!

# specify the packages used:
required_packages = c(
	'crayon' #for coloring terminal output
	, 'bayesplot' #for convenient posterior plots
	, 'github.com/stan-dev/cmdstanr' #for Stan stuff
	, 'tidyverse' #for all that is good and holy
)

#load the helper functions:
source('r/helper_functions.r')
#helper_functions.r defines:
# - halfhelmert_contrasts()
# - get_contrast_matrix()
# - install_if_missing()
# - add_diagnostic_bools()
# - print.stan_summary_tbl()
# - add_stan_summary_tbl_class()
# - is_cor_diag()

#install any required packages not already present
install_if_missing(required_packages)

# define a shorthand for the pipe operator
`%>%` = magrittr::`%>%`


# Read & Prep the data for Stan ----

a = readRDS('rds/a.rds')

#make sure everything's a factor with sensible level orders
a$session = factor(a$session,levels=c(1,2))
a$task = factor(a$task,levels=c('n','x'))
a$warning = factor(a$warning,levels=c('lo','hi'))
a$cuing = factor(a$cuing,levels=c('invalid','valid'))
a$flankers = factor(a$flankers,levels=c('incongruent','congruent','neutral'))
a$sstroop = factor(a$sstroop,levels=c('incongruent','congruent'))

#set contrasts; treatment for session, half-helmert for all others
contrasts(a$session) = contr.treatment
contrasts(a$task) = halfhelmert_contrasts
contrasts(a$warning) = halfhelmert_contrasts
contrasts(a$cuing) = halfhelmert_contrasts
contrasts(a$flankers) = halfhelmert_contrasts
contrasts(a$sstroop) = halfhelmert_contrasts

# Compute inputs to model ----

# specify the contrast matrix
# This is complicated in order to get one set of contrast columns for the
# first session and a second set of contrasts columns for the second session.
# Might be more straightforward to have joined two sets of independently-created
# contrasts, but I wanted to try doing it in one pipe
W =
	(
		a
		#get the contrast matrix (wrapper on stats::model.matrix)
		%>% get_contrast_matrix(
			formula = ~ 0 + session + warning*cuing*sstroop*flankers
		)
		#convert to tibble
		%>% tibble::as_tibble(.name_repair='unique')
		#make sure there's no grouping at outset
		%>% dplyr::ungroup()
		#add variable to keep track of original row order
		%>% dplyr::mutate(obs_row= 1:dplyr::n())
		#rename the session column to session1 for neatness
		%>% dplyr::rename(session1=session)
		#pivot the session contrast columns to long
		%>% tidyr::pivot_longer(
			cols=c('session1','session2')
			, names_to = 'which_session_contrast'
			, values_to = 'session_contrast_value'
		)
		#pivot all the other columns
		%>% tidyr::pivot_longer(
			cols = c(-obs_row,-which_session_contrast,-session_contrast_value)
			, names_to = 'which_other_contrast'
			, values_to = 'other_contrast_value'
		)
		# multiply, yielding 0s or original values
		%>% dplyr::mutate(
			other_contrast_value = other_contrast_value*session_contrast_value
		)
		# don't need this column anymore
		%>% dplyr::select(-session_contrast_value)
		# pivot back to wide
		%>% tidyr::pivot_wider(
			names_from = c(which_session_contrast,which_other_contrast)
			, values_from = other_contrast_value
		)
		# make sure original order is reinstated
		%>% dplyr::arrange(obs_row)
		# don't need this column anymore
		%>% dplyr::select(-obs_row)
		#phew!
	)

#quick glimpse; lots of rows
nrow(W)
View(head(W))

# get the unique entries in W
uW = dplyr::distinct(W) %>% dplyr::arrange_all(dplyr::desc)
nrow(uW) #far fewer!

#double check that contrasts are orthogonal
table(cor(uW)) #should be 0s & 1s only

#Now, for each unique condition specified by uW, the stan model will
# work out values for that condition for each subject, and we'll need to index
# into the resulting subject-by-condition matrix. So we need to create our own
# subject-by-condition matrix and get the indices of the observed data into a
# the array produced when that matrix is flattened.
uW_per_S =
	(
		uW
		#first repeat the matrix so there's a copy for each subject
		%>% dplyr::slice(
			rep(
				dplyr::row_number()
				, length(unique(a$id))
			)
		)
		#now add the subject labels
		%>% dplyr::mutate(
			id = rep(sort(unique(a$id)),each=nrow(uW))
		)
		#add row identifier
		%>% dplyr::mutate(
			row = 1:dplyr::n()
		)
	)

a$obs_index_in_uW_per_S =
	(
		W
		%>% dplyr::mutate(
			id = a$id
		)
		# join to the full contrast matrix W
		%>%	dplyr::left_join(
			uW_per_S
			, by = c(names(uW),'id')
		)
		# pull the row label
		%>%	dplyr::pull(row)
	)


err_summary =
	(
		a
		%>% dplyr::group_by(obs_index_in_uW_per_S)
		%>% dplyr::summarise(
			err_num_obs = dplyr::n()
			, err_num_err = sum(error)
		)
	)


#package for Stan
data_for_stan = tibble::lst(

	# lrt: lrt on each trial
	lrt = a$lrt[a$error==0]

	# err: err on each trial
	, err = a$error

	# err_num_obs: number of error observations per cell of the summary
	, err_num_obs = err_summary$err_num_obs

	# err_num_err: number of errors per cell of the summary
	, err_num_err = err_summary$err_num_err

	# uW: unique entries in the within predictor matrix
	, uW = uW

	# lrt_index: index of each lrt in flattened subject-by-condition value matrix
	, lrt_index = a$obs_index_in_uW_per_S[a$error==0]

	# err_index: index of each err in flattened subject-by-condition value matrix
	, err_index = a$obs_index_in_uW_per_S

	# err_cell_index: index of each err summary in flattened subject-by-condition value matrix
	, err_cell_index = err_summary$obs_index_in_uW_per_S

	# num_subj: number of subj
	, num_subj = length(unique(a$id))

	# things below are computable from the above

	# num_lrt: number of trials
	, num_lrt = length(lrt)

	# num_lrt: number of trials
	, num_err = length(err)

	# num_err_cells: number of error summary cells
	, num_err_cells = nrow(err_summary)

	# num_rows_uW: num rows in uW
	, num_rows_uW = nrow(uW)

	# num_cols_uW: num cols in uW
	, num_cols_uW = ncol(uW)

)

# Sample and inspect the model ----
mod = cmdstanr::cmdstan_model('stan/err_lrt_loc_scale_FAST.stan')
fit = mod$sample(
	data = data_for_stan
	, chains = parallel::detectCores()/2-1
	, parallel_chains = parallel::detectCores()/2-1
	, seed = 1
	, iter_warmup = 1e3
	, iter_sampling = 1e3
	, init = 0
	, refresh = 20
)
beepr::beep()

#gather summary (inc. diagnostics)
fit_summary =
	(
		fit$summary(
			variables = c(
				'mean_coef_'
				, 'sd_coef_'
				, 'cor'
				, 'err_logodds_for_subj_cond'
				, 'lrt_loc_for_subj_cond'
				, 'lrt_scale_for_subj_cond'
			)
		)
		%>% dplyr::select(variable,mean,q5,q95,rhat,contains('ess'))
		%>% dplyr::filter(
			!is_cor_diag_or_lower_tri(variable,prefix='cor')
		)
		%>% sort_by_variable_size()
		%>% add_stan_summary_tbl_class()
		%>% add_diagnostic_bools(fit)
	)
print(fit_summary)
beepr::beep()


#checking that all rhats<1.01 & all ESS>1e3
(
	fit_summary
	%>% dplyr::select(rhat,ess_bulk,ess_tail)
	%>% summary()
)


#quick plot:
(
	fit$draws(variables=c('mean_coef_'))
	%>% bayesplot::mcmc_intervals()
)


#save for later
fit$save_object('rds/fit.rds')
beepr::beep()

