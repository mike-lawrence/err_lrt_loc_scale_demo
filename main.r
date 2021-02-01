#required packages: tidyverse, cmdstanr
library(tidyverse)

a = readRDS('rds/a.rds')

# Prep the data for Stan ----

#sort by subject
a %>% arrange(id) -> a

#make sure everything's a factor with sensible level orders
a$session = factor(a$session,levels=c(1,2))
a$task = factor(a$task,levels=c('n','x'))
a$warning = factor(a$warning,levels=c('lo','hi'))
a$cuing = factor(a$cuing,levels=c('invalid','valid'))
a$flankers = factor(a$flankers,levels=c('incongruent','congruent','neutral'))
a$sstroop = factor(a$sstroop,levels=c('incongruent','congruent'))

#define a function for nicer contrasts
halfhelmert_contrasts = function(...) contr.helmert(...)*.5

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
		%>% ezStan::get_contrast_matrix(
			formula = ~ 0 + session + warning*cuing*sstroop*flankers
		)
		#convert to tibble
		%>% tibble::as_tibble(.name_repair='unique')
		%>% dplyr::mutate(obs_row= 1:n())
		%>% tidyr::pivot_longer(cols=c('session','session2'))
		%>% tidyr::pivot_longer(
			cols = c(-obs_row,-name,-value)
			, names_to = 'name2'
			, values_to = 'value2'
		)
		%>% dplyr::mutate(
			value2 = value*value2
		)
		%>% dplyr::select(-value)
		%>% tidyr::pivot_wider(
			names_from = c(name,name2)
			, values_from = value2
		)
		%>% dplyr::arrange(obs_row)
		%>% dplyr::select(-obs_row)
		# %>% dplyr::distinct()
		# %>% View()
	)

#quick glimpse; lots of rows
nrow(W)
# View(W)

# get the unique entries in W
uW = dplyr::distinct(W) %>% dplyr::arrange_all(desc)
nrow(uW)

#double check that contrasts are orthogonal
table(cor(uW)) #should be 0s & 1s only

#for each unique condition specified by uW, the stan model will
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
			err_num_obs = n()
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

mod = cmdstanr::cmdstan_model('stan/err_lrt_loc_scale_FAST.stan')

#sample the model
fit = mod$sample(
	data = data_for_stan
	, chains = parallel::detectCores()/2-1
	, parallel_chains = parallel::detectCores()/2-1
	, seed = 1
	, iter_warmup = 1e3
	# , iter_sampling = ceiling(num_samples_to_obtain/(phys_cores_minus_one*2))
	, iter_sampling = 1e3
	, init = 0
	# update every 10% (yeah, the progress info sucks; we're working on it)
	, refresh = 20 #(1e3*2)/10
	# , adapt_delta = .99
)

# 1911.8 for the err_summary version
fit$cmdstan_diagnose()
beepr::beep()

fit$summary(variables=c('cor'))


#checking that all rhats<1.01 & all ESS>1e3
(
	fit$summary()
	%>% dplyr::select(rhat,ess_bulk,ess_tail)
	%>% summary()
	#NA's are the _extended variables
)


#quick plot:
(
	fit$draws(variables=c('noise','mean_coef','sd_coef','weights'))
	%>% bayesplot::mcmc_intervals()
)


saveRDS(post,file=filename)
beepr::beep()

