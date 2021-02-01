halfhelmert_contrasts = function(...){
	contr.helmert(...)*.5
}

get_contrast_matrix = function(
	data
	, formula
	, contrast_kind = NULL
){
	if (inherits(data, "tbl_df")) {
		data = as.data.frame(data)
	}
	vars = attr(terms(formula),'term.labels')
	vars = vars[!grepl(':',vars)]
	if(length(vars)==1){
		data = data.frame(data[,vars])
		names(data) = vars
	}else{
		data = data[,vars]
	}
	vars_to_rename = NULL
	for(i in vars){
		if(is.character(data[,i])){
			data[,i] = factor(data[,i])
		}
		if( is.factor(data[,i])){
			if(length(levels(data[,i]))==2){
				vars_to_rename = c(vars_to_rename,i)
			}
			if(!is.null(contrast_kind) ){
				contrasts(data[,i]) = contrast_kind
			}
		}
	}
	mm = model.matrix(data=data,object=formula)
	dimnames(mm)[[2]][dimnames(mm)[[2]]=='(Intercept)'] = '(I)'
	for(i in vars_to_rename){
		dimnames(mm)[[2]] = gsub(paste0(i,1),i,dimnames(mm)[[2]])
	}
	attr(mm,'formula') = formula
	attr(mm,'data') = data
	return(mm)
}



#' Installs any packages not already installed
#' @examples
#' \dontrun{
#' install_if_missing(c('tidyverse','github.com/stan-dev/cmdstanr'))
#' }
install_if_missing = function(pkgs){
	missing_pkgs = NULL

	for(this_pkg in pkgs){
		path = NULL
		try(
			path <-	find.package(basename(this_pkg),quiet=T,verbose=F)
			, silent = T
		)
		if(is.null(path)){
			missing_pkgs = c(missing_pkgs,this_pkg)
		}
	}
	cran_missing = missing_pkgs[!grepl('github.com/',fixed=T,missing_pkgs)]
	if(length(cran_missing)>0){
		message('The following required but uninstalled CRAN packages will now be installed:\n',paste(cran_missing,collapse='\n'))
		install.packages(cran_missing)
	}
	github_missing = missing_pkgs[grepl('github.com/',fixed=T,missing_pkgs)]
	github_missing = gsub('github.com/','',github_missing)
	if(length(github_missing)>0){
		message('The following required but uninstalled Github packages will now be installed:\n',paste(this_pkg,collapse='\n'))
		remotes::install_github(github_missing)
	}
	invisible()
}

#define a function that adds diagnostics as a metadata attribute
add_diagnostic_bools = function(x,fit){
	sink('/dev/null')
	diagnostics = fit$cmdstan_diagnose()$stdout #annoyingly not quiet-able
	sink(NULL)
	diagnostic_bools = list(
		treedepth_maxed = stringr::str_detect(diagnostics,'transitions hit the maximum')
		, ebfmi_low = stringr::str_detect(diagnostics,' is below the nominal threshold')
		, essp_low = stringr::str_detect(diagnostics,'The following parameters had fewer than')
		, rhat_high = stringr::str_detect(diagnostics,'The following parameters had split R-hat greater than')
	)
	attr(x,'meta') = list(diagnostic_bools=diagnostic_bools)
	return(x)
}

#define a custom print method that shows the diagnostics in the metadata attribute
print.stan_summary_tbl = function(x,...) {
	meta = attr(x,'meta')
	if(any(unlist(meta$diagnostic_bools))){
		cat(crayon::bgRed('WARNING:\n'))
	}
	if(meta$diagnostic_bools$treedepth_maxed){
		cat(crayon::bgRed('Treedepth maxed\n'))
	}
	if(meta$diagnostic_bools$ebfmi_low){
		cat(crayon::bgRed('E-BMFI low\n'))
	}
	if(meta$diagnostic_bools$essp_low){
		cat(crayon::bgRed('ESS% low for one or more parameters\n'))
	}
	if(meta$diagnostic_bools$rhat_high){
		cat(crayon::bgRed('R-hat high for one or more parameters\n'))
	}
	NextMethod(x,...)
	invisible(x)
}

#create a new S3 class for custom-printing stanfit summary tables
add_stan_summary_tbl_class = function(x){
	class(x) <- c("stan_summary_tbl",class(x))
	return(x)
}

#function to detect whether a variable name indicates that it's on the diagonal
# of a correlation parameter matrix
has_underscore_suffix = function(x){
	bare_has_underscore_suffix = stringr::str_ends(x,'_')
	has_index_suffix = stringr::str_ends(x,']')
	indexed_has_underscore_suffix =
		(
			tibble::tibble(x = x[has_index_suffix])
			%>% tidyr::separate(x,sep='\\[',into='x',extra='drop')
			%>% dplyr::mutate( out=stringr::str_ends(x,'_') )
			%>% dplyr::pull(out)
		)
	bare_has_underscore_suffix[has_index_suffix] = indexed_has_underscore_suffix
	return(bare_has_underscore_suffix)
}


#function to detect whether a variable name indicates that it's on the diagonal
# or lower-tri element of a correlation parameter matrix
is_cor_diag_or_lower_tri = function(x,prefix){
	has_prefix = stringr::str_starts(x,prefix)
	x = x[has_prefix]
	del = function(x,to_del){gsub(to_del,'',x,fixed=T)}
	to_toss =
		(
			x
			%>% del(prefix)
			%>% del('[')
			%>% del(']')
			%>% tibble::tibble(x = .)
			%>% tidyr::separate(x,into=c('i','j'))
			%>% dplyr::mutate( to_toss = (i==j) | (i>j) )
			%>% dplyr::pull(to_toss)
		)
	has_prefix[has_prefix] = to_toss
	return(has_prefix)
}

#function to sort a stan summary table by size of variables
sort_by_variable_size = function(x){
	x2 =
		(
			x
			%>% tidyr::separate(
				variable
				, sep = '\\['
				, into = 'var'
				, extra = 'drop'
				, remove = F
			)
		)
	(
		x2
		%>% dplyr::group_by(var)
		%>% dplyr::summarise(count = dplyr::n(),.groups = 'drop')
		%>% dplyr::full_join(x2,by='var')
		%>% dplyr::arrange(count,var,variable)
		%>% dplyr::select(-count,-var)
	)
}

