# lascaux canary
# bard's caged origin story
# cryogenesis
pirls_MIcombine <-
	function (results, variances, call = sys.call(), df.complete = Inf, ...) {
		m <- length(results)
		oldcall <- attr(results, "call")
		if (missing(variances)) {
			variances <- suppressWarnings(lapply(results, vcov))
			results <- lapply(results, coef)
		}
		vbar <- variances[[1]]
		cbar <- results[[1]]
		for (i in 2:m) {
			cbar <- cbar + results[[i]]
			vbar <- vbar + variances[[i]]
		}
		cbar <- cbar/m
		vbar <- vbar/m

		# MODIFICATION
		# evar <- var(do.call("rbind", results))
		evar <- sum( ( unlist( results ) - cbar )^2 / 4 )

		
		r <- (1 + 1/m) * evar/vbar
		df <- (m - 1) * (1 + 1/r)^2
		if (is.matrix(df)) df <- diag(df)
		if (is.finite(df.complete)) {
			dfobs <- ((df.complete + 1)/(df.complete + 3)) * df.complete *
			vbar/(vbar + evar)
			if (is.matrix(dfobs)) dfobs <- diag(dfobs)
			df <- 1/(1/dfobs + 1/df)
		}
		if (is.matrix(r)) r <- diag(r)
		rval <- list(coefficients = cbar, variance = vbar + evar *
		(m + 1)/m, call = c(oldcall, call), nimp = m, df = df,
		missinfo = (r + 2/(df + 3))/(r + 1))
		class(rval) <- "MIresult"
		rval
	}
library(httr)

tf <- tempfile()

this_url <- "https://pirls2021.org/data/downloads/P21_Data_R.zip"

GET( this_url , write_disk( tf ) , progress() )

unzipped_files <- unzip( tf , exdir = tempdir() )
library(haven)

# limit unzipped files to those starting with `asg` followed by three letters followed by `r5`
asg_fns <-
	unzipped_files[ 
		grepl( 
			'^asg[a-z][a-z][a-z]r5' , 
			basename( unzipped_files ) , 
			ignore.case = TRUE 
		) 
	]

# further limit asg files to the first ten countries
countries_thru_bulgaria <-
	c("aad", "adu", "alb", "are", "aus", "aut", "aze", "bfl", "bfr", "bgr")

fns_thru_bulgaria <-
	paste0( paste0( '^asg' , countries_thru_bulgaria , 'r5' ) , collapse = "|" )

asg_aad_bgr_fns <-
	asg_fns[ grepl( fns_thru_bulgaria , basename( asg_fns ) , ignore.case = TRUE ) ]

pirls_df <- NULL

for( rdata_fn in asg_aad_bgr_fns ){

	this_tbl_name <- load( rdata_fn )
	
	this_tbl <- get( this_tbl_name ) ; rm( this_tbl_name )
	
	this_tbl <- zap_labels( this_tbl )
	
	this_df <- data.frame( this_tbl )
	
	names( this_df ) <- tolower( names( this_df ) )
	
	pirls_df <- rbind( pirls_df , this_df )
	
}

# order the data.frame by unique student id
pirls_df <- pirls_df[ with( pirls_df , order( idcntry , idstud ) ) , ]
# pirls_fn <- file.path( path.expand( "~" ) , "PIRLS" , "this_file.rds" )
# saveRDS( pirls_df , file = pirls_fn , compress = FALSE )
# pirls_df <- readRDS( pirls_fn )
# identify all columns ending with `01` thru `05`
ppv <- grep( "(.*)0[1-5]$" , names( pirls_df ) , value = TRUE )

# remove those ending digits
ppv_prefix <- gsub( "0[1-5]$" , "" , ppv )

# identify each of the possibilities with exactly five matches (five implicates)
pv <- names( table( ppv_prefix )[ table( ppv_prefix ) == 5 ] )

# identify each of the `01` thru `05` plausible value columns
pv_columns <-
	grep( 
		paste0( "^" , pv , "0[1-5]$" , collapse = "|" ) , 
		names( pirls_df ) , 
		value = TRUE 
	)
pv_wide_df <- pirls_df[ c( 'idcntry' , 'idstud' , pv_columns ) ]

pirls_df[ pv_columns ] <- NULL
pv_long_df <- 
	reshape( 
		pv_wide_df , 
		varying = lapply( paste0( pv , '0' ) , paste0 , 1:5 ) , 
		direction = 'long' , 
		timevar = 'implicate' , 
		idvar = c( 'idcntry' , 'idstud' ) 
	)

names( pv_long_df ) <- gsub( "01$" , "" , names( pv_long_df ) )
pirls_long_df <- merge( pirls_df , pv_long_df )

pirls_long_df <- pirls_long_df[ with( pirls_long_df , order( idcntry , idstud ) ) , ]

stopifnot( nrow( pirls_long_df ) == nrow( pv_long_df ) )

stopifnot( nrow( pirls_long_df ) / 5 == nrow( pirls_df ) )
pirls_list <- split( pirls_long_df , pirls_long_df[ , 'implicate' ] )
weights_df <- pirls_df[ c( 'jkrep' , 'jkzone' ) ]

for( j in 1:75 ){
	for( i in 0:1 ){
		weights_df[ weights_df[ , 'jkzone' ] != j , paste0( 'rw' , i , j ) ] <- 1
		
		weights_df[ weights_df[ , 'jkzone' ] == j , paste0( 'rw' , i , j ) ] <- 
			2 * ( weights_df[ weights_df[ , 'jkzone' ] == j , 'jkrep' ] == i )
	}
}

weights_df[ c( 'jkrep' , 'jkzone' ) ] <- NULL

library(survey)
library(mitools)

pirls_design <- 
	svrepdesign(
		weights = ~totwgt ,
		repweights = weights_df , 
		data = imputationList( pirls_list ) ,
		type = "other" ,
		scale = 0.5 ,
		rscales = rep( 1 , 150 ) ,
		combined.weights = FALSE ,
		mse = TRUE
	)
pirls_design <- 
	update( 
		pirls_design , 
		
		one = 1 ,
		
		countries_thru_bulgaria = 
		
			factor( 
			
				as.numeric( idcntry ) ,
				
				levels = c(7842L, 7841L, 8L, 784L, 36L, 40L, 31L, 956L, 957L, 100L) ,

				labels =
					c("Abu Dhabi, UAE", "Dubai, UAE", "Albania", "UAE", "Australia", "Austria",
					"Azerbaijan", "Belgium (Flemish)", "Belgium (French)","Bulgaria")
				
			) ,
		
		sex = factor( itsex , levels = 1:2 , labels = c( "female" , "male" ) ) ,
		
		always_speak_language_of_test_at_home =
			ifelse( asbg03 %in% 1:4 , as.numeric( asbg03 == 1 ) , NA )

	)
pirls_MIcombine( with( pirls_design , svyby( ~ one , ~ one , unwtd.count ) ) )

pirls_MIcombine( with( pirls_design , svyby( ~ one , ~ sex , unwtd.count ) ) )
pirls_MIcombine( with( pirls_design , svytotal( ~ one ) ) )

pirls_MIcombine( with( pirls_design ,
	svyby( ~ one , ~ sex , svytotal )
) )
pirls_MIcombine( with( pirls_design , svymean( ~ asrrea , na.rm = TRUE ) ) )

pirls_MIcombine( with( pirls_design ,
	svyby( ~ asrrea , ~ sex , svymean , na.rm = TRUE )
) )
pirls_MIcombine( with( pirls_design , svymean( ~ countries_thru_bulgaria ) ) )

pirls_MIcombine( with( pirls_design ,
	svyby( ~ countries_thru_bulgaria , ~ sex , svymean )
) )
pirls_MIcombine( with( pirls_design , svytotal( ~ asrrea , na.rm = TRUE ) ) )

pirls_MIcombine( with( pirls_design ,
	svyby( ~ asrrea , ~ sex , svytotal , na.rm = TRUE )
) )
pirls_MIcombine( with( pirls_design , svytotal( ~ countries_thru_bulgaria ) ) )

pirls_MIcombine( with( pirls_design ,
	svyby( ~ countries_thru_bulgaria , ~ sex , svytotal )
) )
pirls_MIcombine( with( pirls_design ,
	svyquantile(
		~ asrrea ,
		0.5 , se = TRUE , na.rm = TRUE 
) ) )

pirls_MIcombine( with( pirls_design ,
	svyby(
		~ asrrea , ~ sex , svyquantile ,
		0.5 , se = TRUE ,
		ci = TRUE , na.rm = TRUE
) ) )
pirls_MIcombine( with( pirls_design ,
	svyratio( numerator = ~ asrlit , denominator = ~ asrrea )
) )
sub_pirls_design <- subset( pirls_design , idcntry %in% c( 36 , 40 , 31 , 956 ) )
pirls_MIcombine( with( sub_pirls_design , svymean( ~ asrrea , na.rm = TRUE ) ) )
this_result <-
	pirls_MIcombine( with( pirls_design ,
		svymean( ~ asrrea , na.rm = TRUE )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	pirls_MIcombine( with( pirls_design ,
		svyby( ~ asrrea , ~ sex , svymean , na.rm = TRUE )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( pirls_design$designs[[1]] )
pirls_MIcombine( with( pirls_design , svyvar( ~ asrrea , na.rm = TRUE ) ) )
# SRS without replacement
pirls_MIcombine( with( pirls_design ,
	svymean( ~ asrrea , na.rm = TRUE , deff = TRUE )
) )

# SRS with replacement
pirls_MIcombine( with( pirls_design ,
	svymean( ~ asrrea , na.rm = TRUE , deff = "replace" )
) )
# MIsvyciprop( ~ always_speak_language_of_test_at_home , pirls_design ,
# 	method = "likelihood" , na.rm = TRUE )
# MIsvyttest( asrrea ~ always_speak_language_of_test_at_home , pirls_design )
# MIsvychisq( ~ always_speak_language_of_test_at_home + countries_thru_bulgaria , pirls_design )
glm_result <- 
	pirls_MIcombine( with( pirls_design ,
		svyglm( asrrea ~ always_speak_language_of_test_at_home + countries_thru_bulgaria )
	) )
	
summary( glm_result )
australia_design <- subset( pirls_design , countries_thru_bulgaria %in% "Australia" )

stopifnot( nrow( australia_design ) == 5487 )

result <- pirls_MIcombine( with( australia_design , svymean( ~ asrrea ) ) )

stopifnot( round( coef( result ) , 3 ) == 540.134 )

stopifnot( round( SE( result ) , 3 ) == 1.728 )

australia_fn <- unzipped_files[ grepl( 'ASGAUS' , basename( unzipped_files ) ) ]
australia_tbl_name <- load( australia_fn )
australia_tbl <- get( australia_tbl_name ) ; rm( australia_tbl_name )
australia_tbl <- zap_labels( australia_tbl )
australia_df <- data.frame( australia_tbl )
names( australia_df ) <- tolower( names( australia_df ) )

estimate <-
	mean( c(
		with( australia_df , weighted.mean( asrrea01 , totwgt ) ) ,
		with( australia_df , weighted.mean( asrrea02 , totwgt ) ) ,
		with( australia_df , weighted.mean( asrrea03 , totwgt ) ) ,
		with( australia_df , weighted.mean( asrrea04 , totwgt ) ) ,
		with( australia_df , weighted.mean( asrrea05 , totwgt ) )
	) )

stopifnot( round( estimate , 3 ) == 540.134 )

for( k in 1:5 ){

	this_variance <- 0
	
	for( j in 1:75 ){
		for( i in 0:1 ){
			this_variance <- 
				this_variance + 
				( 
					weighted.mean( 
						australia_df[ , paste0( 'asrrea0' , k ) ] , 
						ifelse( 
							j == australia_df[ , 'jkzone' ] , 
							australia_df[ , 'totwgt' ] * 2 * ( australia_df[ , 'jkrep' ] == i ) , 
							australia_df[ , 'totwgt' ] 
						)
					) -
					weighted.mean( 
						australia_df[ , paste0( 'asrrea0' , k ) ] , 
						australia_df[ , 'totwgt' ]
					)
				)^2
		}
	}
	
	assign( paste0( 'v' , k ) , this_variance * 0.5 )

}

sampling_variance <- mean( c( v1 , v2 , v3 , v4 , v5 ) )
stopifnot( round( sampling_variance , 3 ) == 2.653 )

imputation_variance <-
	( 6 / 5 ) * 
	( 
		( ( with( australia_df , weighted.mean( asrrea01 , totwgt ) ) - estimate )^2 / 4 ) +
		( ( with( australia_df , weighted.mean( asrrea02 , totwgt ) ) - estimate )^2 / 4 ) +
		( ( with( australia_df , weighted.mean( asrrea03 , totwgt ) ) - estimate )^2 / 4 ) +
		( ( with( australia_df , weighted.mean( asrrea04 , totwgt ) ) - estimate )^2 / 4 ) +
		( ( with( australia_df , weighted.mean( asrrea05 , totwgt ) ) - estimate )^2 / 4 ) 
	)

stopifnot( round( imputation_variance , 3 ) == 0.333 )

stopifnot( round( sampling_variance + imputation_variance , 3 ) == 2.987 )

