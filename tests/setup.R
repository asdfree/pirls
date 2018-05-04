if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )

options("lodown.cachaca.savecache"=FALSE)

library(lodown)
this_sample_break <- Sys.getenv( "this_sample_break" )
pirls_cat <- get_catalog( "pirls" , output_dir = file.path( getwd() ) )
record_categories <- ceiling( seq( nrow( pirls_cat ) ) / ceiling( nrow( pirls_cat ) / 4 ) )
pirls_cat <- pirls_cat[ record_categories == this_sample_break , ]
pirls_cat <- lodown( "pirls" , pirls_cat )
if( any( pirls_cat$year == 2016 ) ){











library(survey)
library(mitools)

# load the ASG (student background) + ASH (home background) merged design
pirls_design <- readRDS( file.path( getwd() , "2016/asg_design.rds" ) )

# optional step to limit memory usage
variables_to_keep <-
	c( 'idcntry' , 'itsex' , 'asdage' , 'asrrea' , 'asrlit' )
	
pirls_design$designs <-
	lapply( 
		pirls_design$designs ,
		function( w ) {
			w$variables <- w$variables[ variables_to_keep ]
			w
		}
	)

gc()

pirls_design <- 
	update( 
		pirls_design , 
		
		one = 1 ,
		
		idcntry = factor( idcntry ) ,
		
		sex = factor( itsex , labels = c( "male" , "female" ) ) ,
		
		age_ten_or_older = as.numeric( asdage >= 10 )

	)
MIcombine( with( pirls_design , svyby( ~ one , ~ one , unwtd.count ) ) )

MIcombine( with( pirls_design , svyby( ~ one , ~ idcntry , unwtd.count ) ) )
MIcombine( with( pirls_design , svytotal( ~ one ) ) )

MIcombine( with( pirls_design ,
	svyby( ~ one , ~ idcntry , svytotal )
) )
MIcombine( with( pirls_design , svymean( ~ asrrea ) ) )

MIcombine( with( pirls_design ,
	svyby( ~ asrrea , ~ idcntry , svymean )
) )
MIcombine( with( pirls_design , svymean( ~ sex , na.rm = TRUE ) ) )

MIcombine( with( pirls_design ,
	svyby( ~ sex , ~ idcntry , svymean , na.rm = TRUE )
) )
MIcombine( with( pirls_design , svytotal( ~ asrrea ) ) )

MIcombine( with( pirls_design ,
	svyby( ~ asrrea , ~ idcntry , svytotal )
) )
MIcombine( with( pirls_design , svytotal( ~ sex , na.rm = TRUE ) ) )

MIcombine( with( pirls_design ,
	svyby( ~ sex , ~ idcntry , svytotal , na.rm = TRUE )
) )
MIcombine( with( pirls_design ,
	svyquantile(
		~ asrrea ,
		0.5 , se = TRUE 
) ) )

MIcombine( with( pirls_design ,
	svyby(
		~ asrrea , ~ idcntry , svyquantile ,
		0.5 , se = TRUE ,
		keep.var = TRUE , ci = TRUE 
) ) )
MIcombine( with( pirls_design ,
	svyratio( numerator = ~ asrlit , denominator = ~ asrrea )
) )
sub_pirls_design <- subset( pirls_design , idcntry %in% c( 36 , 40 , 31 , 957 ) )
MIcombine( with( sub_pirls_design , svymean( ~ asrrea ) ) )
this_result <-
	MIcombine( with( pirls_design ,
		svymean( ~ asrrea )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	MIcombine( with( pirls_design ,
		svyby( ~ asrrea , ~ idcntry , svymean )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( pirls_design$designs[[1]] )
MIcombine( with( pirls_design , svyvar( ~ asrrea ) ) )
# SRS without replacement
MIcombine( with( pirls_design ,
	svymean( ~ asrrea , deff = TRUE )
) )

# SRS with replacement
MIcombine( with( pirls_design ,
	svymean( ~ asrrea , deff = "replace" )
) )
MIsvyciprop( ~ age_ten_or_older , pirls_design ,
	method = "likelihood" , na.rm = TRUE )
MIsvyttest( asrrea ~ age_ten_or_older , pirls_design )
MIsvychisq( ~ age_ten_or_older + sex , pirls_design )
glm_result <- 
	MIcombine( with( pirls_design ,
		svyglm( asrrea ~ age_ten_or_older + sex )
	) )
	
summary( glm_result )
australia_usa_design <- subset( pirls_design , idcntry %in% c( 36 , 840 ) )

rm( pirls_design ) ; gc()

results <-
	MIcombine( 
		with( 
			australia_usa_design , 
			svyby( 
				~ asrrea , 
				~ idcntry , 
				svymean 
			) 
		) 
	)

stopifnot( round( coef( results )[1] , 2 ) == 544.36 )
stopifnot( round( SE( results )[1] , 2 ) == 2.53 )
stopifnot( round( coef( results )[2] , 2 ) == 549.44 )
stopifnot( round( SE( results )[2] , 2 ) == 3.09 )

}
