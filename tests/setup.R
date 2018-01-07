if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )

library(lodown)
lodown( "pirls" , output_dir = file.path( getwd() ) )
library(lodown)
# examine all available PIRLS microdata files
pirls_cat <-
	get_catalog( "pirls" ,
		output_dir = file.path( getwd() ) )

# 2011 only
pirls_cat <- subset( pirls_cat , year == 2011 )
# download the microdata to your local computer
lodown( "pirls" , pirls_cat )

library(survey)
library(mitools)

# load the ASG (student background) + ASH (home background) merged design
pirls_design <- readRDS( file.path( getwd() , "2011/asg_design.rds" ) )

# optional step to limit memory usage
variables_to_keep <-
	c( 'idcntry' , 'itsex' , 'itbirthy' , 'asrrea' , 'asrlit' )
	
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
		
		born_2001_or_later = as.numeric( itbirthy >= 2001 )

	)
lodown:::pirls_MIcombine( with( pirls_design , svyby( ~ one , ~ one , unwtd.count ) ) )

lodown:::pirls_MIcombine( with( pirls_design , svyby( ~ one , ~ idcntry , unwtd.count ) ) )
lodown:::pirls_MIcombine( with( pirls_design , svytotal( ~ one ) ) )

lodown:::pirls_MIcombine( with( pirls_design ,
	svyby( ~ one , ~ idcntry , svytotal )
) )
lodown:::pirls_MIcombine( with( pirls_design , svymean( ~ asrrea ) ) )

lodown:::pirls_MIcombine( with( pirls_design ,
	svyby( ~ asrrea , ~ idcntry , svymean )
) )
lodown:::pirls_MIcombine( with( pirls_design , svymean( ~ sex , na.rm = TRUE ) ) )

lodown:::pirls_MIcombine( with( pirls_design ,
	svyby( ~ sex , ~ idcntry , svymean , na.rm = TRUE )
) )
lodown:::pirls_MIcombine( with( pirls_design , svytotal( ~ asrrea ) ) )

lodown:::pirls_MIcombine( with( pirls_design ,
	svyby( ~ asrrea , ~ idcntry , svytotal )
) )
lodown:::pirls_MIcombine( with( pirls_design , svytotal( ~ sex , na.rm = TRUE ) ) )

lodown:::pirls_MIcombine( with( pirls_design ,
	svyby( ~ sex , ~ idcntry , svytotal , na.rm = TRUE )
) )
lodown:::pirls_MIcombine( with( pirls_design , svyquantile( ~ asrrea , 0.5 , se = TRUE ) ) )

lodown:::pirls_MIcombine( with( pirls_design ,
	svyby( 
		~ asrrea , ~ idcntry , svyquantile , 0.5 ,
		se = TRUE , keep.var = TRUE , ci = TRUE 
) ) )
lodown:::pirls_MIcombine( with( pirls_design ,
	svyratio( numerator = ~ asrlit , denominator = ~ asrrea )
) )
sub_pirls_design <- subset( pirls_design , idcntry %in% c( 36 , 40 , 31 , 957 ) )
lodown:::pirls_MIcombine( with( sub_pirls_design , svymean( ~ asrrea ) ) )
this_result <-
	lodown:::pirls_MIcombine( with( pirls_design ,
		svymean( ~ asrrea )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	lodown:::pirls_MIcombine( with( pirls_design ,
		svyby( ~ asrrea , ~ idcntry , svymean )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( pirls_design$designs[[1]] )
lodown:::pirls_MIcombine( with( pirls_design , svyvar( ~ asrrea ) ) )
# SRS without replacement
lodown:::pirls_MIcombine( with( pirls_design ,
	svymean( ~ asrrea , deff = TRUE )
) )

# SRS with replacement
lodown:::pirls_MIcombine( with( pirls_design ,
	svymean( ~ asrrea , deff = "replace" )
) )
lodown:::MIsvyciprop( ~ born_2001_or_later , pirls_design ,
	method = "likelihood" , na.rm = TRUE )
lodown:::MIsvyttest( asrrea ~ born_2001_or_later , pirls_design )
lodown:::MIsvychisq( ~ born_2001_or_later + sex , pirls_design )
glm_result <- 
	lodown:::pirls_MIcombine( with( pirls_design ,
		svyglm( asrrea ~ born_2001_or_later + sex )
	) )
	
summary( glm_result )

