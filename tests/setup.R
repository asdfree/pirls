if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )
devtools::install_github("ajdamico/lodown")
library(lodown)
lodown( "pirls" , output_dir = file.path( getwd() ) )
