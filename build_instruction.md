## in order to get Rcpp loaded properly 
## so that headers are discoverable
usethis::use_rcpp()

## to compile cpp files, run
Rcpp::compileAttributes()

## to make man/ documentation, 
## and NAMESPACE
## write in roxygen format,
## add @export and run
devtools::document()