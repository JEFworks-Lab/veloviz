## in order to get Rcpp loaded properly 
## so that headers are discoverable
usethis::use_rcpp()

## to compile cpp files, run
Rcpp::compileAttributes()

## to make imported packages discoverable
## and to make C++ functions discoverable
## need to create package file
## manually created veloviz.R

## to make man/ documentation, 
## and NAMESPACE
## write in roxygen format,
## add @export and run
devtools::document()
