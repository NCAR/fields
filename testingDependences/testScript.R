
setwd("~/Dropbox/Home/Repositories/fields/testingDependences")
library( devtools)
library( fields)
pkg<- "spam"
install.packages(pkg, lib="lib",dependencies=TRUE)
library(spam)
check( "lib/spam", manual = FALSE, vignettes= FALSE)


