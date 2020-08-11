#For building pkdown website:

#Builds reference index
pkgdown::build_reference_index()

#Builds site; do before pushing stable version
pkgdown::build_site()