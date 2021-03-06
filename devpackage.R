library(testthat)
library(devtools)
library(usethis)

use_gpl3_license("Pierre")

use_git()
#use_github("git@github.com:Demiperimetre/PIPScreening.git")
# en ligne de commande depuis les commandes de github

use_readme_rmd() # creer la description

use_lifecycle_badge("experimental")
badgecreatr::badge_license()
badgecreatr::badge_codecov(ghaccount = "demiperimetre",ghrepo="PIPScreening",branch="master")

use_r("toySimulator")

use_r("utilsFunctions")
use_rcpp_armadillo("utils")
use_tidy_description()
#use_test("utils.cpp")
use_test("BridgeIS")
usethis::use_rcpp_armadillo()
#use_r("PIPScreening-package")

check()
load_all() # pour charger les fonctinos
document() #generer les doc
test() # bug si on ne refait pas un buil ?

covr::package_coverage() # test
covr::report() # pour avoir le rapport


#### plantage tres bizarre sur les tests... Recomplier le src ? virer l'environnement ???

?sim3


usethis::use_vignette("ExampleScenario")

# tests
goodpractice::goodpractice()
#use_testthat() # dosseier test
#use_test("toySimulator.R")

Rcpp::compileAttributes()
