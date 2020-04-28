library(testthat)
library(devtools)
library(usethis)

use_gpl3_license("Pierre")

use_git()
#use_github("git@github.com:Demiperimetre/PIPScreening.git")
# en ligne de commande depuis les commandes de github

use_readme_rmd() # creer la description

use_lifecycle_badge("experimental")


use_r("toySimulator")

use_tidy_description()

check()
load_all() # pour charger les fonctinos
document() #generer les doc
test()
?sim3


# tests
goodpractice::goodpractice()
use_testthat() # dosseier test
use_test("toySimulator.R")


