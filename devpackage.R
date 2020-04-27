library(testthat)
library(devtools)
library(usethis)

use_gpl3_license("Pierre")

use_git()
#use_github("git@github.com:Demiperimetre/PIPScreening.git")
# en ligne de commande depuis les commandes de github

use_readme_rmd() # creer la description

use_lifecycle_badge("experimental")



check()
load_all() # pour charger les fonctinos
document() #generer les doc
