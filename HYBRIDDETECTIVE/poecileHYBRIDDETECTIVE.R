#devtools::install_github("rystanley/genepopedit")
#devtools::install_github("bwringe/parallelnewhybrid")
#devtools::install_github("bwringe/hybriddetective")

library(hybriddetective)
library(genepopedit)

getTopLoc(GRP = "genepopHYBRIDDETECTIVE.txt", panel.size = 200, where.PLINK = /data5/K_Feldmann_data/plink, where.PGDSpider = /data5/K_Feldmann_data/PGDSpider_2.1.1.5/PGDSpider2-cli.jar)

