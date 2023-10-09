Supplementary file for the paper “Ultra-high dimensional confounder selection algorithms comparison with application
  to radiomics data” by I. Baldé, and D. Ghosh, for publication in Biometrics (Practice section).

There are Rmarkdown documents and three R examples of code (SIS + GOAL, SIS + OAL and CBS), to help the reader reproduce similar simulation results to the ones obtained in the paper. 
The three R examples code consist of: 

- SIS + GOAL :  R code containing documented functions to run the simulation combinaing SIS and GOAL.

- SIS + OAL :  R code containing documented functions to run the simulation combinaing SIS and OAL.

- CBS       :  R code containing documented functions to run the simulation of CBS.

Note that R packages "devtools", "lqa", "Ball", "glmnet", "grpreg"and "MASS" have to be installed in order to use this code. 


```{r}
# install.packages("devtools") # install "devtool" package first
library(devtools)
# install_github("cran/lqa")  # install "lqa" package from GitHub
library(lqa)
# install.packages("Ball")
library(Ball)
# install.packages("glmnet")
library(glmnet)
# install.packages("grpreg")
library(grpreg)
# install.packages("MASS")
library(MASS) # version 3.3.1

```