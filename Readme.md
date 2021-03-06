Applied Bayesian Inference
================

These are my personal course materials for: Applied Bayesian Inference - with Susumu Shikano, Constance. The repository comprises:

-   `assigns` Couple of Assignments
-   `bib` contains the syllabus and some literature <!-- * `data` pet data -->
-   `labs` scripts and code from lab sessions <!-- * `slides` course presentations -->

Info JAGS
---------

-   Jags another Gibbs Sampler
-   Clone of Bugs
-   BUGS (Bayesian inference Using Gibbs Sampling):
    1.  Bayesian Inference
    2.  Graphical modeling
    3.  simulation-based inference
-   rjags returns a dynamic model and we can draw samples from this model.
-   Bugs is inspired by S
-   for loop in bugs are a macro expansions of single line codes (no control flow statement)

Install JAGS
------------

-   [brew](https://brew.sh/)

``` bash
brew install jags
```

Packages
--------

``` r
# install.packages("pacman")
pacman::p_load(dplyr, ggplot2, rjags, rjags, purrr)
```

Other JAGS interfaces
---------------------

-   `rjags` set up model;
-   `r2jags` set up, burn in and sample.
-   `runjags` long burn in phase (4000)

Define JAGS Models
------------------

``` r
model_string <- "
  model {
    ## priors:
    a ~ dunif(0,5)
    b ~ dunif(-10,10)
    sigma ~ dunif(0,3)
    
    ## structure:
    for (i in 1:N) {
        y[i] ~ dnorm(a * x[i] + b, pow(sigma, -2))
    }
  }
"


lm_code <- function(){
    ## priors:
    a ~ dunif(0,5)
    b ~ dunif(-10,10)
    sigma ~ dunif(0,3)
    
    ## structure:
    for (i in 1:N) {
        y[i] ~ dnorm(a * x[i] + b, pow(sigma, -2))
    }
}

parse_model <- function(x){
  x %>% 
    deparse() %>% 
    glue::glue_collapse(., "\n") %>% 
    stringr::str_replace("function.*?\\(\\)", "model")  %>% 
    textConnection
}

parse_model(lm_code)
```

    ## A connection with                            
    ## description "."             
    ## class       "textConnection"
    ## mode        "r"             
    ## text        "text"          
    ## opened      "opened"        
    ## can read    "yes"           
    ## can write   "no"

``` r
model_string %>% 
    textConnection
```

    ## A connection with                            
    ## description "."             
    ## class       "textConnection"
    ## mode        "r"             
    ## text        "text"          
    ## opened      "opened"        
    ## can read    "yes"           
    ## can write   "no"
