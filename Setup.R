packageLoad <-
  function(x) {
    for (i in 1:length(x)) {
      if (!x[i] %in% installed.packages()) {
        install.packages(x[i])
      }
      library(x[i], character.only = TRUE)
    }
  }


# create a string of package names
packages <- c('tidyverse',
              'rmarkdown', 
              'ggthemes', 
              'RColorBrewer',
              'viridis',
              'ggdark',
              'plotly',
              'lterdatasampler',
              'rstatix', 
              'usethis', 
              'gitcreds', 
              'httr', 
              'jsonlite', 
              'dataRetrieval',
              'sf',
              'terra',
              'mlr3',
              'terra',
              'mgcv',
              'lubridate',
              'appeears',
              'keyring', 
              'mlr3learners',
              'data.table'
)

# use the packageLoad function we created on those packages
packageLoad(packages)

