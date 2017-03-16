library(tidyverse)

data.files <- dir("raw-data", pattern = "*.csv")
test.data.files <- dir("test-raw-data", pattern = "*.csv")

(missing.files <- setdiff(test.data.files, data.files))

xx <- tibble(file = data.files) %>% 
  group_by(file) %>% 
  do(origin.data = read_csv(file.path("raw-data", .$file)),
     test.data = read_csv(file.path("test-raw-data", .$file))) 

xx[unlist(map2(xx$origin.data, xx$test.data, ~ nrow(.x) != nrow(.y) || tryCatch(nrow(anti_join(.x, .y)) != 0, error = function(err) TRUE))), ]

