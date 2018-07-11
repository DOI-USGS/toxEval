library(readr)
library(tidyr)
library(dplyr)

benchstuff <- read_csv("D:/LADData/RCode/toxEval_Archive/other benchmarks/EPAAqLifeBenchmarksPest.csv")

long_bench <- gather(benchstuff, endpoint, value, -CASN, -Year, -Compound) %>%
  mutate(remark = "")

long_bench$remark[grep(">",long_bench$value)] <- ">"
long_bench$remark[grep("<",long_bench$value)] <- "<"

long_bench$value <- gsub(" ","",long_bench$value)
long_bench$value <- gsub("<","",long_bench$value)
long_bench$value <- gsub(">","",long_bench$value)
long_bench$value <- as.numeric(long_bench$value)
long_bench <- rename(long_bench, CAS=CASN)


write.csv(long_bench, file="long_aqu.csv", row.names = FALSE, na = "")
