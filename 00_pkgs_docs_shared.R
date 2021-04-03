suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(pacman))

wkPath <- c('./download', './processData')
for(i in wkPath){
  wkPathi = i
  # wkPathi = paste0(sectionName, '/', i)
  #每一个子项目都含plot、result、input
  if (!dir.exists(wkPathi)) dir.create(wkPathi)
}
rm(list=c('i', 'wkPathi', 'wkPath'))