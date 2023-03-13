library(ggsci)
library("dplyr")
library("tidyverse")
library("survival")
library("survminer")
library("eulerr")
library("ggplot2")
library("reshape")
library("ggfortify")
library("gridExtra")
library("survMisc")
library("rsq")
library("pROC")
library("data.table")
library("RColorBrewer")
library("gt")
library("DT")
library(plotly)
library(dplyr)
library(carData)
library(gapminder)
library(babynames)
col_pce <- rev(brewer.pal(10, "Paired"))
col_p = c(col_pce[1:3],
          col_pce[5],
          col_pce[4],
          col_pce[c(6, 10)],
          "black",
          col_pce[9],
          "turquoise",
          col_pce[7])
tems = theme(
  axis.text.x = element_text(
    color = "grey20",
    size = 20,
    angle = 90,
    hjust = .5,
    vjust = .5,
    face = "plain"
  ),
  axis.text.y = element_text(
    color = "grey20",
    size = 20,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ),
  axis.title.x = element_text(
    color = "grey20",
    size = 20,
    angle = 0,
    hjust = .5,
    vjust = 0,
    face = "plain"
  ),
  axis.title.y = element_text(
    color = "grey20",
    size = 20,
    angle = 90,
    hjust = .5,
    vjust = .5,
    face = "plain"
  )
)


