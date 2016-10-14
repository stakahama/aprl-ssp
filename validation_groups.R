#!/usr/bin/env Rscript

################################################################################
##
## validation_groups.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
## -----------------------------------------------------------------------------
##
## This file is part of APRL-SSP
##
## APRL-SSP is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## APRL-SSP is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with APRL-SSP.  If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

##
## This script compares matched groups against
##   independent estimate of group counts
##

library(reshape2)
library(dplyr)

if(!exists("full_join"))
  full_join <- function(...) merge(...,all=TRUE)

## =============================================================================
###_* read args

argv <- commandArgs(TRUE)

files <- c("compounds"=argv[1], # compounds_{groupfile}.csv
           "groups"=argv[2])     # matched_{groupfile}.csv
prefix <- argv[3]               # prefix for output file

## files <- c("compounds"="compounds_FTIRextra.csv",
##            "groups"="matched_FTIRextra.csv")
## prefix <- "FTIRextra"

## =============================================================================
###_* read files

compounds <- read.csv(files["compounds"],check.names=FALSE)
compounds$SMILES <- NULL
compounds$"source" <- NULL

groups <- read.csv(files["groups"],check.names=FALSE)

## =============================================================================
###_* manipulate tables

Fillna <- function(x,value=0) replace(x,is.na(x),value)

paired <- full_join(
  melt(compounds,id.vars="compound",variable.name="group",value.name="n.ref"),
  melt(groups,id.vars="compound",variable.name="group",value.name="n.matched")
  )
paired$n.matched <- Fillna(paired$n.matched)
paired$group <- factor(paired$group)

## =============================================================================
###_* functions and parameters for plotting

## add perturbations from normal distribution
rjitter <- function(x,sd=.1) unclass(x)+rnorm(length(x),0,sd)

npanels <- n2mfrow(nlevels(paired$group))
pdf(paste0(prefix,"_validation-counts.pdf"),
    width=3*npanels[1],height=3*npanels[2])
par(mfrow=npanels)
par(mar=c(2,2,.2,.2),mgp=c(1.8,.2,0),oma=c(2,2,1.5,1.5),tck=0.02)
par(cex=1)
for(.group in levels(paired$group)) {
  with(paired %>% filter(group==.group),{
    lim <- c(0,max(n.ref))
    plot.new()
    plot.window(lim,lim)
    abline(0,1,col=8)
    set.seed(1)
    points(rjitter(n.ref,.01),rjitter(n.matched,.01),
           col="cornflowerblue",pch=19)
    at <- round(axTicks(1))
    axis(1,at)
    axis(2,at,las=1)
    axis(3,at,FALSE)
    axis(4,at,FALSE)
    box()
  })
  ## text(par("usr")[1],par("usr")[4],adj=c(-.5,1.5),.stype,cex=1.2)
  mtext(.group,3,line=0)
}
mtext("True count",1,outer=TRUE,cex=1.2)
mtext("Matched count",2,outer=TRUE,cex=1.2)
dev.off()

###_* diagnostic

Optionalwrite <- function(subtable,outfile) {
  if(file.exists(outfile)) file.remove(outfile)
  if(nrow(subtable)>0) write.csv(subtable,outfile,row.names=FALSE)
  invisible(NULL)
}

Optionalwrite(filter(paired,n.ref!=n.matched),
              paste0(prefix,"_discrepancy.csv"))
