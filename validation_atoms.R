#!/usr/bin/Rscript

################################################################################
##
## validation_atoms.R
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
## This script compares matched atoms in _atomfulltable.csv against 
##   independent estimate of total atom counts in _commonatoms.csv
##

library(reshape2)
library(dplyr)

if(!exists("full_join"))
  full_join <- function(...) merge(...,all=TRUE)

## =============================================================================
###_* read args

argv <- commandArgs(TRUE)

files <- c("fulltable"=argv[1], # _atomfulltable.csv
           "atoms"=argv[2])     # _commonatoms.csv
prefix <- argv[3]               # prefix for all output files


## files <- c("fulltable"="apinenemech_MCMgroups_atomfulltable.csv",
##            "atoms"="apinenemech_commonatoms.csv")
## prefix <- "apinenemech"

## =============================================================================
###_* read files

fulltable <- unique(read.csv(files["fulltable"]))
fulltable$type <- factor(fulltable$type)
fulltable$shorttype <- factor(substring(fulltable$type,1,1))

atoms <- unique(read.csv(files["atoms"]))
names(atoms)[-1] <- toupper(substring(names(atoms)[-1],1,1))
atoms <- melt(atoms,id.vars="compound",variable.name="shorttype",value.name="n.ref")

## =============================================================================
###_* manipulate tables

Isempty <- function(x) x==""

Fillna <- function(x,value=0) replace(x,is.na(x),value)

atomtypes <- unique(fulltable %>% select(compound,atom,type,shorttype))

numgroups <- fulltable %>% filter(!Isempty(group)) %>%
  group_by(compound,atom) %>%
    summarize(n=length(group))

matched <- fulltable %>% filter(!Isempty(group)) %>%
  select(compound,atom,shorttype) %>% distinct() %>%
    count(compound,shorttype) %>%
      rename(n.matched=n)

unique.df <- full_join(atomtypes,numgroups)
unique.df$n <- Fillna(unique.df$n)

complete.df <- full_join(atoms,matched)
complete.df$shorttype <- factor(complete.df$shorttype)
complete.df$n.matched <- Fillna(complete.df$n.matched)

## remove elements for which there are no entries
complete.df <- complete.df %>% group_by(shorttype) %>%
  do(if(any(.$n.ref > 0)) . else data.frame()) %>%
    ungroup() %>% mutate(shorttype=factor(shorttype))

## =============================================================================
###_* functions and parameters for plotting

## add perturbations from normal distribution
rjitter <- function(x,sd=.1) unclass(x)+rnorm(length(x),0,sd)

## =============================================================================
###_* PLOT: validation (specificity/uniqueness)

pdf(paste0(prefix,"_validation-specificity.pdf"),
    width=7,height=6.5)
par(mar=c(4,4,1,1),mgp=c(1.8,.2,0),tck=0.02,cex=1.2)
with(unique.df,{
  set.seed(1)  
  plot(rjitter(shorttype),rjitter(n),
       pch=19,col=rgb(t(col2rgb("cornflowerblue"))/255,alpha=.3),
       axes=FALSE,ann=FALSE)
  axis(1,seq(1,nlevels(shorttype)),levels(shorttype))
  axis(2,unique(n),las=1)
  axis(3,,FALSE)
  axis(4,,FALSE)
  box()
  box(bty="L")
})
title(xlab="Atom type",ylab="Number of groups",cex.lab=1.2)
dev.off()

## =============================================================================
###_* PLOT: validation (completeness)

pdf(paste0(prefix,"_validation-completeness.pdf"),
    width=3*nlevels(complete.df$shorttype),height=3.5)
par(mfrow=c(1,nlevels(complete.df$shorttype)))
par(mar=c(2,2,.2,.2),mgp=c(1.8,.2,0),oma=c(2,2,1.5,1.5),tck=0.02)
par(cex=1)
for(.stype in levels(complete.df$shorttype)) {
  with(complete.df %>% filter(shorttype==.stype),{
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
  mtext(.stype,3,line=.5)
}
mtext("True count",1,outer=TRUE,cex=1.2)
mtext("Matched count",2,outer=TRUE,cex=1.2)
dev.off()

## =============================================================================
###_* DIAGNOSTIC: zero, multiple matches

Optionalwrite <- function(subtable,outfile) {
  if(file.exists(outfile)) file.remove(outfile)
  if(nrow(subtable)>0) write.csv(subtable,outfile,row.names=FALSE)
  invisible(NULL)
}

Optionalwrite(filter(unique.df,n==0),
              paste0(prefix,"_zeromatches.csv"))

Optionalwrite(filter(unique.df,shorttype!="C" & n>1),
              paste0(prefix,"_multiplematches.csv"))
