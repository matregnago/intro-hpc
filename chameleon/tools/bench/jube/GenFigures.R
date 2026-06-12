####
##
##  @file GenFigures.R
##  @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
##                       Univ. Bordeaux. All rights reserved.
##
##  @version 1.4.0
##  @author Mathieu Faverge
##  @author Florent Pruvost
##  @date 2025-12-19
##
####
library(ggplot2)
library(plyr)
library(reshape2)

##
## Load the file and clean up
##
df <- read.csv( "chameleon.csv", header=TRUE, sep=",", comment.char="#" )

## Factorize the parameter for a better use
df$algorithm = factor(df$algorithm)
df$hostname  = factor(df$hostname)
df$mpivendor = factor(df$mpivendor)
df$ngpu      = factor(df$ngpu)
df$nmpi      = factor(df$nmpi)
df$p         = factor(df$p)
df$q         = factor(df$q)
df$precision = factor(df$precision)

## Get the hostname list
hostnames <- ddply( df, c("hostname"), summarise, N = length(hostname) )$hostname

## Get the mpivendor list
mpinames <- ddply( df, c("mpivendor"), summarise, N = length(mpivendor) )$mpivendor

##
## Generate one figure for CPU times, and one figure for GFlop/s
## mpi - string to caracterize the mpivendors included in df
##
gen_mpi_figure <- function( df, mpi )
{
    ## summary( df )
    if ( nrow( df ) == 0 ) {
        return(NULL);
    }

    df$nmpi = revalue( df$nmpi, c("1" = "1 MPI Process",
                                  "4" = "4 MPI Processes (2x2)",
                                  "9" = "9 MPI Processes (3x3)") )

    df$algorithm = revalue( df$algorithm, c("gemm"      = "Matrix Multiply (gemm)",
                                            "geqrf_hqr" = "QR factorization (geqrf)",
                                            "potrf"     = "Cholesky factorization (potrf)") )

    df$precision = revalue( df$precision, c("s" = "Single real (s)",
                                            "d" = "Double real (d)" ) )

    df$mpivendor <- revalue( df$mpivendor, c("openmpi" = "OpenMPI",
                                             "nmad" = "NewMadeleine") )

    ##
    ##  CPU time version
    ##
    p <- ggplot( df, aes( x=n, y=cputime, colour=mpivendor, linetype=precision ) )

    p <- p + facet_grid( algorithm ~ nmpi, scales="free_x" )
    p <- p + geom_line(  )
    p <- p + scale_colour_discrete( name=element_blank() ) ##name="Number of MPI processes" )
    p <- p + scale_linetype_manual( values=c("solid", "dashed"),
                                   name=element_blank() )## name="Arithmetic",

    p <- p + theme( legend.position="bottom" )

    p <- p + scale_y_log10( )
    p <- p + ylab( "CPU Time (s)" )
    p <- p + xlab( "Matrix size (M = N = K)" )

    filename <- paste( "chameleon", "plafrim", h, "time", mpi, sep="_")
    filename <- paste( filename, ".png", sep="" )

    ggsave( filename, plot=p, width=10, height=10 )

    ##
    ##  GFLop/s version
    ##
    p <- ggplot( df, aes( x=n, y=gflops, colour=mpivendor, linetype=precision ) )

    p <- p + facet_grid( algorithm ~ nmpi, scales="free_x" )
    p <- p + geom_line(  )
    p <- p + scale_colour_discrete( name=element_blank() ) ##name="Number of MPI processes" )
    p <- p + scale_linetype_manual( values=c("solid", "dashed"),
                                   name=element_blank() )## name="Arithmetic",

    p <- p + theme( legend.position="bottom" )

    p <- p + ylab( "GFlop/s" )
    p <- p + xlab( "Matrix size (M = N = K)" )

    filename <- paste( "chameleon", "plafrim", h, "perf", mpi, sep="_")
    filename <- paste( filename, ".png", sep="" )

    ggsave( filename, plot=p, width=10, height=10 )
}

##
## Generate one figure with both CPU times and GFlop/s in facet for single node exeperiments
##
gen_shm_figure <- function( df )
{
    if ( nrow( df ) == 0 ) {
        return(NULL);
    }

    df$algorithm = revalue( df$algorithm, c("gemm"      = "Matrix Multiply (gemm)",
                                            "geqrf_hqr" = "QR factorization (geqrf)",
                                            "potrf"     = "Cholesky factorization (potrf)") )

    df$precision = revalue( df$precision, c("s" = "Single real (s)",
                                            "d" = "Double real (d)" ) )

    df$mpivendor <- revalue( df$mpivendor, c("openmpi" = "OpenMPI",
                                             "nmad" = "NewMadeleine") )

    df$type <- revalue( df$type, c("cputime" = "CPU Time (s)",
                                   "gflops"  = "GFlop/s") )

    p <- ggplot( df, aes( x=n, y=value, colour=mpivendor, linetype=precision ) )

    p <- p + facet_grid( type ~ algorithm, scales="free" )
    p <- p + geom_line(  )
    p <- p + scale_linetype_manual( name=element_blank(),
                                   values=c("solid", "dashed"))
    p <- p + theme( legend.position="bottom" )

    p <- p + xlab( "Matrix size (M = N = K)" )
    p <- p + ylab( element_blank() )
    p <- p + guides( colour="none" )

    filename <- paste( "chameleon", "plafrim", h, sep="_")
    filename <- paste( filename, ".png", sep="" )

    ggsave( filename, plot=p, width=10, height=7 )
}

##
## Generate the distributed figures
##

#Restict to bora
hostnames <- c( 'bora' )

for ( h in hostnames ) {
    for ( mpi in mpinames ) {
        ldf <- df[ (df$hostname == h) & (df$mpivendor == mpi),]
        gen_mpi_figure( ldf, mpi )
    }

    ldf <- df[ (df$hostname == h),]
    gen_mpi_figure( ldf, "both" )
}

##
## Generate the distributed figures
##

## Restrict the list
hostnames <- c( 'sirocco' )
mpinames  <- c( 'openmpi' )

# Restructure the data to have cpu time and gflops on separated lines
shmdf <- melt( df,
              id.vars = c( "algorithm", "hostname", "k", "m", "mpivendor", "n",
                          "ngpu", "nmpi", "nthread", "p", "precision", "q" ),
              variable.name = "type",
              value.name    = "value" )

# loop over hostname and mpivendor
for ( mpi in mpinames ) {
    for ( h in hostnames ) {
        ldf <- shmdf[ (shmdf$hostname == h) & (shmdf$mpivendor == mpi),]
        gen_shm_figure( ldf )
    }
}


