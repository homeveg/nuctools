#!/usr/bin/env Rscript

# load required packages
if (!require("devtools")) install.packages("https://cran.r-project.org/src/contrib/devtools_1.12.0.tar.gz",
                                         repos=NULL, method="libcurl")
if (!require("pacman")) install.packages("https://cran.r-project.org/src/contrib/pacman_0.4.6.tar.gz",
                                         repos=NULL, method="libcurl")
pacman::p_load( ggplot2, grid, plyr, zoo, optparse )

##############################################################################

#####################
# define functions
#####################
# return linear model equation as a character string 
lm_eqn <- function(df) {
  m <- lm(y ~ x, df)
  eq <- substitute( 
    italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
    list( a = format(coef(m)[1], digits = 2), 
          b = format(coef(m)[2], digits = 2), 
          r2 = format( summary(m)$r.squared, digits = 3 ) 
          ) 
    )
  as.character( as.expression(eq) )  

  }

#############
# load data
#############
load_data<-function(wd,sample) {
  path.to_matrix<-paste(wd,sample, sep="/")
  print(path.to_matrix)
  data <- as.data.frame(read.delim(path.to_matrix, header=FALSE, dec=".", sep="\t"))
  data <- data[,colSums(is.na(data))<nrow(data)]
  if (is.null(dim(data))) {
    df <- data.frame( distance=seq(1,length(data)), frequency = data )
  } else {
    df <- data.frame( distance=seq(1,length(rowSums(data))), frequency = rowSums(data))
  }
  #remove first element
  df<-df[-1,]
  print(head(df))
  return(df)
}

##############
# find peaks
##############
argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.predict <- predict(loess(y ~ x, ...), se = TRUE)
  y.max <- rollapply(zoo(y.predict$fit), 2*w+1, max, align="center")
  delta <- y.max - y.predict$fit[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.predict$fit, y.se=y.predict$s )
  
  #plot(x,y)
  #lines(x,y.predict$fit, col="red")
  #lines(x,y.predict$fit+2*y.predict$s, lty=2) #rough & ready CI
  #lines(x,y.predict$fit-2*y.predict$s, lty=2)
}

################################################################
# ggplot loess appriximation and inclusion with regression line 
################################################################
plot_approx<-function(df, sampleID, subtitle, sample_info, peaks_tweek, wd, output, ylim, xlim, span, win) {
  ymin<-ylim[1]
  ymax<-ylim[2]
  if(max(df$frequency)<ymax) { ymax<-max(df$frequency) }
  if(max(df$distance)<xlim) { xlim<-max(df$distance) }
  
  x<-df$distance
  y<-df$frequency
  
  # estimate peaks position from loess fit
  peaks_delta<-c()
  peaks_delta[1]<-0
  peaks <- argmax(x, y, w=win, span=span)
  peak_minus1<-peaks$i[1]
  
  for (k in 2:length(peaks$i)) {
    peaks_delta<-peaks$i[k]-peak_minus1
    if ( peaks_delta <= 2*win ) {
      print(paste("removing peak Nr.",k))
      peaks$i[k-1]<-peaks$i[k-1] + floor(peaks_delta/2)
      peak_minus1<-peaks$i[k-1]
      peaks$i[k]<-NA
    } else {
      peak_minus1<-peaks$i[k]
    }
  }
  peaks$i<-peaks$i[!is.na(peaks$i)]
  
  peaks.df <-data.frame(distance=x[peaks$i], frequency=peaks$y.hat[peaks$i])
  
  df<-cbind(df, loess_line=peaks$y.hat)
  df<-cbind(df, plus_se=peaks$y.hat+2*peaks$y.se)
  df<-cbind(df, minus_se=peaks$y.hat-2*peaks$y.se)
  
  p <- ggplot(df, aes(x=distance, y=frequency)) + 
    #stat_smooth(size = 1, se = TRUE, method = "loess", span = span, aes(outfit=fit<<-..y..)) + 
    geom_point(size=0.05) + 
    geom_line(aes(x=distance, y=loess_line), col="royalblue", size=1) +
    geom_line(aes(x=distance, y=plus_se), lty=2, col="royalblue", alpha=.5) +
    geom_line(aes(x=distance, y=minus_se), lty=2, col="royalblue") +
    geom_ribbon(aes(x=distance, ymax=plus_se, ymin=minus_se), fill="royalblue", alpha=.2)

  p <-p + scale_x_continuous(limits = c(0, xlim))
  p <-p + scale_y_continuous(limits = c(ymin, ymax))
  
  p <-p + ggtitle(paste(sampleID,": ",subtitle,"\n",sample_info, sep="")) +
    xlab("Distance from origine (bp)") + 
    ylab("Counts") + 
    theme(
      plot.title = element_text(colour="black", size = 20, face = "bold", hjust = "0.5"),
      axis.text.x = element_text(colour="black", size = 16),
      axis.text.y = element_text(colour="black", size = 16),
      axis.title.x = element_text(colour="black", size = 17, face = "bold"),
      axis.title.y = element_text(colour="black", size = 17, face = "bold")
    ) + 
    # draw peaks tops
    geom_point(data=peaks.df, col="red", aes(x=distance, y=frequency)) + 
    geom_text(data=peaks.df,aes(label=format(round_any(peaks.df$frequency,1000,f=ceiling), scientific=TRUE)),hjust=0, vjust=0, nudge_x = 5, nudge_y = 5, size=5, angle = 45)
  
  # linear regression on peaks 
  reg.df <- data.frame(x=seq(1:dim(peaks.df)[1]), y=peaks.df$distance)
  reg_fit<-lm(y ~ x, reg.df)
  summary(reg_fit)
  reg.equation<-lm_eqn(reg.df)
  g <- ggplot(reg.df, aes(x=x, y=y)) + geom_smooth(se = FALSE, method = "lm", formula=y~x-1) + 
    geom_point() + geom_text(aes(label=y),hjust=1, vjust=2, size=5)
  g <- g + geom_text(x = 2, y = (9*max(reg.df$y)/10), hjust=0, label = reg.equation, parse = TRUE, size=6) 
  g1 = ggplotGrob(g)
  # draw plot with inlet
  
  p1 <- p + annotation_custom(grob = g1, xmin = xlim/2, xmax = xlim, ymin = ymin+(ymax-ymin)/2, ymax = ymax)
  #save output
  #p1
  return(p1)
}

##############################################################################

############################
# read command line options
############################
peaks_tweek<-1

# read command line arguments
option_list = list(
  make_option("--input", type="character", default=NULL, 
              help="Nucleosome repeat length file name", metavar="NRL.txt"),
  make_option("--out", type="character", default="out.PNG", 
              help="output figure file name [default= %default]", metavar="NRLplot.png"),
  make_option("--dir", type="character", default=NULL, 
              help="path to working directory", metavar="path"),
  make_option("--info", type="character", default=NULL, 
              help="Sample description", metavar="sample info"),
  make_option("--sample", type="character", default="Sample ID", 
              help="Sample name", metavar="sample ID"),
  make_option("--maxX", type="integer", default=1000, 
              help="X-axis maximum [default= %default]", metavar="1000"),
  make_option("--maxY", type="integer", default=10000, 
              help="Y-axis maximum [default= %default]", metavar="10000"),
  make_option("--minY", type="integer", default=1000, 
              help="Y-axis minimum [default= %default]", metavar="1000"),
  make_option("--span", type="double", default=0.05, 
              help="LOESS aproximation span [default= %default]", metavar="[0-1]"),
  make_option("--peakW", type="integer", default=20, 
              help="peak width and minimal disatnce between 2 adjacent peaks [default= %default]", metavar="20")
); 

# parse arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check if required arguments specified
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Please specify input file name.", call.=FALSE)
}

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Please specify path to a working directory.", call.=FALSE)
}
wd <- opt$dir

sampleID <- opt$sample
sample_info <- opt$info
xlim<-opt$maxX
ylim<-c(opt$minY,opt$maxY)
span<-opt$span
w<-opt$peakW

subtitle.NRL<-c("Nucleosome repeat length")
df.NRL<-load_data(wd,opt$input)
dim(df.NRL)

figure.NRL<-plot_approx(df.NRL[1:xlim,], sampleID, subtitle.NRL, sample_info, peaks_tweek, wd,output.NRL, ylim, xlim, span, w)
print(figure.NRL)
ggsave(paste(wd,opt$out,sep="/"), device = "png", width = 10, height = 10, dpi = 300)
print(paste("saving figure to", paste(wd,opt$out,sep="/") ) )
#dev.off()


