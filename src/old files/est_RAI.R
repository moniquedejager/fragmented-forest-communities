est_RAI <- function(n, outputtype) {
  # function to plot estimated rank abundance distribution:
  
  # n is the abundance vector (number of individuals per species)
  # outputtype: 'fit' provides the goodness-of-fit;
  #             'plot' provides the SAD figure;
  #             'estimates' provides the METE estimates of the abundances
  #             'differences' provides the differences per rank between the 
  #                           actual abundance and the METE estimate
  
  df2 <- data.frame(spp = 1:length(n), 
                    abund = sort(n, decreasing = T),
                    rank = 1:length(n)) 
  df2$relAbund <- df2$abund / sum(df2$abund)
  
  esf1 <- meteESF(spp=df2$spp, abund=df2$abund)
  sad1 <- sad(esf1)
  sad1 # print function returns useful summary
  
  df2$estRank <- 1:length(n)

  df2$estAbund <- sad1$q(((length(n)-1):0)/length(n))
  df2$estRelAbund <- df2$estAbund/sum(df2$estAbund)
  
  df2$dif <- df2$abund[order(df2$rank)] - df2$estAbund[order(df2$estRank)]
  
  if (outputtype == 'plot'){
    p2 <- ggplot(df2, aes(x=rank, y=abund )) +
      geom_line(aes(x=estRank, y=estAbund), color='mediumseagreen', linewidth=1.1) +
      geom_point() + 
      scale_y_continuous(trans='log10') + 
      xlab('Rank') + 
      ylab('Abundance')
    print(p2)
  }
  
  if (outputtype == 'fit'){
    return(sum(df2$dif^2))
  }
  
  if (outputtype == 'estimates'){
    return(df2$estAbund)
  }
  
  if (outputtype == 'differences'){
    return(df2$dif)
  }
}
