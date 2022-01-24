# change the TippingPoint.default from the TippingPoint R package to calculate the binary p-values 
# with fisher.exact(., tsmethod="central") from the exact2x2 R package, instead of using prop.test
# with the continuity correction
library(exact2x2)

## Get other functions needed from the TippingPoint R package
#  from the /R/utilities.R file in the package..need to get that from the source files
## format cat output
format_cat<-function(x) format(as.character(sprintf("%.7g", x)),width =8,justify="left")


## check range of values
check.range <- function(q,q.set) {q[q<=max(q.set) & q>=min(q.set)]}


TippingPoint.modified<-function (outcome, treat, group.infor = FALSE, plot.type = c("estimate", 
                                                             "p.value", "both"), summary.type = c("density", 
                                                                                                  "credible.region", "convex.hull"), alpha = 0.95, 
          HistMeanT = NULL, HistMeanC = NULL, ind.values = FALSE, impValuesT = NA, 
          impValuesC = NA, impValuesColor = NA, show.points = TRUE, 
          point.size = 1, point.shape = 19, S = 3, n.grid = 150, ...) 
{
  # Color Modification: Change all colour= to "black" for the options we used
  # not all the color modifications are marked....
  if (is.numeric(outcome)) {
    Yobs <- outcome
  }
  else stop("Outcome should be a numeric vector!", call. = FALSE)
  if (all(!(treat %in% c(0, 1)))) 
    stop("treat should be a numeric vector with value of 0 or 1 !", 
         call. = FALSE)
  stopifnot(is.logical(group.infor), is.logical(ind.values), 
            is.logical(show.points))
  if (alpha > 1 || alpha < 0) 
    stop("alpha must be a single number between 0 and 1", 
         call. = FALSE)
  plot.type <- match.arg(plot.type)
  summary.type <- match.arg(summary.type)
  if (length(unique(Yobs[!is.na(Yobs)])) > 2){ TPtype <- "continuous"
   } else TPtype <- "binary"
  if (TPtype == "continuous") {
    X_t_mis = seq.int(mean(Yobs[treat == 1], na.rm = TRUE) - 
                        S * sd(Yobs[treat == 1], na.rm = TRUE), mean(Yobs[treat == 
                                                                            1], na.rm = TRUE) + S * sd(Yobs[treat == 1], na.rm = TRUE), 
                      length.out = n.grid)
    X_c_mis = seq.int(mean(Yobs[treat == 0], na.rm = TRUE) - 
                        S * sd(Yobs[treat == 0], na.rm = TRUE), mean(Yobs[treat == 
                                                                            0], na.rm = TRUE) + S * sd(Yobs[treat == 0], na.rm = TRUE), 
                      length.out = n.grid)
    Welch.t.test = function(meanYt_mis, meanYc_mis, treat, 
                            Yobs) {
      Nt = sum(treat == 1)
      Nc = sum(treat == 0)
      N = length(treat)
      L = length(meanYt_mis)
      Ntmis = sum(is.na(Yobs) & treat == 1)
      Ntobs = Nt - Ntmis
      meanYobsT = mean(Yobs[treat == 1], na.rm = TRUE)
      Ncmis = sum(is.na(Yobs) & treat == 0)
      Ncobs = Nc - Ncmis
      meanYobsC = mean(Yobs[treat == 0], na.rm = TRUE)
      p.val = array(NA, L)
      for (i in 1:L) {
        d = (meanYobsT * Ntobs + meanYt_mis[i] * Ntmis)/Nt - 
          (meanYobsC * Ncobs + meanYc_mis[i] * Ncmis)/Nc
        Vart = (var(Yobs[treat == 1], na.rm = TRUE) * 
                  (Ntobs - 1) + Ntobs * Ntmis/Nt * (meanYobsT - 
                                                      meanYt_mis[i])^2)/(Ntobs)
        Varc = (var(Yobs[treat == 0], na.rm = TRUE) * 
                  (Ncobs - 1) + Ncobs * Ncmis/Nc * (meanYobsC - 
                                                      meanYc_mis[i])^2)/(Ncobs)
        Welch.DF = (Vart/Nt + Varc/Nc)^2/((Vart/Nt)^2/(Ntobs) + 
                                            (Varc/Nc)^2/(Ncobs))
        t.stat = d/sqrt(Vart/Nt + Varc/Nc)
        p.val[i] = 2 * pt(-abs(t.stat), df = Welch.DF)
      }
      p.val
    }
    effect.size.mean = function(meanYt_mis, meanYc_mis, treat, 
                                Yobs) {
      (meanYt_mis * sum(is.na(Yobs[treat == 1])) + sum(Yobs[treat == 
                                                              1], na.rm = TRUE))/sum(treat == 1) - (meanYc_mis * 
                                                                                                      sum(is.na(Yobs[treat == 0])) + sum(Yobs[treat == 
                                                                                                                                                0], na.rm = TRUE))/sum(treat == 0)
    }
    p.values <- outer(X_t_mis, X_c_mis, FUN = Welch.t.test, 
                      treat = treat, Yobs = Yobs)
    theta <- outer(X_t_mis, X_c_mis, FUN = effect.size.mean, 
                   treat, Yobs)
    colnames(theta) = X_c_mis
    rownames(theta) = X_t_mis
    df <- melt(theta)
    names(df) <- c("MisYt", "MisYc", "value")
    df = data.frame(df)
    df$p.value = melt(p.values)[, 3]
    p <- ggplot(df, aes(x = MisYt, y = MisYc))
    if (!is.null(HistMeanT)) 
      p <- p + geom_vline(xintercept = HistMeanT, colour = "purple", 
                          lwd = 1.01, alpha = 0.7, lty = 1)
    if (!is.null(HistMeanC)) 
      p <- p + geom_hline(yintercept = HistMeanC, colour = "purple", 
                          lwd = 1.01, alpha = 0.7, lty = 1)
    if (plot.type == "estimate") {
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect")
      if (min(df$value) > 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "white", high = "orange", 
                                     space = "Lab")
      if (max(df$value) < 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "steelblue", high = "white", 
                                     space = "Lab")
      if (max(df$value) > 0 & min(df$value) < 0) 
        p <- p + scale_fill_gradient2("estimate", 
                                      low = "steelblue", high = "orange", 
                                      space = "Lab")
    }
    if (plot.type == "p.value") {
      p <- p + geom_tile(aes(fill = p.value)) + labs(title = "P-value from a hypothesis test")
      # Modification 2: Change color to black-gray-white
      #p <- p + scale_fill_gradient2("p.value", low = "white", 
      #                              high = "darkolivegreen", space = "Lab", 
      #                              limits = c(0, 1))
       p <- p + scale_fill_gradient2("p.value", low = "white", 
                                    high = "black", space = "Lab", 
                                    limits = c(0, 1))


      if (sum(df$p.value <= 0.05) > 0) {
        p <- p + geom_tile(data = subset(df, p.value <= 
                                           0.05), fill = "darkred", alpha = 0.1) + 
          stat_contour(aes(z = p.value), breaks = c(0.05), 
                       color = "darkred")
       # Modification 3: Change color to black-gray-white
       #  p <- p + geom_tile(data = subset(df, p.value <= 
       #                                    0.05), fill = "darkred", alpha = 0.1) + 
       #   stat_contour(aes(z = p.value), breaks = c(0.05), 
       #                color = "darkred")
        p <- p + geom_tile(data = subset(df, p.value <= 
                                           0.05), fill = "black", alpha = 0.1) + 
          stat_contour(aes(z = p.value), breaks = c(0.05), 
                       color = "black")

      }
    }
    if (plot.type == "both") {
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect and significant p-value in red grid if any")
      if (min(df$value) > 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "white", high = "orange", 
                                     space = "Lab")
      if (max(df$value) < 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "steelblue", high = "white", 
                                     space = "Lab")
      if (max(df$value) > 0 & min(df$value) < 0) 
        p <- p + scale_fill_gradient2("estimate", 
                                      low = "steelblue", high = "orange", 
                                      space = "Lab")
      if (sum(df$p.value <= 0.05) > 0) {
        p <- p + geom_tile(data = subset(df, p.value <= 
                                           0.05), fill = "darkred", alpha = 0.1) + 
          stat_contour(aes(z = p.value), breaks = c(0.05), 
                       color = "darkred")
      }
    }
    p <- p + labs(x = "Average outcome for nonrespondents in treatment group", 
                  y = "Average outcome for nonrespondents in control group")
    # Modification 4: change color to black
    #p <- p + geom_hline(yintercept = mean(Yobs[treat == 0], 
    #                                      na.rm = TRUE), colour = "darkblue", lty = 2, 
    #                    alpha = 0.7) + geom_vline(xintercept = mean(Yobs[treat == 
    #                                                                       1], na.rm = TRUE), colour = "darkblue", lty = 2, 
    #                                              alpha = 0.7)
    p <- p + geom_hline(yintercept = mean(Yobs[treat == 0], 
                                          na.rm = TRUE), colour = "black", lty = 2, 
                        alpha = 0.7) + geom_vline(xintercept = mean(Yobs[treat == 
                                                                           1], na.rm = TRUE), colour = "black", lty = 2, 
                                                  alpha = 0.7)

    ContMinMax = check.range(range(Yobs[treat == 0], na.rm = TRUE), 
                             X_c_mis)
    if (length(ContMinMax)) 
      # Modification 5: change color to black
      #p <- p + geom_hline(yintercept = range(Yobs[treat == 
      #                                              0], na.rm = TRUE), colour = "dark blue", 
      #                    alpha = 0.5)
      p <- p + geom_hline(yintercept = range(Yobs[treat == 
                                                    0], na.rm = TRUE), colour = "black", 
                          alpha = 0.5)
    TreatMinMax = check.range(range(Yobs[treat == 1], na.rm = TRUE), 
                              X_t_mis)
    if (length(TreatMinMax)) 
      # modification 6: change color to black
      #p <- p + geom_vline(xintercept = range(Yobs[treat == 
      #                                              1], na.rm = TRUE), colour = "dark blue", 
      #                    alpha = 0.5)
      p <- p + geom_vline(xintercept = range(Yobs[treat == 
                                                    1], na.rm = TRUE), colour = "black", 
                          alpha = 0.5)
    if (group.infor) {
      cat("\nGroup Information:\n\n")
      cat("Groups                    \tTreatment \tControl\n")
      cat("Size                     \t", format_cat(sum(treat)), 
          "\t", format_cat(sum(treat == 0)), "\n")
      cat("Number of nonrespondents \t", format_cat(sum(is.na(Yobs[treat == 
                                                                     1]))), "\t", format_cat(sum(is.na(Yobs[treat == 
                                                                                                              0]))), "\n")
      cat("% of nonrespondents      \t", format_cat(sum(is.na(Yobs[treat == 
                                                                     1]))/sum(treat)), "\t", format_cat(sum(is.na(Yobs[treat == 
                                                                                                                         0]))/sum(treat == 0)), "\n")
      cat("Observed average response\t", format_cat(mean(Yobs[treat == 
                                                                1], na.rm = TRUE)), "\t", format_cat(mean(Yobs[treat == 
                                                                                                                 0], na.rm = TRUE)), "\n")
      cat("Observed min response    \t", format_cat(min(Yobs[treat == 
                                                               1], na.rm = TRUE)), "\t", format_cat(min(Yobs[treat == 
                                                                                                               0], na.rm = TRUE)), "\n")
      cat("Observed max response    \t", format_cat(max(Yobs[treat == 
                                                               1], na.rm = TRUE)), "\t", format_cat(max(Yobs[treat == 
                                                                                                               1], na.rm = TRUE)), "\n")
    }
  }
  if (TPtype == "binary") {
    X_t_mis = seq(0, sum(is.na(Yobs) & treat == 1), 1)
    X_c_mis = seq(0, sum(is.na(Yobs) & treat == 0), 1)
    proportion.test = function(sumYt_mis, sumYc_mis, treat, 
                               Yobs) {
      L = length(sumYt_mis)
      p.val = array(NA, L)
      #######################################################
      # HERE IS THE BIG MODIFICATION
      #######################################################
      ### Change the following to use Fisher's exact p-values
      #for (i in 1:L) p.val[i] = prop.test(c(sumYt_mis[i] + 
      #                                        sum(Yobs * treat, na.rm = TRUE), sumYc_mis[i] + 
      #                                        sum(Yobs * (1 - treat), na.rm = TRUE)), c(sum(treat == 
      #                                         1), sum(treat == 0)))$p.value
      FET<-function(x,n){
        fisher.exact(matrix(c(x[1],n[1]-x[1],x[2],n[2]-x[2]),2,2),tsmethod="central")$p.value
      }
      for (i in 1:L) p.val[i] = FET(c(sumYt_mis[i] + 
                                              sum(Yobs * treat, na.rm = TRUE), sumYc_mis[i] + 
                                              sum(Yobs * (1 - treat), na.rm = TRUE)), c(sum(treat == 
                                                                                              1), sum(treat == 0)))
      #######################################################
      # END OF THE MODIFICATION
      #######################################################
      
      p.val
    }
    effect.size.sum = function(sumYt_mis, sumYc_mis, treat, 
                               Yobs) {
      (sumYt_mis + sum(Yobs[treat == 1], na.rm = TRUE))/sum(treat == 
                                                              1) - (sumYc_mis + sum(Yobs[treat == 0], na.rm = TRUE))/sum(treat == 
                                                                                                                           0)
    }
    p.values <- outer(X_t_mis, X_c_mis, FUN = proportion.test, 
                      treat, Yobs)
    theta <- outer(X_t_mis, X_c_mis, FUN = effect.size.sum, 
                   treat, Yobs)
    colnames(theta) = X_c_mis
    rownames(theta) = X_t_mis
    df <- melt(theta)
    names(df) <- c("MisYt", "MisYc", "value")
    df = data.frame(df)
    df$p.value = melt(p.values)[, 3]
    p <- ggplot(df, aes(x = MisYt, y = MisYc))
    if (!is.null(HistMeanT)) {
      HistValueT = HistMeanT * sum(treat == 1) - sum(Yobs[treat == 
                                                            1], na.rm = TRUE)
      HistValueT = check.range(HistValueT, X_t_mis)
      if (!is.null(HistValueT)) 
        # Color Mod: to black
        #p <- p + geom_vline(xintercept = HistValueT, 
        #                    colour = "purple", lwd = 1.01, alpha = 0.7, 
        #                    lty = 1)
        p <- p + geom_vline(xintercept = HistValueT, 
                            colour = "black", lwd = 1.01, alpha = 0.7, 
                            lty = 1)
    }
    if (!is.null(HistMeanC)) {
      HistValueC = HistMeanC * sum(treat == 0) - sum(Yobs[treat == 
                                                            0], na.rm = TRUE)
      HistValueC = check.range(HistValueC, X_c_mis)
      if (!is.null(HistValueC)) 
        # color mod: to black
        #p <- p + geom_hline(yintercept = HistValueC, 
        #                    colour = "purple", lwd = 1.01, alpha = 0.7, 
        #                    lty = 1)
        p <- p + geom_hline(yintercept = HistValueC, 
                            colour = "black", lwd = 1.01, alpha = 0.7, 
                            lty = 1)
    }
    if (plot.type == "estimate") {
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect")
      if (min(df$value) > 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "white", high = "orange", 
                                     space = "Lab")
      if (max(df$value) < 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "steelblue", high = "white", 
                                     space = "Lab")
      if (max(df$value) > 0 & min(df$value) < 0) 
        p <- p + scale_fill_gradient2("estimate", 
                                      low = "steelblue", high = "orange", 
                                      space = "Lab")
    }
    if (plot.type == "p.value") {
      p <- p + geom_tile(aes(fill = p.value)) + labs(title = "P-value from a hypothesis test")
      p <- p + scale_fill_gradient2("p.value", low = "white", 
                                    high = "black", space = "Lab", 
                                    limits = c(0, 1))
      if (sum(df$p.value > 0.05) > 0) 
        p <- p + geom_tile(data = subset(df, p.value > 
                                           0.05), alpha = 0, colour = "white")
      if (sum(df$p.value <= 0.05) > 0) 
        p <- p + geom_tile(data = subset(df, p.value <= 
                                           0.05), color = "black", alpha = 0)
    }
    if (plot.type == "both") {
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect and significant p-value in red grid if any")
      if (min(df$value) > 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "white", high = "orange", 
                                     space = "Lab")
      if (max(df$value) < 0) 
        p <- p + scale_fill_gradient("estimate", 
                                     low = "steelblue", high = "white", 
                                     space = "Lab")
      if (max(df$value) > 0 & min(df$value) < 0) 
        p <- p + scale_fill_gradient2("estimate", 
                                      low = "steelblue", high = "orange", 
                                      space = "Lab")
      if (sum(df$p.value > 0.05) > 0) 
        p <- p + geom_tile(data = subset(df, p.value > 
                                           0.05), alpha = 0, colour = "white")
      if (sum(df$p.value <= 0.05) > 0) 
        p <- p + geom_tile(data = subset(df, p.value <= 
                                           0.05), colour = "darkred", alpha = 0)
    }
    if (ind.values == TRUE && (plot.type != "p.value")) 
      p <- p + geom_text(aes(x = MisYt, y = MisYc, label = round(value, 
                                                                 2)), color = "black", alpha = 0.6, size = 3)
    if (ind.values == TRUE && (plot.type == "p.value")) 
      p <- p + geom_text(aes(x = MisYt, y = MisYc, label = round(p.value, 
                                                                 2)), color = "black", alpha = 0.6, size = 3)
    ################################
    # Modification.. change label
    ###############################
    p <- p + labs(x = "Number of successes among nonrespondents in group 1", 
                  y = "Number of successes among nonrespondents in group 2")
    p <- p + geom_hline(yintercept = mean(Yobs[treat == 0], 
                                          na.rm = TRUE) * sum(is.na(Yobs[treat == 0])), colour = "black", 
                        lty = 2, alpha = 0.7) + geom_vline(xintercept = mean(Yobs[treat == 
                                                                                    1] * sum(is.na(Yobs[treat == 1])), na.rm = TRUE), 
                                                           colour = "black", lty = 2, alpha = 0.7)
    if (group.infor) {
      cat("\nGroup Information:\n\n")
      cat("Groups                    \tTreatment \tControl\n")
      cat("Size                     \t", format_cat(sum(treat)), 
          "\t", format_cat(sum(treat == 0)), "\n")
      cat("Number of nonrespondents \t", format_cat(sum(is.na(Yobs[treat == 
                                                                     1]))), "\t", format_cat(sum(is.na(Yobs[treat == 
                                                                                                              0]))), "\n")
      cat("% of nonrespondents      \t", format_cat(sum(is.na(Yobs[treat == 
                                                                     1]))/sum(treat)), "\t", format_cat(sum(is.na(Yobs[treat == 
                                                                                                                         0]))/sum(treat == 0)), "\n")
      cat("Observed proportion      \t", format_cat(mean(Yobs[treat == 
                                                                1], na.rm = TRUE)), "\t", format_cat(mean(Yobs[treat == 
                                                                                                                 0], na.rm = TRUE)), "\n")
    }
  }
  if (!all(is.na(impValuesT))) {
    if (all(is.na(impValuesC))) 
      stop("No simulations for control group are specified.")
    n.methods = dim(impValuesT)[2]
    MeanTtrx <- vector("list", n.methods)
    if (!all(is.na(impValuesColor))) 
      Col.pal <- impValuesColor
    else {
      Col.pal = brewer.pal(9, "Set1")
      if (n.methods < 3) {
        Col.pal <- Col.pal[c(5, 6)]
      }
      else Col.pal <- Col.pal[1:n.methods]
    }
    for (i in 1:n.methods) {
      if (summary.type == "density") {
        MeanTtrx[[i]] = cbind(impValuesT[, i], impValuesC[, 
                                                          i])
        p = p + geom_density2d(data = data.frame(MeanTtrx[[i]], 
                                                 row.names = NULL), aes(X1, X2), col = Col.pal[i], 
                               alpha = 1, bins = 6, size = 0.7)
      }
      if (summary.type == "convex.hull") {
        RemoveOutliers <- function(Mtrx) {
          Sx <- cov(Mtrx)
          MahDist <- mahalanobis(Mtrx, colMeans(Mtrx), 
                                 Sx)
          Mtrx = Mtrx[MahDist < quantile(MahDist, alpha), 
                      ]
          Mtrx
        }
        MeanTtrx[[i]] = RemoveOutliers(cbind(impValuesT[, 
                                                        i], impValuesC[, i]))
        hpts <- chull(MeanTtrx[[i]])
        hpts <- c(hpts, hpts[1])
        p = p + geom_polygon(data = data.frame(MeanTtrx[[i]][hpts, 
                                                             ], row.names = NULL), aes(X1, X2), col = Col.pal[i], 
                             alpha = 0.3, fill = Col.pal[i])
      }
      if (summary.type == "credible.region") {
        MeanTtrx[[i]] = cbind(impValuesT[, i], impValuesC[, 
                                                          i])
        Rect = bayesSurv::credible.region(MeanTtrx[[i]], 
                                          probs = alpha)[[1]]
        p = p + geom_rect(aes_string(xmin = Rect[1, 1], 
                                     ymin = Rect[1, 2], xmax = Rect[2, 1], ymax = Rect[2, 
                                                                                       2]), fill = Col.pal[i], alpha = 0, color = Col.pal[i], 
                          size = 1)
      }
      if (show.points) {
        if (TPtype == "continuous") 
          p = p + geom_point(size = point.size, shape = point.shape, 
                             data = data.frame(X1 = impValuesT[, i], X2 = impValuesC[, 
                                                                                     i]), aes(X1, X2), col = Col.pal[i])
        if (TPtype == "binary") 
          p = p + geom_jitter(size = point.size, shape = point.shape, 
                              data = data.frame(X1 = impValuesT[, i], X2 = impValuesC[, 
                                                                                      i]), aes(X1, X2), col = Col.pal[i])
      }
    }
  }
  p
}