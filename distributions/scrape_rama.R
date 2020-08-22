normalize <- function(x) {
  #cat('before', summary(x), '\n')
  #x <- scale(x)
  #x <- exp(x)
  xmx <- max(x)
  xmn <- min(x)

  x <- (x - xmn) / (xmx - xmn)
  #print(c(xmx, xmn))
  #cat('after', summary(x), '\n')
  #x <- log(x)
  #x <- exp(x)
  x
}

# Get density of points in 2 dimensions.
# Source: https://slowkow.com/notes/ggplot2-color-by-density/
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(normalize(dens$z[ii]))
  #return(dens$z[ii])
}

load_dat <- function() {
  dat <- read.table('cullpdbs.pp.out')
  names(dat) <- c('aa_full', 'phi', 'psi')
  dat$aa <- substring(dat$aa_full, 0, 3)
  dat$ss_ <- as.numeric(substring(dat$aa_full, 5, 6))
  print(unique(dat$aa))
  print(unique(dat$ss_))
  dat <- dat[!(is.nan(dat$psi) | is.nan(dat$phi) | is.infinite(dat$psi) | is.infinite(dat$phi)),]
  dat <- dat[!(dat$psi == 0 & dat$phi == 0),]

  # only care about sheet and helix (1,3,9+5)
  print(nrow(dat))
  dat <- dat[dat$ss_ %in% c(1,3,5,9),]
  print(nrow(dat))
  dat$ss <- factor(ifelse(dat$ss_ == 1, 'sheet', ifelse(dat$ss_ %in% c(9, 5), 'helix', 'loop')), levels=c('helix','sheet','loop'))
  dat$color = ifelse(dat$ss == 'sheet', 'red', ifelse(dat$ss == 'helix', 'blue', 'green'))
  dat
}

plot_ss_prev <- function(dat, n=5e5) {
  dat <- dat[sample(nrow(dat), n),]
  for (s in c('sheet', 'helix', 'loop')) {
    dh <- dat[dat$ss == s,]
    print(paste(s, nrow(dh)))
    dat[dat$ss == s,'density'] <- get_density(dh$psi, dh$phi, n=100)
  }
  breaks <- c(-pi, -pi/2, 0, pi/2, pi)
  labels <- c(expression(-pi),
              expression(-pi/2),
              expression(0),
              expression(pi/2),
              expression(pi))
  q <- ggplot(dat, aes(phi, psi, color=density)) +
    geom_point(size=0.5) +
    scale_color_distiller(palette = "Spectral", direction=-1, values=c(0, 0.1, 0.2, 0.3, 0.5, 1)) +
    facet_grid(.~ss) +
    xlab(expression(phi)) + ylab(expression(psi)) +
    scale_x_continuous(breaks=breaks, labels=labels) +
    scale_y_continuous(breaks=breaks, labels=labels) +
    theme(legend.position = "none")
  #ggplot(dat, aes(phi, psi, color=ss)) + geom_density2d() +
  #  scale_x_continuous(breaks=breaks, labels=labels, limits=c(-pi, pi)) +
  #  scale_y_continuous(breaks=breaks, labels=labels, limits=c(-pi, pi))
}

plot_ss <- function(dat, n=5e5) {
  breaks <- c(-pi, -pi/2, 0, pi/2, pi)
  labels <- c(expression(-pi),
              expression(-pi/2),
              expression(0),
              expression(pi/2),
              expression(pi))
  qo <- ggplot(dat[sample(nrow(dat), n),], aes(phi, psi)) +
    geom_bin2d(bins=100, aes(fill=..density..)) +
    facet_wrap(.~ss) +
    xlab(expression(phi)) + ylab(expression(psi)) +
    scale_x_continuous(breaks=breaks, labels=labels, expand=c(0, 0)) +
    scale_y_continuous(breaks=breaks, labels=labels, expand=c(0, 0)) +
    theme(legend.position = "none") +
    scale_fill_distiller(palette = "Spectral", direction=-1, trans='log10',
                         values=c(0, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1))
                         #values=c(0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.6, 0.7, 0.9, 1))
  ggsave('phi_psi.pdf', qo, width=8, height=3)

  # Separate??
  plot_one <- function(ss) {
    print(paste('plotting', ss))
    q <- ggplot(dat[dat$ss==ss,], aes(phi, psi)) +
      geom_bin2d(bins=100, aes(fill=..density..)) +
      scale_x_continuous(breaks=breaks, labels=labels, expand=c(0, 0)) +
      scale_y_continuous(breaks=breaks, labels=labels, expand=c(0, 0)) +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
      labs(title=ss, x=expression(phi), y=expression(psi)) +
      scale_fill_distiller(palette = "Spectral", direction=-1, trans='log10',
                           values=c(0, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1))
                           #values=c(0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.6, 0.7, 0.9, 1))
    ggsave(paste('phi_psi_', ss, '.pdf', sep=''), q, width=4, height=4)
  }
  for (ss in unique(dat$ss)) {
    plot_one(ss)
  }

  qo
}

StatRasa <- ggplot2::ggproto("StatRasa", ggplot2::Stat,
                             compute_group = function(data, scales, fun, fun.args) {
                               # Change default arguments of the function to the 
                               # values in fun.args
                               args <- formals(fun)
                               for (i in seq_along(fun.args)) {
                                 if (names(fun.args[i]) %in% names(fun.args)) {
                                   args[[names(fun.args[i])]] <- fun.args[[i]]
                                 } 
                               }
                               formals(fun) <- args

                               # Apply function to data
                               fun(data)
                             })

# stat function used in ggplot
stat_rasa <- function(mapping = NULL, data = NULL,
                      geom = "point", 
                      position = "identity",
                      fun = NULL,
                      ...,
                      show.legend = NA,
                      inherit.aes = TRUE) {
  # Check arguments 
  if (!is.function(fun)) stop("fun must be a function")

  # Pass dotted arguments to a list
  fun.args <- match.call(expand.dots = FALSE)$`...`

  ggplot2::layer(
                 data = data,
                 mapping = mapping,
                 stat = StatRasa,
                 geom = geom,
                 position = position,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 check.aes = FALSE,
                 check.param = FALSE,
                 params = list(
                               fun = fun, 
                               fun.args = fun.args,
                               na.rm = FALSE,
                               ...
                 )
  )
}

IrregularContour <- function(data, breaks = scales::fullseq,
                             binwidth = NULL,
                             bins = 10) {
  require(contoureR)
  if (is.function(breaks)) {
    # If no parameters set, use pretty bins to calculate binwidth
    if (is.null(binwidth)) {
      binwidth <- diff(range(data$z)) / bins
    }

    breaks <- breaks(range(data$z), binwidth)
  }
  breaks <- c(0.03, 0.8, 0.1, 0.15, 0.2, 0.25, 0.4, 0.5, 0.6, 1)

  cl <- contoureR::getContourLines(x = data$x, y = data$y, z = data$z,
                                   levels = breaks)

  if (length(cl) == 0) {
    warning("Not possible to generate contour data", call. = FALSE)
    return(data.frame())
  }
  cl <- cl[, 3:7]
  colnames(cl) <- c("piece", "group", "x", "y", "level")
  return(cl)
}

stat_contour_irregular <- function(...) {
  stat_rasa(fun = IrregularContour, geom = "path", ...)
}

plot_rasa <- function(dat) {
  breaks <- c(-pi, -pi/2, 0, pi/2, pi)
  labels <- c(expression(-pi),
              expression(-pi/2),
              expression(0),
              expression(pi/2),
              expression(pi))
  color_values=c(0, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1)
  ggplot(dat[sample(nrow(dat), 50000),], aes(phi, psi)) +
    geom_bin2d(bins=100, aes(fill=..density..)) +
    stat_contour_irregular(aes(z=density, color=..level..)) +
    facet_wrap(.~ss) +
    theme(legend.position = "none") +
    xlab(expression(phi)) + ylab(expression(psi)) +
    scale_x_continuous(breaks=breaks, labels=labels, expand=c(0, 0)) +
    scale_y_continuous(breaks=breaks, labels=labels, expand=c(0, 0)) +
    scale_color_distiller(palette = "Spectral", direction=-1, trans='log10',
                          values=color_values) +
    scale_fill_distiller(palette = "Spectral", direction=-1, trans='log10',
                         values=color_values)
                          #values=c(0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.6, 0.7, 0.9, 1))

}

foobar <- function() {
library(ggplot2)
dat <- load_dat()
  for (s in c('sheet', 'helix', 'loop')) {
    dh <- dat[dat$ss == s,]
    print(paste(s, nrow(dh)))
    dat[dat$ss == s,'density'] <- get_density(dat[dat$ss==s,'psi'], dat[dat$ss==s,'phi'], n=100)
  }
print('writing plot...')
ggsave("rasa_plot.pdf", plot_rasa(dat), width=8, height=3)
}
