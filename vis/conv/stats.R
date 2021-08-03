deg2rad <- function(d) d * pi / 180

# Everything in degrees.
bivariate_vonmises_deg <- function(phi, psi, mu, nu, sigma) {
  exp(sigma * cos(deg2rad(phi-mu)) + sigma * cos(deg2rad(psi - nu)))
}
bivariate_vonmises <- function(phi, psi, mu, nu, sigma) {
  exp(sigma * cos(phi-mu) + sigma * cos(psi - nu))
}

#conv_vm_i <- function(pdf, vm, xx, phi, psi) {
conv_vm_i <- function(xy, phi, psi) {
  sum(bivariate_vonmises(xy[,1], xy[,2], 0, -0.5*pi, 2) *
      bivariate_vonmises(phi-xy[,1], psi-xy[,2], -pi/2, pi/2, 10))
      #bivariate_vonmises(phi-xy[,1], psi-xy[,2], -pi, pi, 2))
}

conv_vm_test <- function(pdf) {
  #require(movMF)
  # Convolution is defined as:
  # f(x) = \int_{\theta=0}^{2\pi} f_1(\theta)f_2(x-\theta)d\theta
  # f(x, y) = \int_{\phi=0}^{2\pi} \int_{\psi=0}^{2\pi} f_1(\phi, \psi)f_2(x-\phi, y-\psi)
  xx=seq(-pi, pi, by=pi/50)
  xy=expand.grid(xx, xx)
  sapply(xx, function(x) {
     print(x);
    sapply(xx, function(y) {
      conv_vm_i(xy, x, y)
        })
      })
}
conv_vm_dat <- function(xx, pdf, x, y) {
  sum(sapply(xx, function(phi) {
    #print(paste('phi', phi));
    sum(sapply(xx, function(psi) pdf[phi, psi] *
               bivariate_vonmises_deg(x-phi, y-psi, 180, 180, 5)))
        }))
}
conv_vm_long <- function(pdf) {
  xx = 1:360
  sapply(xx, function(x) {
    sapply(xx, function(y) {
      print(paste(x, y));
      conv_vm_dat(xx, pdf, x, y)
    })
  })
}
conv_vm <- function(pdf) {
  #require(imagine)
  xx = 1:360
  vm = matrix(sapply(1:(360^2), function(i)
                     bivariate_vonmises_deg(i/360, i %% 360, 180, 180, 5)),
              ncol=360, byrow=T)
  print(dim(vm))
  print(dim(pdf))
  #convolution2D(pdf, vm)
  convolveMat(pdf, vm)
}

# Works great if you have a maximum of 1 value for each entry.
# Otherwise, you'll need to do something like:
#   expandDF(dat %>% group_by(i, j) %>% summarize(total=sum(prob), .groups='keep'), 'total')
expandDF2 <- function(df) {
  mm <- expandDF(df %>% group_by(i, j) %>% summarize(total=sum(prob), .groups='keep'),
                 'total')
  mm/sum(mm)
}
expandDF <- function(df, col='prob', nrow=360, ncol=360) {
  # 360x360 matrix
  mm <- matrix(0, nrow=nrow, ncol=ncol)

  dt <- as_tibble(df)
  mm[as.matrix(dt[,c('i','j')])] = dt[[col]]
  mm
}
load_dat <- function(fl) {
  dat <- read.table(fl, skip=3)
  names(dat) <- c('phi', 'psi', 'prob')
  dat$i <- dat$phi + 180 + 1
  dat$j <- dat$psi + 180 + 1
  dat
}
dostuff <- function() {
  dat <- read.table('ARG_3_3_3_r1_pdf.txt', skip=3)
  names(dat) <- c('phi', 'psi', 'prob')
  dat$i <- dat$phi + 180 + 1
  dat$j <- dat$psi + 180 + 1
  dat
}


convolveMat <- function(X, kernel) {
  # From https://discourse.mc-stan.org/t/discrete-convolution-by-direct-summation-can-this-be-made-more-efficient/969/3
  nr = dim(X)[1]
  nc = dim(X)[2]
  kl = dim(X)[1]
  kl1 = kl-1
  kernelp = t(kernel)
  Y = matrix(0, nrow=nr-kl1, ncol=nc)

  for (i in 1:(nr-kl1)) {
    Y[i,] = kernelp * X[i:(i+kl1), ]
  }
  Y
}

cdf_sampling <- function(mm, n=1000) {
  nrow <- dim(mm)[1]
  ncol <- dim(mm)[2]

  cdf_p <- cumsum(c(mm))
  samps = runif(n=n)
  # Want the index before it's > the value.
  samp_idxs <- sapply(samps, function(x) max(which(cdf_p <= x)))

  data.frame(phi=samp_idxs %% ncol - 180, psi=samp_idxs %/% ncol - 180)
}

cdf_sampling_v2 <- function(mm, n=1000, by_col=TRUE) {
  if (by_col) {
    cdf_x <- cumsum(colSums(mm)) / sum(mm)
    cdf_xy <- sapply(1:nrow(mm), function(i) cumsum(mm[,i])/sum(mm[,i]))
  } else {
    cdf_x <- cumsum(rowSums(mm)) / sum(mm)
    cdf_xy <- sapply(1:nrow(mm), function(i) cumsum(mm[i,])/sum(mm[i,]))
  }

  samps <- cbind(runif(n=n), runif(n=n))
  samp_idxs <- sapply(1:n, function(i) {
                        # This may be the first thing, so adjust accordingly.
                        col=max(which(cdf_x <= samps[i,1]), 1)
                        row=max(which(cdf_xy[,col] <= samps[i,2]), 1)
                        #print(paste(col, row))
                        c(col, row)
                 })
  if (by_col) {
    data.frame(phi=samp_idxs[2,]-180, psi=samp_idxs[1,]-180)
  } else {
    data.frame(phi=samp_idxs[1,]-180, psi=samp_idxs[2,]-180)
  }
}

# How's our CDF sampling approach work?
plot_diffs <- function() {
  dat <- rbind(load_dat('ASP_3_r1_pdf.txt'), load_dat('ASP_1_r1_pdf.txt'), load_dat('ASP_5_r1_pdf.txt'))
  mm <- expandDF2(dat)
  # Plot the original
  dat %>% group_by(phi, psi) %>% summarize(total=sum(prob)) %>% ggplot(aes(phi,psi, fill=log(total))) + geom_tile()
  samps <- cdf_sampling(mm, n=10000)
  # Plot the new samps
  ggplot(samps, aes(phi, psi)) + geom_bin2d()

  # I should note that the graphs look pretty good.
  p_q <- plot_margs(mm, samps)
}

plot_margs <- function(orig_mm, samps) {
  marg.col <- colSums(orig_mm)
  marg.row <- rowSums(orig_mm)

  mm.samp <- expandDF2(samps %>% group_by(phi,psi) %>%
                                 mutate(i=phi + 180 + 1,j=psi + 180 + 1, prob=n()))
  marg_s.col <- colSums(mm.samp)
  marg_s.row <- rowSums(mm.samp)
  print(length(marg.row))
  print(length(marg_s.row))

  df <- rbind(data.frame(type='original', marg='row', idx=seq(1, length(marg.row)), val=marg.row),
              data.frame(type='original', marg='col', idx=seq(1, length(marg.col)), val=marg.col),
              data.frame(type='samps',    marg='row', idx=seq(1, length(marg_s.row)), val=marg_s.row),
              data.frame(type='samps',    marg='col', idx=seq(1, length(marg_s.col)), val=marg_s.row))
  p <- ggplot(df, aes(idx, val, color=marg, lty=type)) + geom_line()
  df <- rbind(data.frame(marg='row', idx=seq(1, length(marg.row)), val=marg.row, sval=marg_s.row),
              data.frame(marg='col', idx=seq(1, length(marg.col)), val=marg.col, sval=marg_s.col))
  q <- ggplot(df, aes(idx, val-sval, lty=marg)) + geom_line()
  list(p, q)
}

test3 <- function() {
  xx <- seq(-50,49, length.out=100)
  yy <- seq(-30,29, length.out=100)
  fun <- function(x, y) floor(15*exp(-((x-10)^2+y^2)/100)+15*exp(-((x+25)^2+y^2)/100));
  df <- expand.grid(xx,yy)
  names(df) <- c('xx','yy')
  df$total <- sapply(1:nrow(df), function(i) fun(df[i,'xx'], df[i,'yy']))
  df$i <- df$xx + 51
  df$j <- df$yy + 31
  mat <- expandDF(df, col='total', nrow=100, ncol=100)
  mat <- mat/sum(mat)

  samps1 <- cdf_sampling(mat, n=1000)
  samps2 <- cdf_sampling_v2(mat, n=1000)
}
