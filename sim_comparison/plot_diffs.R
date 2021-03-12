require(ggplot2)

options(warn=1)

# Either specify a data frame for `dd` or a single column.
# If you do the former, you must also supply pdb and lr.
num_outside <- function(dd, pdb=NULL, lr=NULL, m_iqr=5) {
  if (!is.null(pdb)) {
    ecol <- dd[dd$j==0 & dd$pdb==pdb & dd$lr == lr & dd$types %in% c('p p', 'h h'),'i.energy']
  } else {
    ecol = dd
  }
  Q <- quantile(ecol, na.rm=T)
  iqr <- IQR(ecol, na.rm=T)
  n_over <- sum(ecol > Q[4]+m_iqr*iqr)
  n_under <- sum(ecol < Q[2]-m_iqr*iqr)
  #data.frame(n_under=n_under, n_over=n_over, under=Q[2]+m_iqr*iqr, over=Q[4]+m_iqr*iqr);
  c(Q[2]-m_iqr*iqr, Q[4]+m_iqr*iqr, n_over, n_under)
}
remove_outliers <- function(dat, col) {
  dat <- dat[!is.na(dat[,col]),]
  outs <- num_outside(dat[,col], m_iqr=5)
  dat[dat[,col] > outs[1] & dat[,col] < outs[2],]
}

normalize <- function(x) {
  cat('before', summary(x), '\n')
  #x <- scale(x)
  #x <- exp(x)
  xmx <- max(x)
  xmn <- min(x)

  x <- (x - xmn) / (xmx - xmn)
  #print(c(xmx, xmn))
  cat('after', summary(x), '\n')
  #x <- log(x)
  x <- exp(x)
  x
}
# Get density of points in 2 dimensions.
# Source: https://slowkow.com/notes/ggplot2-color-by-density/
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n=100, ...) {
  dens <- MASS::kde2d(x, y, n=n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(normalize(dens$z[ii]))
  #return(dens$z[ii])
}

read_amber <- function(pdb, lr) {
  fn = paste(pdb, '_', lr, '_u_stats.txt', sep='')
  print(paste('AMBER: reading from', fn))
  dat.amber <- read.table(fn, na.strings=c('na', 'NA'))
  names(dat.amber) <- c('pdb', 'stdev', 'lr', 'type', 'conf', 'a', 'b', 'c', 'd', 'e', 'crmsd', 'ucrmsd', 'amber.energy', 'cvc.energy')
  dat.amber$type = 'Amber'
  dat.amber$stdev = ''
  dat.amber <- remove_outliers(dat.amber, 'cvc.energy')
  print(summary(dat.amber))
  dat.amber$density = get_density(dat.amber$cvc.energy, dat.amber$crmsd)
  dat.amber$udensity = get_density(dat.amber$cvc.energy, dat.amber$ucrmsd)
  dat.amber
}
read_prody <- function(pdb, lr, r) {
  fn = paste('r',r,'_', pdb, '_', lr, '_u.prot_en_crmsd.txt', sep='')
  print(paste('PRODY: reading from', fn))
  dat.prody <- read.table(fn)
  names(dat.prody) <- c('conf', 'cvc.energy', 'amber.energy', 'crmsd', 'ucrmsd', 'pd.rmsd')
  dat.prody$lr = lr
  dat.prody$type = 'ProDy'
  dat.prody$stdev = paste('rmsd=', r, sep='')
  dat.prody$pdb = pdb
  dat.prody <- remove_outliers(dat.prody, 'cvc.energy')
  print(summary(dat.prody))
  dat.prody$density = get_density(dat.prody$cvc.energy, dat.prody$crmsd)
  dat.prody$udensity = get_density(dat.prody$cvc.energy, dat.prody$ucrmsd)
  dat.prody
}
read_mine <- function(pdb, lr, stdev=5) {
  fn = paste('hp.5.s',stdev,'-',pdb,'_', lr, '_u-conv.prot_en_ll.txt', sep='')
  internal_read_mine(pdb, lr, stdev, fn)
}
read_lc <- function(pdb, lr, stdev=5) {
  if (pdb != '1H1V' && pdb != '1F6M') {
    data.frame(type=character(),pdb=character(),lr=character(),conf=numeric(),crmsd=numeric(),ucrmsd=numeric(),cvc.energy=numeric(),amber.energy=numeric(), density=numeric(), udensity=numeric(), stdev=numeric())
  } else {
    fn = paste('lc_hp.5.s',stdev,'-',pdb,'_', lr, '_u-conv.prot_en_ll.txt', sep='')
    dat <- internal_read_mine(pdb, lr, stdev, fn)
    dat$type = 'vMRSHD-lc'
    dat
  }
}
internal_read_mine <- function(pdb, lr, stdev=5, fn) {
  print(paste('vMRSHD: reading from', fn))
  dat.mine <- read.table(fn)
  names(dat.mine) <- c('conf', 'a','b','c', 'crmsd', 'ucrmsd', 'amber.energy', 'cvc.energy')
  dat.mine$type = 'vMRSHD'
  dat.mine$pdb = pdb
  dat.mine$lr = lr
  dat.mine$stdev = paste('stdev=', stdev, sep='')
  print(summary(dat.mine))

  dat.mine <- dat.mine[dat.mine$conf <= 1000 | dat.mine$conf == 9999,]
  dat.mine <- remove_outliers(dat.mine, 'cvc.energy')
  dat.mine$density = get_density(dat.mine$cvc.energy, dat.mine$crmsd)
  dat.mine$udensity = get_density(dat.mine$cvc.energy, dat.mine$ucrmsd)
  dat.mine
}

load_dat <- function() {
  dat <- NULL
  for (pdb in c('1H1V', '1F6M')) {
    for (lr in c('r', 'l')) {
      print(paste(pdb, lr))

      print(paste(' vmrshd'))
      dat.mine <- rbind(read_mine(pdb, lr, 1),
                        read_mine(pdb, lr, 5),
                        read_mine(pdb, lr, 50))
      print(paste(' vmrshd-lc'))
      dat.mine_lc <- read_lc(pdb, lr, 5)

      print(paste(' amber'))
      dat.amber <- read_amber(pdb, lr)

      print(paste(' prody'))
      # 5 isn't super great :(
      dat.prody <- rbind(read_prody(pdb, lr, 1),
                         read_prody(pdb, lr, 5),
                         read_prody(pdb, lr, 50))

      this.dat <- rbind(dat.amber[,c('type','pdb','lr','conf','crmsd','ucrmsd','cvc.energy','amber.energy', 'density', 'udensity', 'stdev')],
                        dat.prody[,c('type','pdb','lr','conf','crmsd','ucrmsd','cvc.energy','amber.energy', 'density', 'udensity', 'stdev')],
                        dat.mine [,c('type','pdb','lr','conf','crmsd','ucrmsd','cvc.energy','amber.energy', 'density', 'udensity', 'stdev')],
                        dat.mine_lc [,c('type','pdb','lr','conf','crmsd','ucrmsd','cvc.energy','amber.energy', 'density', 'udensity', 'stdev')])
      if (is.null(dat)) {
        dat <- this.dat
      } else {
        dat <- rbind(dat, this.dat)
      }
    }
  }
  dat$norm <- 'sample'
  dat[dat$conf == "9999","norm"] <- "original"
  dat$lr = ifelse(dat$lr=='l', 'ligand', 'receptor')
  dat$pdb = as.character(dat$pdb)

  # Add some data for the bound energy values.
  en.singles <- read.table("singles.energy.txt")
  names(en.singles) <- c('pdb.full', 'cvc.energy')
  en.singles$pdb <- substr(en.singles$pdb.full, 1, 4)
  en.singles$lr <- ifelse(substr(en.singles$pdb.full, 6, 6) == "l", "ligand", "receptor")
  en.singles$norm <- ifelse(substr(en.singles$pdb.full, 8, 8) == "u", "unbound", "bound")
  en.singles$crmsd <- ifelse(en.singles$norm == "bound", 0, NA)
  en.singles$ucrmsd <- ifelse(en.singles$norm == "bound", NA, 0)
  en.singles$conf <- -1
  dextra <- do.call(rbind, apply(unique(dat[,c("type","stdev","pdb","lr")]), 1,
                       function(x) {
                         d=dat[dat$type==x[1] & dat$stdev==trimws(x[2]) & dat$pdb==x[3] & dat$lr==x[4] & dat$norm=='original',];
                         d$norm = 'bound'
                         d$density = -1
                         d$udensity = -1
                         d$ucrmsd = d$crmsd
                         d$crmsd = 0
                         d$conf = -1
                         d$cvc.energy = en.singles[en.singles$pdb==x[3] & en.singles$lr==x[4] & en.singles$norm=='bound','cvc.energy']
                         d
                      }
                       ))

  rbind(dat, dextra)
}

en_difs <- function() {
  do.call(rbind, lapply(unique(en.singles$pdb), function(pdb) {
    do.call(rbind, lapply(c('ligand', 'receptor'), function(lr) {
        data.frame(pdb=as.character(pdb),
                   lr=lr,
                   diff.1=en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'bound','cvc.energy'],
                   diff.2=en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'unbound','cvc.energy'],
                   diff=en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'bound','cvc.energy'] - en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'unbound','cvc.energy'])
    }))
  }))

  ggplot(do.call(rbind, lapply(unique(en.singles$pdb), function(pdb) {
    do.call(rbind, lapply(c('ligand', 'receptor'), function(lr) {
        data.frame(pdb=as.character(pdb),
                   lr=lr,                                                                                    diff.1=en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'bound','cvc.energy'],
                   diff.2=en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'unbound','cvc.energy'],
                   diff=en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'bound','cvc.energy'] - en.singles[en.singles$pdb == pdb & en.singles$lr == lr & en.singles$norm == 'unbound','cvc.energy'])
    }))
  }))
  , aes(pdb, diff/1e3)) + geom_col() + facet_grid(.~lr)
}

plots <- function() {
  max_centroid_disc <- read.table('max_centroid_disc.txt', h=T)
  levels(max_centroid_disc$stdev) <- c(levels(max_centroid_disc$stdev), '')
  max_centroid_disc[max_centroid_disc$stdev == '-','stdev'] = ''
  mcd <- merge(max_centroid_disc, dat)
  mcd <- mcd[mcd$idx == mcd$conf,]

  ggplot(dat, aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_wrap(pdb+lr~type, scales='free', ncol=3) + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=4, col='red') + scale_y_log10(lim=c(4, 10))

  ggplot(dat[dat$crmsd < 10,], aes(crmsd, fill=interaction(type, as.factor(stdev)))) + geom_density(aes(y=..scaled..), alpha=0.2) + facet_grid(lr+pdb~., scales='free') + geom_vline(dat=dat[dat$norm == "original",], aes(xintercept=crmsd), lty=2, col='red') + geom_vline(xintercept=5, lty=3, col='black')

  ggplot(dat[dat$conf != 9999,], aes(crmsd, ucrmsd)) + stat_bin2d(binwidth=c(1,1)) + facet_grid(pdb+lr~type+stdev, labeller = label_wrap_gen(multi_line=FALSE)) + scale_fill_gradient(low='lightblue', high='darkred', trans='log10') + geom_point(dat=dat[dat$pdb == '1H1V' & dat$norm == "original",], aes(crmsd, ucrmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD')

  ggplot(dat[dat$crmsd < 15 & dat$cvc.energy > -4e5 & dat$cvc.energy < 1e3,], aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_grid(pdb+lr~type+stdev, scales='free') + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=4, col='red')

  ggplot(dat[dat$crmsd < 15 & dat$cvc.energy > -3e4 & dat$cvc.energy < 1e3,], aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_wrap(pdb+lr~type+stdev, scales='free', ncol=4, labeller = label_wrap_gen(multi_line=FALSE)) + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red')

  ggplot(dat[dat$pdb=='1H1V' & dat$crmsd < 15 & dat$cvc.energy > -2.5e4 & dat$cvc.energy < 1e3,], aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_grid(pdb+lr~type+stdev, scales='free') + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$pdb == '1H1V' & dat$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + geom_point(dat=mcd[mcd$pdb == '1H1V' & mcd$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD')
  ggplot(dat[dat$pdb=='1H1V' & dat$ucrmsd < 15 & dat$cvc.energy > -2.5e4 & dat$cvc.energy < 1e3,], aes(cvc.energy/1000, ucrmsd)) + geom_point(aes(color=density), size=0.5) + facet_grid(pdb+lr~type+stdev, scales='free') + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$pdb == '1H1V' & dat$norm == "original",], aes(cvc.energy/1000, ucrmsd), pch=4, size=2, stroke=1, col='red') + geom_point(dat=mcd[mcd$pdb == '1H1V' & mcd$norm == "original",], aes(cvc.energy/1000, ucrmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD')

  ggplot(dat[dat$pdb=='1H1V' & dat$lr=='ligand',], aes(ucrmsd, crmsd)) + geom_point(size=0.5) + facet_grid(pdb+lr~type+stdev, scales='free') + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$pdb == '1H1V' & dat$norm == "original" & dat$lr == 'ligand',], aes(ucrmsd, crmsd), pch=4, size=2, stroke=1, col='red') + geom_point(dat=mcd[mcd$pdb == '1H1V' & mcd$lr == 'ligand',], aes(ucrmsd, crmsd), pch=21, size=2, stroke=1, col='blue', fill=NA) + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD') + xlim(c(NA, 15)) + ylim(c(NA, 15))

  ggplot(dat, aes(cvc.energy/1000, ucrmsd)) + geom_point(aes(color=udensity), size=0.5) + facet_wrap(pdb+lr~type+stdev, scales='free', ncol=7) + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$norm == "original",], aes(cvc.energy/1000, ucrmsd), pch=4, size=2, stroke=1, col='red') + geom_point(dat=dat[dat$norm == "bound",], aes(cvc.energy/1000, ucrmsd), pch=4, size=2, stroke=1, col='blue') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD')
  ggplot(dat, aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_wrap(pdb+lr~type+stdev, scales='free', ncol=7) + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + geom_point(dat=dat[dat$norm == "bound",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='blue') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD')

  q <- ggplot(dat[dat$pdb=='1H1V' & dat$crmsd < 15 & dat$cvc.energy > -2.5e4 & dat$cvc.energy < 1e3,], aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_grid(pdb+lr~type+stdev, scales='free') + scale_color_distiller(palette = "Spectral", direction=-1) + geom_point(dat=dat[dat$pdb == '1H1V' & dat$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD')

  #datpub <- dat %>% filter((pdb=='1H1V' & lr=='ligand') | (pdb=='1F6M' & lr=='receptor'), stdev %in% c('', 'rmsd=5', 'rmsd=50', 'stdev=5', 'stdev=50'))
  datpub <- dat %>% filter(pdb=='1H1V', stdev %in% c('', 'rmsd=5', 'rmsd=50', 'stdev=5', 'stdev=50'))
  my_label_parsed <- function (variable, value) {
    require(plyr)
    if (variable == "lr" || variable == 'pdb') {
      return(as.character(value))
    } else {
      llply(as.character(value), function(x) parse(text = x))
    }
  }

  # Change the label here so we can use 'labeller=label_parsed' and it will give
  # us actual sigmas.
  datpub$label <- factor(paste(datpub$type, datpub$stdev),
                         levels=c(
                    'Amber ',
                    'ProDy rmsd=5',
                    'ProDy rmsd=50',
                    'vMRSHD stdev=5',
                    'vMRSHD stdev=50',
                    'vMRSHD-lc stdev=5'),
                         labels=c(
                    'Amber',
                    'ProDy ~~ rmsd==5',
                    'ProDy ~~ rmsd==50',
                    'vMRSHD ~~ sigma==5',
                    'vMRSHD ~~ sigma==50',
                    'vMRSHD-lc ~~ sigma==5'))
  datpub$label2=paste(datpub$pdb, datpub$lr, sep=': ')
  datpub$pdb2 = paste('"', datpub$pdb, '"', sep='')
  #q <- datpub %>% filter(pdb=='1H1V', !(dat$conf %in% c(-1, 9999))) %>% ggplot(aes(cvc.energy/1000, crmsd)) + geom_bin2d(bins=50) + facet_wrap(lr~label, scales='free', ncol=5) + scale_fill_distiller(palette = "Spectral", direction=-1, trans='log') + geom_point(dat=datpub[datpub$pdb == '1H1V' & datpub$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD') + scale_x_continuous(n.breaks=4)
  q <- datpub %>% filter(pdb=='1H1V', !(conf %in% c(-1, 9999))) %>% ggplot(aes(cvc.energy/1000, crmsd)) + geom_point(aes(color=density), size=0.5) + facet_wrap(lr~label, scales='free', ncol=5, labeller=label_parsed) + scale_color_distiller(palette = "Spectral", direction=-1, trans='log') + geom_point(dat=datpub[datpub$pdb == '1H1V' & datpub$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD') + scale_x_continuous(n.breaks=4)
  ggsave("1H1V_vs_amber_prody.pdf", q, width=10, height=5)

  datpub <- dat %>% filter((pdb=='1H1V' & lr=='ligand') | (pdb=='1F6M' & lr=='receptor'), stdev %in% c('', 'rmsd=5', 'rmsd=50', 'stdev=5', 'stdev=50'))
  datpub$label2=paste(datpub$pdb, datpub$lr, sep=': ')
  q <- datpub %>% filter(!(conf %in% c(-1, 9999))) %>% ggplot(aes(cvc.energy/1000, crmsd)) +
    geom_point(aes(color=density), size=0.5) +
    facet_wrap(lr~label, scales='free', ncol=5, labeller=label_parsed) +
    scale_color_distiller(palette = "Spectral", direction=-1, trans='log') +
    geom_point(dat=datpub[datpub$norm == "original",], aes(cvc.energy/1000, crmsd), pch=4, size=2, stroke=1, col='red') + theme(legend.position = "none") + xlab('Energy (kJ)') + ylab('iRMSD') + scale_x_continuous(n.breaks=4)
  ggsave("both_vs_amber_prody.pdf", q, width=8, height=4)

  do.call(rbind, apply(unique(dat[,c("type", "stdev", 'pdb','lr')]), 1,
                       function(x) {
                         d=dat[dat$type==x[1] & dat$stdev==trimws(x[2]) & dat$pdb==x[3] & dat$lr==x[4],];
                         minx=which.min(d[,'crmsd']);
                         cbind(nrow=nrow(d), d[minx,c('type','stdev','pdb','lr','conf','crmsd','cvc.energy')])
                       }))
  do.call(rbind, apply(unique(dat[,c("type", "stdev", 'pdb','lr')]), 1,
                       function(x) {
                         d=dat[dat$type==x[1] & dat$stdev==trimws(x[2]) & dat$pdb==x[3] & dat$lr==x[4],];
                         minx=which.min(d[d$conf != 9999,'crmsd'])
                         minc=min(d[d$conf != 9999,'crmsd'])
                         maxc=max(d[d$conf != 9999,'crmsd'])
                         avgc=mean(d[d$conf != 9999,'crmsd'])
                         cc = cov(d$cvc.energy, d$crmsd)
                         cdiff = d[d$conf==9999,'crmsd'] - d[minx,'crmsd']
                         cbind(d[minx,c('type','stdev','pdb','lr','crmsd','cvc.energy')], maxc=maxc, avgc=avgc, cdiff=cdiff,ocrmsd=d[d$conf==9999,'crmsd'])
                       }))

  new_disc <- function(dat) {
    startt = Sys.time()
    mydat <- dat
    mydat$crmsd = mydat$crmsd / sd(mydat$crmsd)
    mydat$ucrmsd = mydat$ucrmsd / sd(mydat$ucrmsd)
    mydat$cvc.energy = mydat$cvc.energy / sd(mydat$cvc.energy)
    min(sapply(1:(nrow(mydat)-1), function(i) {
                 if (i %% 100 == 0) {print(paste(i, '/', nrow(mydat), '- time', Sys.time() - startt))}
                 qmn <- min(sapply((i+1):nrow(mydat), function(j) {
                         #qm <- max(sqrt((mydat[i,c('crmsd','ucrmsd','cvc.energy')]-
                         #                mydat[j,c('crmsd','ucrmsd','cvc.energy')])**2), na.rm=T)
                         qm <- max(#sqrt((mydat[i,'crmsd']-mydat[j,'crmsd'])**2),
                                   sqrt((mydat[i,'ucrmsd']-mydat[j,'ucrmsd'])**2),
                                   sqrt((mydat[i,'cvc.energy']-mydat[j,'cvc.energy'])**2), na.rm=T)
                        }) , na.rm=T)
                 #print(paste(i, qmn))
                 qmn
                       }))
  }
  do.call(rbind, apply(unique(dat[,c("type", "stdev", 'pdb','lr')]), 1,
                       function(x) {
                         print(x)
                         d=dat[dat$type==x[1] & dat$stdev==trimws(x[2]) & dat$pdb==x[3] & dat$lr==x[4] & dat$conf <= 1000,];
                         new_disc_d <- new_disc(d)
                         minx=which.min(d[,'crmsd'])
                         cbind(d[minx,c('type','stdev','pdb','lr','crmsd','cvc.energy')], new_disc=new_disc_d, sd.c=sd(d$crmsd), sd.u=sd(d$ucrmsd), sd.e=sd(d$cvc.energy))
                       }))

  t(apply(unique(dat[,c("type", "stdev", 'pdb','lr')]), 1,
                       function(x) {
                         d=dat[dat$type==x[1] & dat$stdev==trimws(x[2]) & dat$pdb==x[3] & dat$lr==x[4] & dat$conf <= 1000,];
                         c(x, crmsd=min(d$crmsd), sd.c=sd(d$crmsd), sd.u=sd(d$ucrmsd), sd.e=sd(d$cvc.energy))
                       }))
}

do_some_summarize <- function() {
  library(dplyr)
  require(tidyr)
  # For table in paper with min/max/avg/Delta
  dat.full %>% group_by(pdb, lr, type, stdev) %>% mutate(orig.crmsd=crmsd[conf==9999]) %>% filter(!(conf %in% c(-1, 9999)), (pdb=='1H1V' & lr=='ligand') | (pdb=='1F6M' & lr=='receptor')) %>% summarize(min=min(crmsd), max=max(crmsd), avg=mean(crmsd), delta=min(orig.crmsd)-min(crmsd)) %>% print(n=30)
  # For table of prob certificates.
  dat.full %>% group_by(pdb, lr, type, stdev) %>% filter(!(conf %in% c(-1, 9999)), (pdb=='1H1V' & lr=='ligand') | (pdb=='1F6M' & lr=='receptor')) %>%
    summarize(cert.10 =sum(crmsd<10)/n(),
              cert.9  =sum(crmsd<9)/n(),
              cert.8  =sum(crmsd<8)/n(),
             # cert.7.5=sum(crmsd<7.5)/n(),
              cert.7  =sum(crmsd<7)/n(),
              cert.6  =sum(crmsd<6)/n(),
              cert.5  =sum(crmsd < 5)/n(),
             # cert.2.5=sum(crmsd<2.5)/n(),
             # cert.1  =sum(crmsd<1)/n(),
        ) %>% print(n=30)

  # An even better plot.
  options(width = 160)
  p <- c(0.90, 0.50, 0.10, 0.05, 0.01, 0.005)
  dat.full %>% group_by(pdb, lr, type, stdev) %>% filter(!(conf == -1), (pdb=='1H1V' & lr=='ligand') | (pdb=='1F6M' & lr=='receptor')) %>% summarize(orig=min(crmsd[conf == 9999]), quants=list(quantile(crmsd[conf != 9999], p)), min=min(crmsd[conf != 9999])) %>% unnest_wider(quants) %>% print(width=100)
} 


# Generate something like this:
#for rl in r l; do for s in 1 5 50; do echo $rl $s; ./prob_diff hp.5.s${s}-1F6M_${rl}_u-conv.prot_en_ll.txt 0 4 1 > min_diff/vmrshd_1F6M_s${s}_${rl}.crmsd; done; done
#for rl in r l; do for r in 1 5 50; do echo $rl $r; ./prob_diff r${r}_1F6M_${rl}_u.prot_en_crmsd.txt 0 3 1 > min_diff/prody_1F6M_r${r}_${rl}.crmsd; done; done
#for rl in r l; do echo $rl; ./prob_diff 1F6M_${rl}_u_stats.txt 4 10 1 > min_diff/amber_1F6M_${rl}.crmsd; done
plot_min_approach <- function() {
  require(dplyr)

  dat <- data.frame(lr=factor(levels=c('r', 'l')), prog=character(), prog_var=character(), size=numeric(), min_val=numeric(), diff=numeric())
  for (lr in c('r', 'l')) {
    #mine
    for (s in c(1, 5, 50)) {
      ff <- read.table(paste("min_diff/vmrshd_1F6M_s", s, "_", lr, ".crmsd", sep=''))
      names(ff) <- c('size', 'min_val', 'diff')
      ff$lr=lr
      ff$prog="vmrshd"
      ff$prog_var=s
      dat <- rbind(dat, ff)
    }
    #prody
    for (r in c(1, 5, 50)) {
      ff <- read.table(paste("min_diff/prody_1F6M_r", r, "_", lr, ".crmsd", sep=''))
      names(ff) <- c('size', 'min_val', 'diff')
      ff$lr=lr
      ff$prog="prody"
      ff$prog_var=r
      dat <- rbind(dat, ff)
    }
    #amber
    ff <- read.table(paste("min_diff/amber_1F6M_", lr, ".crmsd", sep=''))
    names(ff) <- c('size', 'min_val', 'diff')
    ff$lr=lr
    ff$prog="amber"
    ff$prog_var=''
    dat <- rbind(dat, ff)
  }

  #means <- dat %>% group_by(lr, prog, prog_var, size) %>% summarize(diff=mean(diff), min_val=mean(min_val))
  #ggplot(dat %>% group_by(lr, prog, prog_var) %>% summarize(nsteps=max(size), min_val=mean(min_val[which(size==max(size))])), aes(nsteps, min_val, color=prog, pch=prog_var)) + geom_point(size=3) + facet_wrap(.~lr)

  dat
}
plot_chernoff <- function() {
  require(dplyr)

  dat <- data.frame(lr=factor(levels=c('r', 'l')), prog=character(), prog_var=character(), size=numeric(), min_val=numeric(), diff=numeric())
  for (lr in c('r', 'l')) {
    #mine
    for (s in c(1, 5, 50)) {
      ff <- read.table(paste("min_diff/vmrshd_1F6M_s", s, "_", lr, ".chernoff.crmsd", sep=''))
      names(ff) <- c('size', 'chern_next', 'chern_full')
      ff$lr=lr
      ff$prog="vmrshd"
      ff$prog_var=s
      dat <- rbind(dat, ff)
    }
    #prody
    for (r in c(1, 5, 50)) {
      ff <- read.table(paste("min_diff/prody_1F6M_r", r, "_", lr, ".chernoff.crmsd", sep=''))
      names(ff) <- c('size', 'chern_next', 'chern_full')
      ff$lr=lr
      ff$prog="prody"
      ff$prog_var=r
      dat <- rbind(dat, ff)
    }
    #amber
    ff <- read.table(paste("min_diff/amber_1F6M_", lr, ".chernoff.crmsd", sep=''))
    names(ff) <- c('size', 'chern_next', 'chern_full')
    ff$lr=lr
    ff$prog="amber"
    ff$prog_var=''
    dat <- rbind(dat, ff)
  }

  #means <- dat %>% group_by(lr, prog, prog_var, size) %>% summarize(diff=mean(diff), min_val=mean(min_val))
  #ggplot(dat.chern[dat.chern$lr=='r',] %>% pivot_longer(cols=starts_with('chern_'), names_to='type'), aes(size, value, color=type)) + geom_point() + scale_y_log10() + facet_grid(lr~prog+prog_var)

  dat
}
