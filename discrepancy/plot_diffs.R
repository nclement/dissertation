dat.r.halton <- read.table('hp.5.s5-1F6M_r_u-prng.stats')
names(dat.r.halton) <- c('i','j','dist', 'd.energy', 'd.crmsd')
dat.r.halton$method = 'halton'
dat.r.halton$lr = 'receptor'

dat.r.prng <- read.table('hp.5.s5-1F6M_r_u-prng.stats')
names(dat.r.prng) <- c('i','j','dist', 'd.energy', 'd.crmsd')
dat.r.prng$method = 'prng'
dat.r.prng$lr = 'receptor'

dat.l.halton <- read.table('hp.5.s5-1F6M_l_u-prng.stats')
names(dat.l.halton) <- c('i','j','dist', 'd.energy', 'd.crmsd')
dat.l.halton$method = 'halton'
dat.l.halton$lr = 'ligand'

dat.l.prng <- read.table('hp.5.s5-1F6M_l_u-prng.stats')
names(dat.l.prng) <- c('i','j','dist', 'd.energy', 'd.crmsd')
dat.l.prng$method = 'prng'
dat.l.prng$lr = 'ligand'

dat <- rbind(dat.l.prng, dat.r.prng, dat.l.halton, dat.r.halton)
ggplot(dat[dat$dist > 0,], aes(x=dist, y=abs(d.energy))) + stat_density2d(aes(fill = ..level..), geom = "polygon") + scale_fill_distiller(palette="Spectral", direction=1) + facet_wrap(~method)


read_stats_table <- function(pdb, lr, hp, s=5) {
  dat <- read.table(paste('hp.', hp, '.s', s, '-', pdb, '_', lr, '_u-both.stats', sep=''), h=T)
  dat$lr = lr
  dat$pdb = pdb
  dat$hp = hp
  dat
}
read_new_stats <- function(pdb, lr, hp, type, s) {
  dat <- read.table(paste('hp.', hp, '.s', s, '-', pdb, '_', lr, '_u-', type, '.prot_en_ll.txt', sep=''))
  names(dat) <- c('conf', 'll', 'a','b', 'crmsd', 'ucrmsd', 'amberen', 'cvc.energy')
  dat$lr = lr
  dat$pdb = pdb
  dat$hp = hp
  dat$s = s
  dat$type = type
  dat
}
## dat.full.1BKD.l <- read_stats_table('1BKD', 'l', 5)
## dat.full.1BKD.r <- read_stats_table('1BKD', 'r', 5)
## dat.full.1ATN.l <- read_stats_table('1ATN', 'l', 5)
## dat.full.1ATN.r <- read_stats_table('1ATN', 'r', 5)
## dat.full.1ATN.18.l <- read_stats_table('1ATN', 'l', 18)
## dat.full.1ATN.18.r <- read_stats_table('1ATN', 'r', 18)
## dat.full.1F6M.l <- read_stats_table('1F6M', 'l', 5)
## dat.full.1F6M.r <- read_stats_table('1F6M', 'r', 5)
## dat.full <- rbind(dat.full.1F6M.r, dat.full.1F6M.l, dat.full.1ATN.l, dat.full.1ATN.r, dat.full.1BKD.l, dat.full.1BKD.r, dat.full.1ATN.18.r, dat.full.1ATN.18.l)
dat.full <- NULL
#for (s in c(-1, 5, 50)) {
for (s in c(50)) {
  for (i in c(5, 10, 15, 20, 25, 30)) {
    for (lr in c('l', 'r')) {
      for (type in c('halton', 'prng')) {
        print(paste('reading', s, i))
        if (is.null(dat.full)) {
          #dat.full <- read_stats_table('1ATN', 'r', i)
          dat.full <- read_new_stats('1ATN', lr, i, type, s)
        } else {
          #dat.full <- rbind(dat.full, read_stats_table('1ATN', 'r', i))
          dat.full <- rbind(dat.full, read_new_stats('1ATN', lr, i, type, s))
        }
      }
    }
  }
}
dofs <- data.frame(hp=c(5, 10, 15, 20, 25, 30),
                   dofs=c(21, 44, 78, 110, 137, 153))
dat.full <- merge(dat.full, dofs)

add_swapped <- function(dd) {
  dswap <- data.frame(i=dd$j, j=dd$i,
                      ii=dd$jj, jj=dd$ii,
                      i.energy=dd$j.energy, j.energy=dd$i.energy,
                      i.crmsd=dd$j.crmsd, j.crmsd=dd$i.crmsd,
                      dist=dd$dist, # same distance
                      d.energy=-dd$d.energy, # negated distances here
                      d.crmsd=-dd$d.crmsd,
                      type.1=dd$type.2, type.2=dd$type.1,
                      pdb=dd$pdb, lr=dd$lr, hp=dd$hp, dofs=dd$dofs)
  if ('types' %in% names(dd)) {
    dswap$types <- paste(as.character(dswap$type.1), as.character(dswap$type.2))
  }
  print(sort(names(dd)))
  print(sort(names(dswap)))
  rbind(dd, dswap)
}
print(paste('before swap, rows are', nrow(dat.full)))
dat.full <- add_swapped(dat.full)
print(paste('after swap, rows are', nrow(dat.full)))
dat.full$types <- paste(as.character(dat.full$type.1),as.character(dat.full$type.2))

dat.9 <- add_swapped(dat.full[dat.full$i == 9999 | dat.full$j == 9999,])

# hexbin, not super great.
ggplot(dat.9[dat.9$types %in% c('p p', 'h h') & dat.9$i.energy < 0 & dat.9$j.energy < 0,], aes(j.energy, dist)) + geom_hex() + facet_grid(types ~ pdb + lr, scales='free') + geom_point(data=dat.9[dat.9$i == 9999 & dat.9$types %in% c('p p', 'h h'),], aes(x=i.energy, y=0), color='red', size=2)
ggplot(dat.9[dat.9$types %in% c('p p', 'h h') & dat.9$i.energy < 0 & dat.9$j.energy < 0,], aes(j.energy, dist)) + geom_hex() + facet_grid(types ~ pdb + lr, scales='free') + geom_vline(data=dat.9[dat.9$i == 9999 & dat.9$types %in% c('p p', 'h h'),], aes(xintercept=i.energy), color='red', lty=2)
# Not a bad hexbin
ggplot(dat.9[dat.9$i == 9999 & dat.9$j.energy < 1e13 & dat.9$types %in% c('p p', 'h h'),], aes(j.energy, dist)) + geom_hex(bins=30) + facet_grid(types ~ pdb + lr, scales='free') + geom_point(data=dat.9[dat.9$i==9999 & dat.9$types %in% c('p p', 'h h'),], aes(x=i.energy, y=0), color='red', size=2) + scale_fill_distiller(palette='Spectral', direction=-1)
ggplot(dat.9[dat.9$i == 9999 & dat.9$j.energy < 1e5 & dat.9$j.energy > -7e6 & dat.9$types %in% c('p p', 'h h'),], aes(j.energy, dist)) + geom_hex(bins=30) + facet_grid(hp+types ~ pdb + lr, scales='free_x') + geom_vline(data=dat.9[dat.9$i==9999 & dat.9$types %in% c('p p', 'h h'),], aes(xintercept=i.energy), color='red', lty=2) + scale_fill_distiller(palette='Spectral', direction=-1)

# density is more clear
ggplot(dat.9[dat.9$types %in% c('p p', 'h h') & dat.9$i.energy < 0 & dat.9$j.energy < 0,], aes(j.energy, dist, color=types)) + geom_density2d() + facet_wrap(pdb ~ lr, scales='free') + geom_point(data=dat.9[dat.9$i == 9999,], aes(x=i.energy, y=0), color='blue', size=2)



# Either specify a data frame for `dd` or a single column.
# If you do the former, you must also supply pdb and lr.
num_outside <- function(dd, pdb=NULL, lr=NULL, m_iqr=5) {
  if (!is.null(pdb)) {
    ecol <- dd[dd$j==0 & dd$pdb==pdb & dd$lr == lr & dd$types %in% c('p p', 'h h'),'i.energy']
  } else {
    ecol = dd
  }
  Q <- quantile(ecol)
  iqr <- IQR(ecol)
  n_over <- sum(ecol > Q[4]+m_iqr*iqr)
  n_under <- sum(ecol < Q[2]-m_iqr*iqr)
  #data.frame(n_under=n_under, n_over=n_over, under=Q[2]+m_iqr*iqr, over=Q[4]+m_iqr*iqr);
  c(Q[2]-m_iqr*iqr, Q[4]+m_iqr*iqr, n_over, n_under)
}


# Let's try and compute the max/min from every point to every other point.
dfapply <- function(c, f) {
  do.call(rbind, lapply(c, f))
}

dat.by.i <- dfapply(unique(dat.full$pdb), function(pdb) {
  dh.pdb <- dat.full[dat.full$pdb == pdb & dat.full$i.energy < 0 & dat.full$j.energy < 0,];
  dfapply(unique(dh.pdb$lr), function(lr) {
    dh.lr <- dh.pdb[dh.pdb$lr == lr,];
    dfapply(unique(dh.lr$types), function(type) {
      print(paste(pdb, lr, type))
      dh.type <- dh.lr[dh.lr$types == type,];
      # Remove outliers
      outs <- num_outside(dh.type[dh.type$j == 0, 'i.energy'], m_iqr=5)
      dh.type <- dh.type[dh.type$i.energy > outs[1] & dh.type$i.energy < outs[2] &
                         dh.type$j.energy > outs[1] & dh.type$j.energy < outs[2],]
      dfapply(unique(dh.type$i), function(i) {
        data.frame(dist=c(min(abs(dh.type[dh.type$i == i, 'dist'])),
                          max(abs(dh.type[dh.type$i == i, 'dist']))),
                   energy=c(min(abs(dh.type[dh.type$i == i, 'd.energy'])),
                            max(abs(dh.type[dh.type$i == i, 'd.energy']))),
                   crmsd=c(min(abs(dh.type[dh.type$i == i, 'd.crmsd'])),
                           max(abs(dh.type[dh.type$i == i, 'd.crmsd']))),
                   which=c('min', 'max'),
                   dist.9=rep(dh.type[dh.type$i == i & dh.type$j == 9999, 'dist'], 2),
                   d.energy.9=rep(dh.type[dh.type$i == i & dh.type$j == 9999, 'd.energy'], 2),
                   d.crmsd.9=rep(dh.type[dh.type$i == i & dh.type$j == 9999, 'd.crmsd'], 2),
                   crmsd=rep(unique(dh.type[dh.type$i == i, 'i.crmsd'], 2)),
                   energy=rep(unique(dh.type[dh.type$i == i, 'i.energy'], 2)),
                   types=c(type, type), lr=c(lr, lr), pdb=c(pdb, pdb), i=c(i, i))
      })
    })
  })
})
hp_limiters <- c(5, 15, 30)
dat.by.i.maxmin <- dfapply(unique(dat.full$pdb), function(pdb) {
  dh.pdb <- add_swapped(dat.full[dat.full$pdb == pdb &
                                 dat.full$hp %in% hp_limiters &
                                 dat.full$types %in% c('p p', 'h h'),])
  print(paste("--", pdb, nrow(dh.pdb), "--"))
  dfapply(unique(dh.pdb$hp), function(hp) {
  dfapply(unique(dh.pdb$lr), function(lr) {
    dh.lr <- dh.pdb[dh.pdb$lr == lr & dh.pdb$hp == hp,];
    dfapply(unique(dh.lr$types), function(type) {
      dh.type <- dh.lr[dh.lr$types == type,];
      # Remove outliers
      outs <- num_outside(dh.type[dh.type$j == 0, 'i.energy'], m_iqr=5)
      print(paste(pdb, lr, type, hp, 'rows', nrow(dh.type), 'outliers', paste(outs, collapse=' ')))
      dh.type <- dh.type[dh.type$i.energy > outs[1] & dh.type$i.energy < outs[2] &
                         dh.type$j.energy > outs[1] & dh.type$j.energy < outs[2],]
      dfapply(unique(dh.type$i), function(i) {
        data.frame(dist.min=min(abs(dh.type[dh.type$i == i, 'dist'])),
                   dist.max=max(abs(dh.type[dh.type$i == i, 'dist'])),
                   energy.min=min(abs(dh.type[dh.type$i == i, 'd.energy'])),
                   energy.max=max(abs(dh.type[dh.type$i == i, 'd.energy'])),
                   crmsd.min=min(abs(dh.type[dh.type$i == i, 'd.crmsd'])),
                   crmsd.max=max(abs(dh.type[dh.type$i == i, 'd.crmsd'])),
                   dist.9=ifelse(i == 9999, 0, dh.type[dh.type$i == i & dh.type$j == 9999, 'dist']),
                   d.energy.9=ifelse(i == 9999, 0, dh.type[dh.type$i == i & dh.type$j == 9999, 'd.energy']),
                   crmsd=unique(dh.type[dh.type$i == i, 'i.crmsd']),
                   energy=unique(dh.type[dh.type$i == i, 'i.energy']),
                   types=type, lr=lr, pdb=pdb, i=i, hp=hp)
      })
    })
  })
  })
})
dat.by.i.maxmin <- merge(dat.by.i.maxmin, dofs)
dat.by.i.maxmin$ltypes <- ifelse(dat.by.i.maxmin$types == 'p p', 'prng', 'halton')

e_boltzman <- function(df, pdb, lr, type) {
  k <- 1.380649e-23
  T <- 300
  en <- df[df$j == 0 & 
           df$pdb == pdb & df$lr == lr & df$types == type, 'i.energy']
  sum(exp(-en / k / T))
}

# Then can do the following plot:
ggplot(dat.by.i.maxmin[dat.by.i.maxmin$types %in% c('p p', 'h h'),], aes(log(energy.min), energy.max, color=types)) + geom_density2d() + geom_point(alpha=0.2) + facet_wrap(lr~pdb, scales='free')
# Actually, this is pretty good... shows high spread for halton at small
# dimensions, but higher for prng at higher dimensions. (note: this is just for
# 1ATN)
ggplot(dat.by.i.maxmin[dat.by.i.maxmin$types %in% c('p p', 'h h'),], aes(energy.min, energy.max, color=types)) + geom_density2d() + geom_point(alpha=0.2) + facet_wrap(hp+lr~pdb, scales='free') + scale_x_log10() + geom_point(data=dat.by.i.maxmin[dat.by.i.maxmin$i==9999 & dat.by.i.maxmin$types %in% c('p p', 'h h'),], pch=17, size=5)
q <- ggplot(dat.by.i.maxmin[dat.by.i.maxmin$types %in% c('p p', 'h h') & dat.by.i.maxmin$lr=='r',], aes(energy.min/1000, energy.max/1000, color=types)) + geom_density2d() + geom_point(alpha=0.2) + facet_wrap(.~dofs, scales='free') + ylab('max energy diff (kJ)') + xlab('min energy diff (kJ)') + scale_x_log10() + labs(color='Generator') + scale_color_discrete(labels=c('prng', 'halton'))
ggsave("discr_min_v_max_energy_1ATN.pdf", q, width=9, height=5)
dof_lab <- function(string) {
  paste("dim:", string)
}
q <- ggplot(dat.by.i.maxmin[dat.by.i.maxmin$types %in% c('p p', 'h h') & dat.by.i.maxmin$hp %in% c(5, 15, 30),], aes(dist.max, energy.max/1000, color=ltypes)) + geom_density2d() + geom_point(alpha=0.2) + facet_wrap(.~dofs, scales='free_x', labeller=labeller(dofs=dof_lab)) + ylab('max energy dist (kJ)') + xlab('max L2 norm') + labs(color='Generator')
ggsave('discr_dist_energy_1ATN.pdf', q, width=10, height=4)
library(dplyr)
dat.by.i.maxmin %>% group_by(hp, ltypes) %>% summarize(mean_en=mean(energy), sd_en=sd(energy), mean_max_en=mean(energy.max), sd_max_en=sd(energy.max), mean_dist=mean(dist.9), sd_dist=sd(dist.9), mean_max_dist=mean(dist.max), sd_max_dist=sd(dist.max))


df_stats <- data.frame(lr=character(), hp=numeric(), type=character(),
                       energy.max=numeric(), energy.mean=numeric(), energy.min=numeric(),
                       crmsd.max=numeric(), crmsd.mean=numeric(), crmsd.min=numeric(),
                       dist.max=numeric(), dist.mean=numeric(), dist.min=numeric())
# Only applies for 1ATN.
for (hp in c(5, 18)) {
  for (lr in c('l', 'r')) {
    for (type in c('p p', 'h h')) {
      df.here <- dat.by.i.maxmin[dat.by.i.maxmin$lr==lr &
                                 dat.by.i.maxmin$hp==hp &
                                 dat.by.i.maxmin$types==type,]
      print(paste(hp, lr, type))
      df_stats.here <- data.frame(lr=lr, hp=hp, type=type,
                                  energy.max=max(df.here$energy.max, na.rm=T),
                                  energy.min=min(df.here$energy.min, na.rm=T),
                                  energy.mean=mean(
                                                   dat.full[dat.full$pdb=='1ATN' &
                                                            dat.full$hp==hp &
                                                            dat.full$lr==lr &
                                                            dat.full$types==type,'d.energy'], na.rm=T),
                                  crmsd.max=max(df.here$crmsd.max, na.rm=T),
                                  crmsd.min=min(df.here$crmsd.min, na.rm=T),
                                  crmsd.mean=mean(
                                                  dat.full[dat.full$pdb=='1ATN' &
                                                           dat.full$hp==hp &
                                                           dat.full$lr==lr &
                                                           dat.full$types==type,'d.crmsd'], na.rm=T),
                                  dist.max=max(df.here$dist.max, na.rm=T),
                                  dist.min=min(df.here$dist.min, na.rm=T),
                                  dist.mean=mean(
                                                 dat.full[dat.full$pdb=='1ATN' &
                                                          dat.full$hp==hp &
                                                          dat.full$lr==lr &
                                                          dat.full$types==type,'dist'], na.rm=T))
      df_stats <- rbind(df_stats, df_stats.here)
    }
  }
}
print(df_stats)


# Another statistic we've made up. It's kind of like dispersion, but normalized
# by the variance in each direction.
new_disp_funs <- function() {
  # Let's not keep the whole thing, just compute it on the fly. Takes too much
  # memory otherwise.
  dat.half <- unique(dat.9[dat.9$j==9999,c('pdb','dofs','lr','i','i.energy','i.crmsd','hp','type.1','dist')])
  write.table(dat.half, file='1ATN-full-half.stats', row.names=F)

  dat <- read.table('1ATN-full-half.stats', h=T)

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

  new_disc <- function(dat) {
    startt = Sys.time()
    mydat <- dat
    mydat$i.crmsd  = mydat$i.crmsd  / sd(mydat$i.crmsd)
    mydat$i.energy = mydat$i.energy / sd(mydat$i.energy)
    mydat$dist     = mydat$dist     / sd(mydat$dist)
    min(sapply(1:(nrow(mydat)-1), function(i) {
                 if (i %% 100 == 0) {print(paste(i, '/', nrow(mydat), '- time', Sys.time() - startt))}
                 qmn <- min(sapply((i+1):nrow(mydat), function(j) {
                         qm <- max(sqrt((mydat[i,'i.crmsd' ]-mydat[j,'i.crmsd' ])**2),
                                   sqrt((mydat[i,'dist'    ]-mydat[j,'dist'    ])**2),
                                   sqrt((mydat[i,'i.energy']-mydat[j,'i.energy'])**2), na.rm=T)
                        }) , na.rm=T)
                 #print(paste(i, qmn))
                 qmn
                       }))
  }

  qq <- do.call(rbind, apply(unique(dat[,c('pdb','dofs','lr','type.1')]), 1, function(x) {
    print(x)
    d <- remove_outliers(dat[dat$pdb==x[1] & dat$dofs==trimws(x[2]) & dat$lr==x[3] & dat$type.1==x[4],],
                         'i.energy')
    new_disc_d <- new_disc(d)
    data.frame(pdb=x[1], dofs=x[2], lr=x[3], type=x[4],
               disc=new_disc_d, sd.crmsd=sd(d$i.crmsd), sd.energy=sd(d$i.energy), sd.dist=sd(d$dist))
  }))
}

dat.s50 <- NULL
for (h in c(5, 15, 30)) {
  for (type in c('prng', 'halton')) {
    thisdat <- read.table(paste('hp.', h, '.s50-1ATN_r_u-', type, '.prot_en_ll.txt', sep=''))
    names(thisdat) <- c('conf', 'll', 'a','b', 'crmsd', 'ucrmsd', 'amberen', 'cvcen')
    thisdat$hp = h
    thisdat$type = type
    thisdat$s = 50

    if (is.null(dat.s50)) {
      dat.s50 <- thisdat
    } else {
      dat.s50 <- rbind(dat.s50, thisdat)
    }
  }
}
dat.s50 <- merge(dat.s50, dofs)
