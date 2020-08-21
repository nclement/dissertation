read_hp_stats <- function(pdb, lr, hp, type, s) {
  fn = paste('hp.', hp, '.s', s, '-', pdb, '_', lr, '_u-', type, '.hps', sep='')
  if (file.exists(fn)) {
    dat <- read.table(fn)
    names(dat) <- c('fl', 'pca1', 'pca2', 'pca3')
    dat$lr = lr
    dat$pdb = pdb
    dat$hp = hp
    dat$s = s
    dat$type = type
    dat
  } else {
    NULL
  }
}

get_dat <- function() {
  dat.full <- NULL
  for (s in c(-1, 5, 50)) {
  #for (s in c(50)) {
    for (i in c(5, 10, 15, 20, 25, 30)) {
      for (lr in c('l', 'r')) {
        for (type in c('halton', 'prng')) {
          print(paste('reading', s, i))
          if (is.null(dat.full)) {
            dat.full <- read_hp_stats('1ATN', lr, i, type, s)
          } else {
            dat.full <- rbind(dat.full, read_hp_stats('1ATN', lr, i, type, s))
          }
        }
      }
    }
  }
  dat.full$original = "sampled"
  dat.full[grepl("benchmark5", dat.full$fl), "original"] <- "original"

  dofs <- data.frame(hp=c(5, 10, 15, 20, 25, 30),
                     dofs=c(21, 44, 78, 110, 137, 153))
  dat.full <- merge(dat.full, dofs)

  dat.full
}

plot_funs <- function() {
  do.call(rbind, lapply(unique(dat$hp),
        function(hp) {
          print(hp)
          dtest <- dat[dat$hp==hp & dat$lr=='l' & dat$s==-1,];
          do.call(rbind, lapply(1:nrow(dtest),
                function(i) cbind(dtest[i,],
                                  min=min(dist(i, dtest)),
                                  max=max(dist(i, dtest)))))
        }))
}

