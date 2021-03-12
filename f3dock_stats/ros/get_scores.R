require(ggplot2)

load_file <- function(score_fn, rmsd_fn, pdb, tpy) {
  print(score_fn)
  dat.pdb.scores <- read.table(score_fn)
  print(rmsd_fn)
  dat.pdb.rmsd <- read.table(rmsd_fn)
  names(dat.pdb.scores) <- c('score', 'init_rms', 'lig_conf', 'rec_conf', 'conf')
  dat.pdb.scores$lig_conf.n <- as.integer(dat.pdb.scores$lig_conf)
  dat.pdb.scores$rec_conf.n <- as.integer(dat.pdb.scores$rec_conf)
  dat.pdb.scores$score <- as.numeric(dat.pdb.scores$score)
  dat.pdb.scores$rank <- rank(dat.pdb.scores$score)
  names(dat.pdb.rmsd) <- c('conf', 'rmsd')
  dat.pdb.rmsd$conf <- gsub(paste("output_files_", pdb, "_out_tog", sep=""),
                            paste(pdb, ".out.tog", sep=""),
                            dat.pdb.rmsd$conf)

  dat.pdb <- merge(dat.pdb.scores, dat.pdb.rmsd)
  dat.pdb$conf_num <- as.numeric(gsub(".*.tog_", "", dat.pdb$conf))
  dat.pdb$pdb <- pdb
  dat.pdb$type <- tpy
  print(summary(dat.pdb))
  dat.pdb
}

load_dat <- function() {
  datl <- list()

  for (pdb in c('1BKD', '1ATN', '1H1V', '1RKE', '1Y64', '1ZLI')) {
    for (tpy in c('atom', 'f3d')) {
      score_fn = paste(pdb, '.out.', tpy, '.scores.txt', sep='')
      rmsd_fn = paste(pdb, '.out.', tpy, '.rmsd.txt', sep='')
      dat.pdb <- load_file(score_fn, rmsd_fn, pdb, tpy)
      # dat.pdb.scores <- read.table(paste(pdb, '.out.', tpy, '.scores.txt', sep=''))
      # dat.pdb.rmsd <- read.table(paste(pdb, '.out.', tpy, '.rmsd.txt', sep=''))
      # names(dat.pdb.scores) <- c('score', 'init_rms', 'lig_conf', 'rec_conf', 'conf')
      # dat.pdb.scores$lig_conf.n <- as.integer(dat.pdb.scores$lig_conf)
      # dat.pdb.scores$rec_conf.n <- as.integer(dat.pdb.scores$rec_conf)
      # dat.pdb.scores$score <- as.numeric(dat.pdb.scores$score)
      # dat.pdb.scores$rank <- rank(dat.pdb.scores$score)
      # names(dat.pdb.rmsd) <- c('conf', 'rmsd')
      # dat.pdb.rmsd$conf <- gsub(paste("output_files_", pdb, "_out_tog", sep=""),
      #                           paste(pdb, ".out.tog", sep=""),
      #                           dat.pdb.rmsd$conf)

      # dat.pdb <- merge(dat.pdb.scores, dat.pdb.rmsd)
      # dat.pdb$conf_num <- as.numeric(gsub(".*.tog_", "", dat.pdb$conf))

      datl[[paste(pdb, tpy)]] <- dat.pdb
    }
  }
  dat.pdb <- load_file("1RKE.out.f3d.2.scores.txt", "1RKE.out.f3d.2.rmsd.txt", '1RKE', 'f3d.2')
  datl[['1RKE f3d.2']] <- dat.pdb

  do.call(rbind, datl)
}

do_stats <- function(dat) {
  require(dplyr)
  dat %>% group_by(pdb, type) %>% mutate(type_rmin=min(rmsd), type_rmax=max(rmsd), type_scmin=min(score), type_scmax=max(score)) %>% ungroup() %>% group_by(pdb, type, lig_conf.n, rec_conf.n) %>% summarize(norm_rmsd=(min(rmsd)-type_rmin)/(type_rmax-type_rmin), m.rmsd=min(rmsd), norm_score=(min(score)-type_scmin)/(type_scmax-type_scmin), m.score=min(score), m.rank=max(rank)) %>% ggplot(aes(lig_conf.n, rec_conf.n, color=norm_score)) + geom_point() + facet_wrap(pdb~type)

}

