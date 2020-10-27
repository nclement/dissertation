require(dplyr)
require(ggplot2)

# Grab the rest of the things.
read_f3d <- function(pdb) {
  require(data.table)

  print(paste('reading scores for', pdb))
  dd <- fread(paste('f3dock_stats/', pdb, '_all_scores.txt', sep=''), h=T)
  # Do some renaming
  setnames(dd, c('class', 'rmsd', 'score_normalized'),
               c('difficulty', 'badRmsd', 'realScore'))
  dd$conf_type = 'vMRSHD'
  dd$pdb = pdb

  print('  reading rmsd')
  dd.rmsd <- fread(paste('f3dock_stats/', pdb, '.rmsd.txt', sep=''))
  names(dd.rmsd) <- c('name', 'rank', 'rmsd')
  dd.rmsd$a = NA

  dd <- merge(dd, dd.rmsd, by=c('name', 'rank'))
  dd
}

load_dat <- function() {
  require(data.table)
  print('reading atom...')
  dat.atom <- fread("all_scores.atom.txt", h=T)
  print(names(dat.atom))
  dat.atom$conf_type = "atomic"
  dat.atom.rmsd <- fread("all_rmsds.atom.txt", h=T)
  dat.atom <- merge(dat.atom, dat.atom.rmsd, by=c('pdb', 'name', 'rank'))

  print('reading hinge...')
  dat.hinge <- fread("all_scores.hinge.txt", h=T)
  dat.hinge$conf_type = "hinge"
  dat.hinge.rmsd <- fread("all_rmsds.hinge.txt", h=T)
  dat.hinge <- merge(dat.hinge, dat.hinge.rmsd, by=c('pdb', 'name', 'rank'))

  print('reading vmrshd...')
  dat.vmrshd <- fread("all_scores.vmrshd.txt", h=T)
  dat.vmrshd$conf_type = "vmrshd"
  dat.vmrshd.rmsd <- fread("all_rmsds.vmrshd.txt", h=T)
  dat.vmrshd <- merge(dat.vmrshd, dat.vmrshd.rmsd, by=c('pdb', 'name', 'rank'))

  dat.full <- rbind(dat.atom, dat.hinge, dat.vmrshd)
  print(unique(dat.full$pdb))
  #dat.full <- dat.full[dat.full$pdb %in% c('3AAD', '1Y64', '2I9B', '2OT3', '1RKE', '1R8S', '1H1V'),]
  #dat.full <- dat.full[dat.full$pdb %in% c('1ATN', '1BKD', '3FN1'),]

  #for (pdb in unique(dat.full$pdb)) {
  #  dd <- read_f3d(pdb)
  #  dat.full <- rbind(dat.full, dd)
  #}
  dat.full
}

lm_things <- function() {
  # Some simple predictions.
  lm_names <- names(dat.full)[!(names(dat.full) %in% c('rank', 'difficulty', 'pdb', 'name', 'conf_type', 'rmsd', 'a', 'badRmsd', 'ssr', 'ccr', 'hbondEn', 'hbond', 'num_hbond', 'deldispe', paste('mat', 1:12, sep='')))]
  dat.atom.lm <- lm(dat.atom, formula=paste("rmsd ~ ", paste(lm_names, collapse=" + ")))
  dat.atom$predict <- predict(dat.atom.lm, dat.atom)
  dat.hinge.lm <- lm(dat.hinge, formula=paste("rmsd ~ ", paste(lm_names, collapse=" + ")))
  dat.hinge$predict <- predict(dat.hinge.lm, dat.hinge)
}

minBy <- function(x, rank, by) {
  sapply(by, function(t) min(x[rank <= t]))
}
minByQ <- function(x, rank, by) {
  tibble(x=minBy(x, rank, by), r=by)
}

plots <- function() {
  dat.full %>% mutate(sampled=ifelse(name=='bound', 'bound', ifelse(name=='unbound', 'unbound', 'sampled'))) %>% group_by(pdb, conf_type, sampled) %>% summarize(rcor=cor(rmsd, rank), min.1=min(rmsd[rank==1]), min.10=min(rmsd[rank<=10]), min.50=min(rmsd[rank<=50]), min.100=min(rmsd[rank<=100]), min.1000=min(rmsd[rank<=1000]), min.full=min(rmsd)) %>% print(n=20)
  dat.full %>% mutate(sampled=ifelse(name=='bound', 'bound', ifelse(name=='unbound', 'unbound', 'sampled'))) %>% group_by(pdb, conf_type, sampled) %>% summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100, 1000))) %>% print(n=20)
  dat.full %>% mutate(sampled=ifelse(name=='bound', 'bound', ifelse(name=='unbound', 'unbound', 'sampled'))) %>% group_by(pdb, conf_type, sampled) %>% summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100, 1000))) %>% ggplot(aes(paste(pdb,sampled), x, col=as.factor(r))) + geom_point() + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept=5, lty=2, color='red') + facet_grid(conf_type~.)

  dat.full %>% mutate(sampled=case_when(name=='bound' ~ 'bound',
                                        name=='unbound' ~ 'unbound',
                                        TRUE ~ 'sampled'),
                      sampled=factor(sampled, levels=c('bound','unbound','sampled'))) %>%
  group_by(pdb, conf_type, sampled) %>%
  summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100, 1000))) %>%
  ggplot(aes(sampled, x, col=as.factor(r))) + geom_point() + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept=5, lty=2, color='red') + facet_grid(conf_type~pdb)

  dat.full %>% mutate(sampled=case_when(name=='bound' ~ 'bound',
                                        name=='unbound' ~ 'unbound',
                                        TRUE ~ conf_type),
                      sampled=factor(sampled,
                                     levels=c('bound','unbound','atomic', 'hinge', 'vMRSHD'))) %>%
                filter(rank < 50) %>% 
                group_by(pdb, sampled) %>%
                summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100, 1000))) %>%
      ggplot(aes(sampled, x, col=as.factor(r))) + geom_point() + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept=5, lty=2, color='red') + facet_grid(.~pdb)


  dat.full %>% mutate(sampled=case_when(name=='bound' ~ 'bound',
                                        name=='unbound' ~ 'unbound',
                                        TRUE ~ conf_type),
                      sampled=factor(sampled,
                                     levels=c('bound','unbound','atomic', 'hinge', 'vMRSHD'))) %>%
        filter(rank < 50) %>% 
        group_by(pdb, sampled, name) %>%
        summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100))) %>%
        ggplot(aes(name, x, col=as.factor(r))) +
            geom_point() + theme(axis.text.x = element_text(angle = 90)) +
            geom_hline(yintercept=5, lty=2, color='red') + facet_grid(sampled~pdb)

  dat.full %>% mutate(sampled=case_when(name=='bound' ~ 'bound',
                                   name=='unbound' ~ 'unbound',
                                   TRUE ~ conf_type),
                 sampled=factor(sampled, 
                                levels=c('bound','unbound','atomic', 'hinge', 'vMRSHD'))) %>%
          filter(rank < 50) %>%
    ggplot(aes(rmsd, score)) + geom_bin2d() + facet_wrap(sampled~pdb, scales='free') + geom_vline(xintercept=5, col='red', lty=2)


  dat.full %>% mutate(sampled=case_when(name=='bound' ~ 'bound',
                                        name=='unbound' ~ 'unbound',
                                        TRUE ~ conf_type),
                      sampled=factor(sampled,
                                     levels=c('bound','unbound','atomic', 'hinge', 'vMRSHD')),
                      name=case_when(name=='bound' ~ paste('bound', conf_type),
                                     name=='unbound' ~ paste('unbound', conf_type),
                                     TRUE ~ name),
                      ) %>%
        filter(pdb %in% c('1Y64', '1H1V', '1BKD', '1R8S')) %>%
        group_by(pdb, sampled, name) %>%
        summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100, 1000))) %>%
        ggplot(aes(as.factor(r), x, fill=sampled)) + geom_boxplot() + facet_grid(.~pdb)

  # Look at individual runs for each different type.
  dat.full %>% mutate(sampled=case_when(name=='bound' ~ 'bound',
                                        name=='unbound' ~ 'unbound',
                                        TRUE ~ conf_type),
                      sampled=factor(sampled,
                                     levels=c('bound','unbound','atomic', 'hinge', 'vMRSHD')),
                      name=case_when(name=='bound' ~ paste('bound', conf_type),
                                     name=='unbound' ~ paste('unbound', conf_type),
                                     TRUE ~ name),
                      ) %>%
        filter(pdb %in% c('1Y64', '1H1V', '1BKD', '1R8S', '2OT3')) %>%
        group_by(pdb, sampled, name) %>%
        summarize(rcor=cor(rmsd, rank), minByQ(rmsd, rank, by=c(1, 5, 100, 1000))) %>%
        ggplot(aes(name, x, col=as.factor(r))) + geom_point() + facet_grid(pdb~sampled, scale='free_y') + geom_hline(yintercept=5, lty=2, col='red')


}
