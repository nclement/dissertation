require(ggplot2)
dat <- rbind(read.table('test_loop_164_200.txt'),
             read.table('test_loop_163_200.txt'),
             read.table('test_loop_164_201.txt'))
names(dat) <- c('flex', 'res', 'phi','psi','nsol')

breaks <- c(-pi, -pi/2, 0, pi/2, pi)
labels <- c(expression(-pi),
            expression(-pi/2),
            expression(0),
            expression(pi/2),
            expression(pi))

for (ress in list(c(164, 200), c(163, 200), c(164, 201))) {
  q <- ggplot(dat[dat$flex==ress[1] & dat$res==ress[2] & dat$nsol > 0,], aes(phi,psi, fill=nsol)) + 
    scale_color_distiller(palette = "Spectral", direction=-1, values=c(0, 0.1, 0.2, 0.3, 0.5, 1)) +
    geom_tile() + 
    theme(legend.position = "none") +
    scale_x_continuous(breaks=breaks, labels=labels, limits=c(-pi, pi)) +
    scale_y_continuous(breaks=breaks, labels=labels, limits=c(-pi, pi)) +
    labs(x=expression(phi), y=expression(psi))

  ggsave(paste('loop_closure', ress[1], ress[2], '.pdf', sep='_'), q, width=6, height=5)
}
