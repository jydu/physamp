# SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <dutheil@evolbio.mpg.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Created on 09/10/14 by jdutheil
library(RColorBrewer)
cols<-brewer.pal(6, name="Accent")

read.profile<-function(file) {
  dat <- read.table(file, header=TRUE)
  dat$relSeq <- cumsum(-dat$DiffSequences) / dat[1, "NbSequences"]
  dat$relSit <- cumsum(dat$DiffSites) / dat[nrow(dat), "NbSites"]
  return(dat)
}

dat<-read.profile("bppalnoptim_fastml.log")
dat.a1<-read.profile("bppalnoptim_auto_complete.log")
dat.a2<-read.profile("bppalnoptim_auto_single.log")
dat.a3<-read.profile("bppalnoptim_auto_average.log")
dat.a4<-read.profile("bppalnoptim_auto_median.log")
dat.a5<-read.profile("bppalnoptim_auto_centroid.log")
dat.a6<-read.profile("bppalnoptim_auto_ward.log")
re<-range(dat$NbSequences, dat.a1$NbSequences)
ri<-range(dat$NbSites, dat.a1$NbSites)
rt<-range(dat$Iteration, dat.a1$Iteration)

par(mar=c(4,4,4,4)+0.1)
plot(dat$NbSequences~dat$Iteration, type="b", col="black", pch=4, ylim=re, xlim=rt, ylab="# sequences", xlab="iteration", lwd=2)
lines((dat$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat$Iteration, type="b", col="black", pch=1, lwd=2)
step<-round(diff(ri)/5)
axis(side = 4, at = seq(re[1], re[2], by=step*diff(re)/diff(ri)), labels = seq(ri[1], ri[2], by=step))
mtext(side = 4, line = 3, text = "# sites")

lines(dat.a1$NbSequences~dat.a1$Iteration, type="b", col=cols[1], pch=4, lwd=2)
lines((dat.a1$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat.a1$Iteration, type="b", col=cols[1], pch=1, lwd=2)

lines(dat.a2$NbSequences~dat.a2$Iteration, type="b", col=cols[2], pch=4, lwd=2)
lines((dat.a2$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat.a2$Iteration, type="b", col=cols[2], pch=1, lwd=2)

lines(dat.a3$NbSequences~dat.a3$Iteration, type="b", col=cols[3], pch=4, lwd=2)
lines((dat.a3$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat.a3$Iteration, type="b", col=cols[3], pch=1, lwd=2)

lines(dat.a4$NbSequences~dat.a4$Iteration, type="b", col=cols[4], pch=4, lwd=2)
lines((dat.a4$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat.a4$Iteration, type="b", col=cols[4], pch=1, lwd=2)

lines(dat.a5$NbSequences~dat.a5$Iteration, type="b", col=cols[5], pch=4, lwd=2)
lines((dat.a5$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat.a5$Iteration, type="b", col=cols[5], pch=1, lwd=2)

lines(dat.a6$NbSequences~dat.a6$Iteration, type="b", col=cols[6], pch=4, lwd=2)
lines((dat.a6$NbSites-ri[1])*diff(re)/diff(ri)+re[1]~dat.a6$Iteration, type="b", col=cols[6], pch=1, lwd=2)

# Tradeoff curves:
myplot<-function(dat, col) {
  lines(relSit~relSeq, dat, pch=19, col=col, lwd=2, type="b")
}
plot(relSit~relSeq, dat, pch=19, col="black", lwd=2, type="b", xlab="% sequences removed", ylab="% site gains", ylim=c(0, 0.4))
myplot(dat.a1, col=cols[1])
myplot(dat.a2, col=cols[2])
myplot(dat.a3, col=cols[3])
myplot(dat.a4, col=cols[4])
myplot(dat.a5, col=cols[5])
myplot(dat.a6, col=cols[6])
legend(0, 0.4, col=c("black", cols), pch=19, lwd=2, lty="solid", legend=c("FastTree", "complete", "single", "average", "median", "centroid", "Ward"))

