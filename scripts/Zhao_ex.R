pacman::p_load(
  tidyverse
  , sf
  , nimble
  , rstan
  #, 
  #, 
)

ex <- "https://raw.githubusercontent.com/QingZhaoQE/MEE-IMBCR-BBS-eBird-data-integration/refs/heads/main/"

imbcr <- read_csv(paste0(ex, "data_1.%20imbcr.csv"))
bbs <- read_csv(paste0(ex, "data_2.%20bbs.csv"))
ebird <- read_csv(paste0(ex, "data_3.%20ebird.csv"))
gra <- read_csv(paste0(ex, "data_x1.%20grass.csv"))
tem <- read_csv(paste0(ex, "data_x2.%20temperature.csv"))

write_csv(imbcr, "data/ready/imbcr.csv")
write_csv(bbs, "data/ready/bbs.csv")
write_csv(ebd, "data/ready/ebird.csv")
write_csv(gra, "data/ready/grass.csv")
write_csv(tem, "data/ready/temperature.csv")

### Code from Zhao's github repo: 
### Standardize covariates
x1 <- as.matrix(gra[,-c(1,2)])
x2 <- as.matrix(tem[,-c(1,2)])
x1 <- (x1 - 0.5) * 4
x2 <- (x2 - mean(x2)) / sd(x2)
gra[,-c(1,2)] <- x1
tem[,-c(1,2)] <- x2

### Filter 0's in eBird data
efish1 <- ebird$fish_id[which(ebird$count > 0)]
efish0 <- ebird$fish_id[which(ebird$count == 0)]
efish0 <- efish0[which(!(efish0 %in% efish1))]
ebird1 <- ebird[which(ebird$fish_id %in% efish1),]
ebird0 <- ebird[which(ebird$fish_id %in% efish0),]
t1 <- unique(ebird0[,c('fish_id','year')])
t2 <- by(t1$year, t1$fish_id, length)
t2 <- data.frame(names(t2), cbind(t2))
names(t2) <- c('fish_id', 'nyear')
t3 <- t2[which(t2$nyear >= 11),]
efish_keep <- sort(unique(c(ebird1$fish_id, t3$fish_id)))
ebird <- ebird[which(ebird$fish_id %in% efish_keep),]
ebird1 <- ebird[which(ebird$count > 0),]
ebird0 <- ebird[which(ebird$count == 0),]
eid <- paste(ebird0$fish_id, ebird0$year, sep='-')
eid <- unique(eid)
for (i in 1:length(eid)) {
  fish_t <- strsplit(eid[i], split='-')[[1]][1]
  year_t <- strsplit(eid[i], split='-')[[1]][2]
  if (length(which(ebird0$fish_id == fish_t & ebird0$year == year_t)) == 1) {
    tt <- ebird0[which(ebird0$fish_id == fish_t & ebird0$year == year_t),]
  } else {
    t1 <- ebird0[which(ebird0$fish_id == fish_t & ebird0$year == year_t),]
    ss <- sample(1:dim(t1)[1], 1)
    tt <- t1[ss,]
  }
  if (i == 1) {
    out <- tt
  } else {
    out <- rbind(out, tt)
  }
} #i
ebird <- rbind(ebird1, out)

### Prepare data
gra <- gra[which(gra$fish_id %in% sort(unique(c(imbcr$fish_id, bbs$fish_id, ebird$fish_id)))),]
tem <- tem[which(tem$fish_id %in% sort(unique(c(imbcr$fish_id, bbs$fish_id, ebird$fish_id)))),]

info <- data.frame(fish_id = gra$fish_id, site = 1:dim(gra)[1])
gra <- as.matrix(gra[,-c(1:2)])
tem <- as.matrix(tem[,-c(1:2)])
nsite <- dim(gra)[1]
nyear <- dim(gra)[2]

imbcr_t <- merge(info, imbcr)
imbcr_t <- imbcr_t[order(imbcr_t$year, imbcr_t$site),]
nobs_imbcr <- dim(imbcr_t)[1]
years_imbcr <- imbcr_t$year
sites_imbcr <- imbcr_t$site
imbcr_cnt <- as.matrix(imbcr_t[,-c(1:4)])
imbcr_cnt_sum <- rowSums(imbcr_cnt)
ndist <- 20    # number of distance bins
ntime <- 6    # number of time bins
cutoff <- 200 # maximum distance for distance sampling
breaks <- seq(0, cutoff, length.out=ndist + 1)

bbs_t <- merge(info, bbs)
bbs_t <- bbs_t[order(bbs_t$year, bbs_t$site),]
nobs_bbs <- dim(bbs_t)[1]
years_bbs <- bbs_t$year
sites_bbs <- bbs_t$site
route_bbs <- as.numeric(as.factor(bbs_t$route))
nroute_bbs <- length(unique(route_bbs))
obser_bbs <- as.numeric(as.factor(bbs_t$obs))
nobser_bbs <- length(unique(obser_bbs))
new_bbs <- bbs_t$new
bbs_cnt_ori <- as.matrix(bbs_t[,-c(1:6)])
bbs_cnt <- matrix(0, nobs_bbs, 10)
for (i in 1:10) {
  bbs_cnt[,i] <- rowSums(bbs_cnt_ori[,(i-1)*5+1:5])
} # i
nstop_bbs <- dim(bbs_cnt)[2]
stops_bbs <- seq(-2, 2, length.out=nstop_bbs)

ebird_t <- merge(info, ebird)
ebird_t <- ebird_t[order(ebird_t$year, ebird_t$site),]
nobs_ebird <- dim(ebird_t)[1]
years_ebird <- ebird_t$year
sites_ebird <- ebird_t$site
locat_ebird <- as.numeric(as.factor(ebird_t$Location))
nlocat_ebird <- length(unique(locat_ebird))
obser_ebird <- as.numeric(as.factor(ebird_t$observer))
nobser_ebird <- length(unique(obser_ebird))
type_ebird <- ifelse(ebird_t$type == 'Stationary', 0, 1)
duration_ebird <- log(ebird_t$duration)
distance_ebird <- ebird_t$distance
n_obsver_ebird <- log(ebird_t$n_obsver)
ebird_cnt <- ebird_t$count

year0 <- min(c(years_imbcr, years_bbs, years_ebird))
years_imbcr <- years_imbcr - year0 + 1
years_bbs <- years_bbs - year0 + 1
years_ebird <- years_ebird - year0 + 1

