# distinguishability analysis

library(tidyverse)
library(ape)
library(phangorn)
library(lubridate)
library(treedater)
library(Biostrings)
library(patchwork)
library(mgcv)
library(pROC)

# load in and process data by serotype
for (k in 1:3){

  # nt alignment
  aln <- read.dna(paste("Sabin",k,".ntP1.alignment.final.phylip.reduced",sep=''), format = "sequential")
  
  # select complete VP1 with no indels only
  if (k==1){
    gnl=906
    aln = aln[,2643+seq(-gnl+1,0,by=1)]
  } else if (k==2){
    gnl=903
    aln = aln[,2637+seq(-gnl+1,0,by=1)]
  } else if (k==3){
    gnl=900
    aln = aln[,2634+seq(-gnl+1,0,by=1)]
  }
  missing_idx = apply(!as.character(aln)!='-',1,any)
  which(missing_idx)
  
  # filter to complete VP1 with no indels
  aln = aln[!missing_idx,]  
  
  # metadata
  df = read.csv(paste('Sabin',k,'.ntP1.alignment.final.attributes.csv',sep=''))
  df$classification = sub('Sabin[1-3]','Sabin',df$classification) # remove the type label
  df = df |> filter(accession %in% labels(aln))
  
  # filter to cVDPV, WPV, iVDPV, SabinLike, VDPV, Sabin only
  keep = df |> filter(classification %in% c('cVDPV','iVDPV','wild', 'SabinLike','VDPV','Sabin')) |> select(accession)
  df = df |> filter(accession %in% keep$accession)
  
  aln = aln[labels(aln) %in% keep$accession,] 
  
  # translate AA
  aln_aa = as.AAbin(translate(DNAStringSet(apply(as.character(aln), 1, paste0, collapse = "")), 
                              if.fuzzy.codon = "X"))
  
  # counts different from Sabin
  for (CLASS in c('cVDPV','iVDPV','wild', 'SabinLike','VDPV')){
    keep = df |> filter(classification == CLASS | classification=='Sabin') |> select(accession)
    
    nt_d <- as.matrix(dist.ml(aln[labels(aln) %in% keep$accession,] , model = "F81"))[1,] * gnl
    aa_d <- as.matrix(dist.ml(aln_aa[labels(aln_aa) %in% keep$accession] , model = "FLU"))[1,]*gnl/3
    ns_d = (nt_d-aa_d)
    ns_d = ns_d - min(ns_d) # clean up mle slop
    
    if (CLASS == 'iVDPV'){
      # rescale iVDPV clock
      # https://www.sciencedirect.com/science/article/abs/pii/S0264410X15000377
      dt = nt_d/(gnl*0.0072)
    } else {
      # https://journals.asm.org/doi/full/10.1128/jvi.02354-07
      dt = ns_d/(gnl*0.01)
    }

    if (CLASS == 'cVDPV'){
      tmp = data.frame(accession = names(dt),
                       classification=c('Sabin',rep(CLASS,length(dt)-1)),
                       serotype = k,
                       years = dt, 
                       nt_diff=nt_d, 
                       d_n = aa_d,
                       d_s=ns_d)
    } else {
      tmp = tmp |> rbind(data.frame(accession = names(dt),
                                    classification=c('Sabin',rep(CLASS,length(dt)-1)),
                                    serotype = k,
                                    years = dt, 
                                    nt_diff=nt_d, 
                                    d_n = aa_d,
                                    d_s=ns_d) |> filter(classification!='Sabin'))
    }
  }
  # join
  df = df |> left_join(tmp,by=c('accession','classification'))
  
  
  # clean up classification
    mislabeled=df |> filter(nt_diff==0 & !grepl('Sabin',classification)) |> select(accession)
    df$classification[df$accession %in% mislabeled$accession] = 'SabinLike'
  
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC165252/
    mislabeled=df |> filter(classification=='cVDPV' & years>12) |> select(accession)
    df$classification[df$accession %in% mislabeled$accession] = 'wild'
    
    # https://www.ncbi.nlm.nih.gov/nuccore/KC769998
    mislabeled=df |> filter(classification=='VDPV' & d_s>200) |> select(accession)
    df$classification[df$accession %in% mislabeled$accession] = 'wild'
    
  # put all the data together
  if (k==1){
    meta = df
  } else {
    meta = meta |> rbind(df)
  }
}

# organize factor levels
meta = meta |> 
  filter(classification!='Sabin') |> 
  mutate(serotype=factor(serotype, levels=c('1','2','3'))) |>
  mutate(classification = factor(classification,
                                 levels=c('SabinLike','VDPV','cVDPV','iVDPV','wild')))

ggplot(meta) +
  geom_histogram(aes(x=years,fill=classification),binwidth = 1, boundary=0)

ggplot(meta) +
  geom_point(aes(x=years,y=d_n,color=classification))

p1=ggplot(meta,aes(x=d_s,y=d_n,color=classification,shape=serotype)) +
  geom_jitter(width=0.3,height=0.3,alpha=0.2) +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  guides(shape=guide_legend('serotype'))
p1

df = meta |> filter(classification %in% c('cVDPV','iVDPV'))
p2=ggplot(df,aes(x=years,y=d_n,color=classification,fill=classification)) +
  geom_jitter(aes(shape=serotype),width=0.3,height=0.3,alpha=0.2) +
  geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('approx. years since OPV') +
  scale_x_continuous(breaks=seq(0,30,by=5)) +
  guides(shape=guide_legend('serotype'))
p2

p3=ggplot(df,aes(x=d_s,y=d_n,color=classification,fill=classification)) +
  geom_jitter(aes(shape=as.character(serotype)),width=0.3,height=0.3,alpha=0.2) +
  geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  guides(shape=guide_legend('serotype'))
p3

p2 + p3 + p1
ggsave('iVDPV_vs_cVDPV.png',units='in',width=12,height=3.5)


# check consistency with Zhao 2017 fig 2  https://www.csbj.org/article/S2001-0370(17)30081-8/fulltext
# this is comparing slightly different things (distance from sabin vs distance between isolates), 
# but it largely checks out
zhao = read.csv('zhao_fig2_extracted_vpdv_points.csv')
names(zhao) = c('nt_diff','d_n','classification')
zhao = zhao |>
  mutate(nt_diff = (1-nt_diff/100)*903,
         d_n = (1-d_n/100)*301) |>
  mutate(d_s = nt_diff-d_n) |>
  mutate(d_s = d_s - min(d_s))

ggplot() +
  geom_point(data=df,aes(x=1-nt_diff/903,y=1-d_n/301)) +
  geom_point(data=zhao,aes(x=1-nt_diff/903,y=1-d_n/301,color=classification)) 

# checking the dating
# dates of late iVPDVs tend to be off by a few years, 
# generally reflecting under-dispersion relative to the true individual-case variability.
# I think this is good enough for our purposes, but peer review would ask for more careful dating.
df |> filter(classification=='iVDPV' & years>8 ) |>
  group_by(paper) |>
  summarize(n=n(),
            accessions = paste(unique(accession),collapse=":")) |>
  arrange(desc(n))

df |> filter(classification=='iVDPV' & years>8) |>
  ggplot(aes(x=years,y=d_n,color=paper)) +
  geom_point() 


## classifier models

# cVDPV vs iVDPV
df = meta |> filter(classification %in% c('cVDPV','iVDPV')) |> droplevels()
df |> group_by(serotype,classification) |> summarize(n=n())

names(meta)

mod1 = gam(formula = classification ~ d_n + d_s,
           data =df,
           family='binomial')
mod2 = gam(formula = classification ~ d_n + d_s + serotype,
           data =df,
           family='binomial')
summary(mod2)
# mod1 is good enough since no sig diff between types 1 and 2 (and no data for 3)
summary(mod1)
plot(mod1)

# nonlinear model is better than linear model
mod3 = gam(formula = classification ~ s(d_n) + s(d_s),
           data =df,
           family='binomial')
summary(mod3)
c(AIC(mod1),AIC(mod3))
plot(mod3)


# interaction model looks pretty overfit
mod4 = gam(formula = classification ~ ti(d_n,d_s),
           data =df,
           family='binomial')
summary(mod4)
c(AIC(mod3),AIC(mod4))
plot(mod4)

# mod3=mod1 # wow mod1 is way worse, as one can see looking at the classifier plot p6 below
df$iVDPVness = predict(mod3,newdata = df,type='response')

# ROC analysis
mod3_roc = roc(df$classification,df$iVDPVness)
# best thresholds
  youden = coords(mod3_roc,'best',best.method	= 'youden')
  youden

  # weighted = coords(mod3_roc,'best',best.method	= 'youden',best.weights=c(0.1,1e-7))
  weighted = coords(mod3_roc,'best',best.method	= 'youden',best.weights=c(0.1,0.5))
  weighted

opt_df = data.frame(specificity = c(youden$specificity,weighted$specificity),
                    sensitivity = c(youden$sensitivity,weighted$sensitivity),
                    threshold = c(youden$threshold,weighted$threshold),
                    cost=c('equal','c/i = 10x'))
p4=ggroc(mod3_roc) + 
  geom_point(data=opt_df,aes(x=specificity,y=sensitivity,color=cost)) +
  theme_bw() +
  annotate('text',x=0.2,y=0.9,label=paste('auc=',round(auc(df$classification,df$iVDPVness),2),sep=''))

p5=ggplot() +
  stat_ecdf(data=df,aes(x=iVDPVness,color=classification)) +
  geom_hline(data=opt_df,aes(yintercept = threshold, linetype=cost)) +
  theme_bw() +
  guides()

p4+p5

# level sets

# fit is your mgcv::gam with z ~ s(x) + s(y)
  
  # 1) grid over the domain
  d_s_range = seq(0,140,by=1)
  d_n_range = seq(0,50,by=1)
  grid <- expand.grid(d_s = d_s_range, d_n = d_n_range)
  
  # 2) predict on the desired scale
  #   - use type="response" if C is on the response scale
  #   - if C is on the link scale, use type="link"
  zmat <- matrix(predict(mod3, newdata = grid, type = "response"),
                 nrow = length(d_s_range), ncol = length(d_n_range))
  
  # 3) extract level set(s) z = C
  cl <- contourLines(x = d_s_range, y = d_n_range, z = zmat, levels = opt_df$threshold)
  
  # cl is a list of components of the level set; each has $x and $y
  # Example: bind all polylines into one data.frame
  lev_df <- do.call(rbind, lapply(seq_along(cl), function(i)
    data.frame(comp = i, d_s = cl[[i]]$x, d_n = cl[[i]]$y)
  ))
  lev_df$cost = unique(opt_df$cost)[lev_df$comp]


p6=ggplot(df,aes(x=d_s,y=d_n)) +
  geom_jitter(aes(color=classification,fill=classification),width=0.3,height=0.3,alpha=0.2) +
  geom_line(data =lev_df,aes(linetype = cost)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) 
p6


p7=ggplot(df,aes(x=d_s,y=d_n)) +
  geom_jitter(aes(color=classification,fill=classification),width=0.3,height=0.3,alpha=0.2) +
  geom_line(data =lev_df,aes(linetype = cost)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) 
p7

(p4 + p5) / p6
ggsave('classifier.png',units='in',width=12,height=3.5)


## classify vdpvs
meta_pred = meta |> filter(classification %in% c('cVDPV','iVDPV','VDPV'))
meta_pred$iVDPVness = predict(mod3,newdata = meta_pred,type='response')

meta_pred |> group_by(classification) |> summarize(mean(iVDPVness))



p4=ggplot(meta_pred,aes(x=d_s,y=d_n,color=classification,fill=classification)) +
  geom_jitter(aes(shape=as.character(serotype)),width=0.3,height=0.3,alpha=0.2) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  guides(shape=guide_legend('serotype'))
p4

p5=ggplot(meta_pred,aes(x=d_s,y=d_n,color=iVDPVness,fill=iVDPVness)) +
  geom_jitter(aes(shape=as.character(serotype)),width=0.3,height=0.3,alpha=0.2) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  guides(shape=guide_legend('serotype'))
p5

p3 + p4 + p5 + plot_layout(guides='collect')
ggsave('iVDPVness.png',units='in',width=18,height=3.5)
