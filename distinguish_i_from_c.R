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
library(ggbeeswarm)


theme_set(theme_bw()) # Optional: Set a default theme


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
    
    # https://www.ncbi.nlm.nih.gov/nuccore/JF313098
    mislabeled=df |> filter(classification=='wild' & years<12) |> select(accession)
    df$classification[df$accession %in% mislabeled$accession] = 'cVDPV'
    
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
                                 levels=c('SabinLike','VDPV','iVDPV','cVDPV','wild')))

meta = meta |> mutate(circulating = if_else(classification %in% c('cVDPV','wild'),'yes',
                                            if_else(classification %in% 'VDPV','?','no'))) |>
  mutate(circulating = factor(circulating, levels=c('no','yes','?')))

cmap = c("#ae7766","#00a1e4","#dc0073")
names(cmap)=levels(meta$circulating)



## plots

p0=ggplot(meta,aes(x=d_s,y=d_n,color=classification,shape=serotype)) +
  geom_point(alpha=0.2) +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') 
p0

p1=ggplot(meta,aes(x=d_s,y=d_n,color=circulating,shape=serotype)) +
  geom_point(alpha=0.2) +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') 
p1


df = meta |> dplyr::filter(circulating %in% c('no','yes')) |> ungroup() |> droplevels()
p2=ggplot(df,
          aes(x=years,y=d_n,color=circulating,fill=circulating)) +
  geom_point(alpha=0.2) +
  geom_smooth(method='gam',fill='black',alpha=0.3) +
  ylab('non-synonymous substitutions') +
  xlab('approx. years since OPV') +
  scale_x_continuous(breaks=seq(0,30,by=5)) +
  scale_color_manual(values=cmap)
p2

p3=ggplot(df,aes(x=d_s,y=d_n,color=circulating,fill=circulating)) +
  geom_jitter(aes(shape=serotype),width=0.3,height=0.3,alpha=0.2) +
  geom_smooth(method='gam',fill='black',alpha=0.3) +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_color_manual(values=cmap)
p3

p2 + p3 + p1
ggsave('iVDPV_vs_cVDPV.png',units='in',width=12,height=3.5)


## classifier models

# cVDPV vs iVDPV
df = meta |> filter(circulating %in% c('no','yes')) |> droplevels()
df |> group_by(serotype,circulating) |> summarize(n=n())

names(meta)

mod1 = gam(formula = circulating ~ d_n + d_s,
           data =df,
           family='binomial')
mod2 = gam(formula = circulating ~ d_n + d_s + serotype,
           data =df,
           family='binomial')
summary(mod2)
# mod1 is good enough since no sig diff between types 1 and 2 (and no data for 3)
summary(mod1)
plot(mod1)

# nonlinear model is better than linear model
mod3 = gam(formula = circulating ~ s(d_n) + s(d_s),
           data =df,
           family='binomial')
summary(mod3)
c(AIC(mod1),AIC(mod3))
plot(mod3)


# interaction model looks pretty overfit
mod4 = gam(formula = circulating ~ ti(d_n,d_s),
           data =df,
           family='binomial')
summary(mod4)
c(AIC(mod3),AIC(mod4))
plot(mod4)

# years instead of d_s is way worse
mod5 = gam(formula = circulating ~ s(d_n) + s(years),
           data =df,
           family='binomial')
summary(mod5)
c(AIC(mod3),AIC(mod5))
plot(mod5)

# mod3=mod1 # wow mod1 is way worse, as one can see looking at the classifier plot p6 below
df$cVDPVness = predict(mod3,newdata = df,type='response')
df$cVDPVness_years = predict(mod5,newdata = df,type='response')

# ROC analysis
mod3_roc = roc(df$circulating,df$cVDPVness)
mod5_roc = roc(df$circulating,df$cVDPVness_years)

# best thresholds
  youden = coords(mod3_roc,'best',best.method	= 'youden')
  youden
  
  youden_years = coords(mod5_roc,'best',best.method	= 'youden')
  youden_years

  # weighted = coords(mod3_roc,'best',best.method	= 'youden',best.weights=c(10,1e-7))
  weighted = coords(mod3_roc,'best',best.method	= 'youden',best.weights=c(10,0.5))
  weighted
  
  weighted_years = coords(mod5_roc,'best',best.method	= 'youden',best.weights=c(10,0.5))
  weighted_years

opt_df = data.frame(specificity = c(youden$specificity,weighted$specificity,youden_years$specificity,weighted_years$specificity),
                    sensitivity = c(youden$sensitivity,weighted$sensitivity,youden_years$sensitivity,weighted_years$sensitivity),
                    threshold = c(youden$threshold,weighted$threshold,youden_years$threshold,weighted_years$threshold),
                    cost=c('equal','c/i = 10x','equal','c/i = 10x'),
                    model=c('s(dn) + s(ds)','s(dn) + s(ds)','s(dn) + s(years)','s(dn) + s(years)'))
p4=ggroc(mod3_roc) + 
  geom_point(data=opt_df |> filter(model=='s(dn) + s(ds)'),aes(x=specificity,y=sensitivity,color=cost)) +
  theme_bw() +
  annotate('text',x=0.2,y=0.9,label=paste('auc=',round(auc(df$circulating,df$cVDPVness),2),sep=''))

p5=ggplot() +
  stat_ecdf(data=df,aes(x=cVDPVness,color=circulating)) +
  geom_hline(data=opt_df|> filter(model=='s(dn) + s(ds)'),aes(yintercept = threshold, linetype=cost)) +
  theme_bw() +
  guides() +
  scale_color_manual(values=cmap)


p6 =ggplot(data=df,aes(x=circulating,y=log(cVDPVness/(1-cVDPVness)),fill=circulating,group=circulating)) +
  geom_violin(scale='width') +
  geom_quasirandom(varwidth = TRUE,width=0.45,alpha=0.3)+
  geom_hline(data=opt_df|> filter(model=='s(dn) + s(ds)'),
             aes(yintercept = log(threshold/(1-threshold)), linetype=cost)) +
  theme_bw() +
  scale_fill_manual(values=cmap) +
  guides() +
  ylab('logit(cVDPVness)') + xlab('')

p4+p5 +p6

p4.1=ggroc(mod5_roc) + 
  geom_point(data=opt_df |> filter(model!='s(dn) + s(ds)'),aes(x=specificity,y=sensitivity,color=cost)) +
  theme_bw() +
  annotate('text',x=0.2,y=0.9,label=paste('auc=',round(auc(df$circulating,df$cVDPVness_years),2),sep=''))

p5.1=ggplot() +
  stat_ecdf(data=df,aes(x=cVDPVness_years,color=circulating)) +
  geom_hline(data=opt_df|> filter(model!='s(dn) + s(ds)'),aes(yintercept = threshold, linetype=cost)) +
  theme_bw() +
  guides()+
  scale_color_manual(values=cmap)

p6.1 =ggplot(data=df,aes(x=circulating,y=log(cVDPVness_years/(1-cVDPVness_years)),fill=circulating,group=circulating)) +
  geom_violin(scale='width') +
  geom_quasirandom(varwidth = TRUE,width=0.45,alpha=0.3)+
  geom_hline(data=opt_df|> filter(model!='s(dn) + s(ds)'),
             aes(yintercept = log(threshold/(1-threshold)), linetype=cost)) +
  theme_bw() +
  guides() +
  ylab('logit(cVDPVness)') + xlab('')+
  scale_fill_manual(values=cmap)

p4.1+p5.1+p6.1

# level sets

# fit is your mgcv::gam with z ~ s(x) + s(y)
  
  # 1) grid over the domain
  d_s_range = seq(0,2600,by=1)
  d_n_range = seq(0,50,by=1)
  grid <- expand.grid(d_s = d_s_range, d_n = d_n_range)
  
  # 2) predict on the desired scale
  #   - use type="response" if C is on the response scale
  #   - if C is on the link scale, use type="link"
  zmat <- matrix(predict(mod3, newdata = grid, type = "response"),
                 nrow = length(d_s_range), ncol = length(d_n_range))
  
  # 3) extract level set(s) z = C
  cl <- contourLines(x = d_s_range, y = d_n_range, z = zmat, levels = opt_df$threshold[1:2])
  
  # cl is a list of components of the level set; each has $x and $y
  # Example: bind all polylines into one data.frame
  lev_df <- do.call(rbind, lapply(seq_along(cl), function(i)
    data.frame(comp = i, d_s = cl[[i]]$x, d_n = cl[[i]]$y)
  ))
  lev_df$cost = unique(opt_df$cost)[lev_df$comp]


p7=ggplot(df,aes(x=d_s,y=d_n)) +
  geom_jitter(aes(color=circulating,fill=circulating),width=0.3,height=0.3,alpha=0.2) +
  geom_line(data =lev_df,aes(linetype = cost)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_color_manual(values=cmap)
p7



# fit is your mgcv::gam with z ~ s(x) + s(y)

# 1) grid over the domain
years_range = seq(0,30,by=0.1)
d_n_range = seq(0,50,by=1)
grid <- expand.grid(years = years_range, d_n = d_n_range)

# 2) predict on the desired scale
#   - use type="response" if C is on the response scale
#   - if C is on the link scale, use type="link"
zmat <- matrix(predict(mod5, newdata = grid, type = "response"),
               nrow = length(years_range), ncol = length(d_n_range))

# 3) extract level set(s) z = C
cl <- contourLines(x = years_range, y = d_n_range, z = zmat, levels = opt_df$threshold[3:4])

# cl is a list of components of the level set; each has $x and $y
# Example: bind all polylines into one data.frame
lev_df_years <- do.call(rbind, lapply(seq_along(cl), function(i)
  data.frame(comp = i, years = cl[[i]]$x, d_n = cl[[i]]$y)
))
lev_df_years$cost = unique(opt_df$cost)[lev_df_years$comp]


p8=ggplot(df,aes(x=years,y=d_n)) +
  geom_jitter(aes(color=circulating,fill=circulating),width=0.3,height=0.3,alpha=0.2) +
  geom_line(data =lev_df_years,aes(linetype = cost)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('approx. years since OPV') +
  scale_x_continuous(breaks=seq(0,30,by=5)) +
  scale_color_manual(values=cmap)
p8


(p4 + p5 + p6) / (p7 + p8)
ggsave('classifier_dnds.png',units='in',width=12,height=6)
(p4.1 + p5.1 + p6.1)
ggsave('classifier_dnyears.png',units='in',width=12,height=3.5)

## classify vdpvs
meta_pred = meta 
meta_pred$cVDPVness = predict(mod3,newdata = meta_pred,type='response')

meta_pred |> group_by(circulating) |> summarize(mean(cVDPVness))



p9=ggplot(meta_pred,aes(x=d_s,y=d_n,color=circulating,fill=circulating)) +
  geom_jitter(width=0.3,height=0.3,alpha=0.2) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_color_manual(values=cmap)
p9

p10=ggplot(meta_pred,aes(x=d_s,y=d_n)) +
  geom_jitter(aes(color=cVDPVness,fill=cVDPVness),width=0.3,height=0.3,alpha=0.2) +
  geom_line(data =lev_df,aes(linetype = cost)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  guides(shape=guide_legend('serotype')) +
  facet_grid('.~circulating')
p10

p9 /p10 + plot_layout(guides='collect')
ggsave('cVDPVness.png',units='in',width=12,height=6)

