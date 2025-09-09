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

#indel handling?
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
  aln_tmp = gap_as_state(as.phyDat(aln))
  
  # translate AA
  aln_aa = as.AAbin(translate(DNAStringSet(apply(as.character(aln), 1, paste0, collapse = "")), 
                              if.fuzzy.codon = "X"))
  aln_aa_tmp = gap_as_state(as.phyDat(aln_aa))
  
  # counts different from Sabin
  for (CLASS in c('cVDPV','iVDPV','wild', 'SabinLike','VDPV')){
    keep = df |> filter(classification == CLASS | classification=='Sabin') |> select(accession)
    
    nt_d <- as.matrix(dist.ml(aln_tmp[labels(aln) %in% keep$accession,] , model = "F81"))[1,] * gnl
    # aa_d <- as.matrix(dist.ml(aln_aa_tmp[labels(aln_aa) %in% keep$accession] , model = "FLU"))[1,]*gnl/3
    aa_d <- as.matrix(dist.hamming(aln_aa_tmp[labels(aln_aa) %in% keep$accession]))[1,]*gnl/3
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

cmap = c("#04e762","#ae7766","#00a1e4","#dc0073","#f5b700")
names(cmap)=levels(meta$classification)


meta |> 
  filter(classification == 'cVDPV') |>
  filter(years >=1) |> 
  summarize(mean_d_n = mean(d_n),
            median_d_n = median(d_n),
            harmonic_d_n = 1/(mean(1/d_n)))



meta |> 
  filter(classification == 'cVDPV') |>
  filter(years >=1) |>
  mutate(d_n = round(d_n)) |>
  group_by(d_n) |> summarize(n=n()) |> arrange(desc(n))




meta = meta |> mutate(d_n_from_consensus = abs(d_n-3.7)) |>
  mutate(relative_fitness = 0.95^(2200/900*d_n_from_consensus)) |>
  mutate(relative_fitness = relative_fitness/mean(relative_fitness[classification=='cVDPV' & years >=1]))

summary(lm(d_n_from_consensus ~ years,
           data = meta |> filter(years >=1) |> 
             filter(classification == 'cVDPV'))) 

sabin_df = data.frame(serotype = c('1','2','3'),
                      relative_fitness = 2.3/c(14,8,18))
ggplot(meta |> 
         # dplyr::filter(classification %in% c('cVDPV','iVDPV')) |> ungroup() |> 
         droplevels(),
       aes(x=years,y=relative_fitness,color=classification,fill=classification)) +
  geom_point(aes(shape=serotype),alpha=0.2) +
  geom_smooth(method='gam',fill='black',alpha=0.3) +
  geom_hline(data=sabin_df,aes(yintercept = relative_fitness,linetype=serotype)) +
  ylab('relative fitness') + 
  xlab('approx. years since OPV') +
  ylim(c(0,1.2)) +
  # scale_x_continuous(breaks=seq(0,30,by=5)) +
  scale_color_manual(values=cmap)

# consider consensus for wild
tmp=meta |>
  filter(!(grepl('finland',paper,ignore.case = TRUE) | grepl('2010',paper,ignore.case = TRUE))) |>
  group_by(classification,serotype) |>
  summarize(mean_d_n = mean(d_n[years>=0]),
            mean_d_s = mean(d_s[years>=0]),n=n())
tmp
meta = meta |> 
  left_join(tmp) |>
  mutate(adjusted_d_n = if_else(classification == 'wild',d_n-(mean_d_n-3.8),d_n),
                      adjusted_d_s = if_else(classification == 'wild',d_s-(mean_d_s-25.5),d_s))
  
summary(wpv_mod <- glm(d_n ~ d_s + classification, family = 'gaussian',
                       data = meta |> filter(classification %in% c('cVDPV','wild')) |> mutate(d_n = round(d_n), d_s=round(d_s))))
summary(wpv_mod <- glm(adjusted_d_n ~ adjusted_d_s + classification, family = 'gaussian',
                       data = meta |> filter(classification %in% c('cVDPV','wild')) |> mutate(d_n = round(d_n), d_s=round(d_s))))


meta = meta |> mutate(adjusted_d_n_from_consensus = abs(adjusted_d_n-3.7)) |>
  mutate(relative_fitness = 0.95^(2200/900*adjusted_d_n_from_consensus)) |>
  mutate(relative_fitness = relative_fitness/mean(relative_fitness[classification=='cVDPV' & years >=1]))

ggplot(meta |> 
         # dplyr::filter(classification %in% c('cVDPV','iVDPV')) |> ungroup() |> 
         droplevels(),
       aes(x=years,y=relative_fitness,color=classification,fill=classification)) +
  geom_point(aes(shape=serotype),alpha=0.2) +
  geom_smooth(method='gam',fill='black',alpha=0.3) +
  geom_hline(data=sabin_df,aes(yintercept = relative_fitness,linetype=serotype)) +
  ylab('relative fitness') + 
  xlab('approx. years since OPV') +
  ylim(c(0,1.2)) +
  # scale_x_continuous(breaks=seq(0,30,by=5)) +
  scale_color_manual(values=cmap)



# auspice output
rel_fit = round(meta$relative_fitness,3)
meta$relative_fitness__colour <- scales::col_numeric(viridisLite::viridis(length(rel_fit)), domain = range(rel_fit, na.rm = TRUE), na.color = NA)(rel_fit)
write.csv(meta,'processed_metadata.csv',row.names = FALSE)

## plots

ggplot(meta,aes(x=years,fill=classification)) +
  geom_histogram(binwidth = 0.25,boundary=0) +
  facet_grid(. ~ classification, scales='free') +
  xlab('approx. years since OPV') +
  # scale_x_continuous(breaks=seq(0,30,by=5)) +
  scale_fill_manual(values=cmap)

ggplot(meta,aes(x=d_n,fill=classification)) +
  geom_histogram(binwidth = 1,boundary=-0.5) +
  facet_grid(. ~ classification) +
  xlab('non-synonymous substitutions') +
  scale_fill_manual(values=cmap)

ggplot(meta,aes(x=d_s,fill=classification)) +
  geom_histogram(binwidth = 1,boundary=-0.5) +
  facet_grid(. ~ classification, scales='free') +
  # xlim(c(0,20)) +
  xlab('synonymous substitutions') +
  scale_fill_manual(values=cmap)


p1=ggplot(meta,aes(x=d_s,y=d_n,color=classification,shape=serotype)) +
  geom_point(alpha=0.2) +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_color_manual(values=cmap)
p1

p1.5=ggplot(meta,aes(x=adjusted_d_s,y=adjusted_d_n,color=classification,shape=serotype)) +
  geom_point(alpha=0.2) +
  ylab('adjusted non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_color_manual(values=cmap)
p1.5


df = meta |> dplyr::filter(classification %in% c('cVDPV','iVDPV')) |> ungroup() |> droplevels()
p2=ggplot(df,
          aes(x=years,y=d_n,color=classification,fill=classification,group = classification)) +
  geom_point(alpha=0.2) +
  geom_smooth(method='gam',fill='black',alpha=0.3) +
  ylab('non-synonymous substitutions') +
  xlab('approx. years since OPV') +
  scale_x_continuous(breaks=seq(0,30,by=5)) +
  scale_color_manual(values=cmap)
p2

p3=ggplot(df,aes(x=d_s,y=d_n,color=classification,fill=classification)) +
  geom_point(aes(shape=serotype),alpha=0.2) +
  geom_smooth(method='gam',fill='black',alpha=0.3) +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_color_manual(values=cmap)
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
  geom_point(data=df,   aes(x=100*(1-nt_diff/903), y=1-d_n/301, color=cmap[3])) +
  geom_point(data=zhao, aes(x=100*(1-nt_diff/903), y=1-d_n/301, color=cmap[4])) +
  facet_grid(. ~ classification) +
  xlab("nucleotide identity (%)") +
  ylab("amino acid identity (%)") +
  scale_color_identity(
    name   = "source",
    breaks = c(cmap[3], cmap[4]),
    labels = c("Distance from Sabin", "Zhao 2017\ninter-isolate distance"),
    guide  = "legend"
  )
ggsave('supp_Zhoa2017_consistency.png',units='in',width=8, height=3.5)

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

mod1 = gam(formula = classification ~ adjusted_d_n + adjusted_d_s,
           data =df,
           family='binomial')
mod2 = gam(formula = classification ~ adjusted_d_n + adjusted_d_s + serotype,
           data =df,
           family='binomial')
summary(mod2)
# mod1 is good enough since no sig diff between types 1 and 2 (and no data for 3)
summary(mod1)
plot(mod1)

# nonlinear model is better than linear model
mod3 = gam(formula = classification ~ s(adjusted_d_n) + s(adjusted_d_s),
           data =df,
           family='binomial')
summary(mod3)
c(AIC(mod1),AIC(mod3))
plot(mod3)


# interaction model looks pretty overfit
mod4 = gam(formula = classification ~ ti(adjusted_d_n,adjusted_d_s),
           data =df,
           family='binomial')
summary(mod4)
c(AIC(mod3),AIC(mod4))
plot(mod4)

# years instead of d_s is way worse
mod5 = gam(formula = classification ~ s(adjusted_d_n) + s(years),
           data =df,
           family='binomial')
summary(mod5)
c(AIC(mod3),AIC(mod5))
plot(mod5)

# mod3=mod1 # wow mod1 is way worse, as one can see looking at the classifier plot p6 below
df$cVDPVness = predict(mod3,newdata = df,type='response')
df$cVDPVness_years = predict(mod5,newdata = df,type='response')

meta$cVDPVness = predict(mod3,newdata = meta,type='response')
meta$cVDPVness_years = predict(mod5,newdata = meta,type='response')

# ROC analysis
mod3_roc = roc(df$classification,df$cVDPVness)
mod5_roc = roc(df$classification,df$cVDPVness_years)

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
  annotate('text',x=0.2,y=0.9,label=paste('auc=',round(auc(df$classification,df$cVDPVness),2),sep=''))

p5=ggplot() +
  stat_ecdf(data=df,aes(x=cVDPVness,color=classification)) +
  geom_hline(data=opt_df|> filter(model=='s(dn) + s(ds)'),aes(yintercept = threshold, linetype=cost)) +
  theme_bw() +
  guides() +
  scale_color_manual(values=cmap)


p6 =ggplot(data=df,aes(x=classification,y=log(cVDPVness/(1-cVDPVness)),fill=classification,group=classification)) +
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
  annotate('text',x=0.2,y=0.9,label=paste('auc=',round(auc(df$classification,df$cVDPVness_years),2),sep=''))

p5.1=ggplot() +
  stat_ecdf(data=df,aes(x=cVDPVness_years,color=classification)) +
  geom_hline(data=opt_df|> filter(model!='s(dn) + s(ds)'),aes(yintercept = threshold, linetype=cost)) +
  theme_bw() +
  guides()+
  scale_color_manual(values=cmap)

p6.1 =ggplot(data=df,aes(x=classification,y=log(cVDPVness_years/(1-cVDPVness_years)),fill=classification,group=classification)) +
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
  d_s_range = seq(0,140,by=1)
  d_n_range = seq(0,50,by=1)
  grid <- expand.grid(adjusted_d_s = d_s_range, adjusted_d_n = d_n_range)
  
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
    data.frame(comp = i, adjusted_d_s = cl[[i]]$x, adjusted_d_n = cl[[i]]$y)
  ))
  lev_df$cost = unique(opt_df$cost)[lev_df$comp]


p7=ggplot(df,aes(x=adjusted_d_s,y=adjusted_d_n)) +
  geom_jitter(aes(color=classification,fill=classification),width=0.3,height=0.3,alpha=0.2) +
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
grid <- expand.grid(years = years_range, adjusted_d_n = d_n_range)

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
  data.frame(comp = i, years = cl[[i]]$x, adjusted_d_n = cl[[i]]$y)
))
lev_df_years$cost = unique(opt_df$cost)[lev_df_years$comp]


p8=ggplot(df,aes(x=years,y=adjusted_d_n)) +
  geom_jitter(aes(color=classification,fill=classification),width=0.3,height=0.3,alpha=0.2) +
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
meta |> group_by(classification) |> summarize(mean(cVDPVness))



p9=ggplot(meta,aes(x=d_s,y=d_n,color=classification,fill=classification)) +
  geom_jitter(width=0.3,height=0.3,alpha=0.2) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_color_manual(values=cmap)
p9

p9.5=ggplot(meta,aes(x=adjusted_d_s,y=adjusted_d_n,color=classification,fill=classification)) +
  geom_jitter(width=0.3,height=0.3,alpha=0.2) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_color_manual(values=cmap)
p9.5

p10=ggplot(meta,aes(x=adjusted_d_s,y=adjusted_d_n)) +
  geom_point(aes(color=cVDPVness,fill=cVDPVness),alpha=0.2) +
  geom_point(data=meta |> filter(accession %in% c('AF233112','JN812668')), color='red') +
  geom_point(data=meta |> filter(accession %in% c('AF233098')), color='green') + 
  geom_point(data=meta |> filter(accession %in% c('FJ842177')), color='magenta') + 
  geom_line(data =lev_df,aes(linetype = cost)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('non-synonymous substitutions') +
  xlab('synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  guides(shape=guide_legend('serotype')) +
  facet_grid('.~classification')
p10

p9 /p10 + plot_layout(guides='collect')
ggsave('cVDPVness.png',units='in',width=12,height=6)



# dig into cVDPVness scores and sources

meta |> filter(classification=='cVDPV') |>
  group_by(paper) |> 
  summarize(mean_cVDPVness = mean(cVDPVness)) |> 
  arrange(desc(mean_cVDPVness)) |> 
  View()

meta |> filter(classification=='cVDPV') |>
  ggplot() +
  geom_point(aes(x=years,y=cVDPVness,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))

meta |> filter(classification=='cVDPV') |>
  ggplot() +
  geom_point(aes(x=years,y=d_n,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))


meta |> filter(classification=='iVDPV') |>
  group_by(paper) |> 
  summarize(mean_cVDPVness = mean(cVDPVness)) |> 
  arrange(desc(mean_cVDPVness)) |> 
  View()

meta |> filter(classification=='iVDPV') |>
  ggplot() +
  geom_point(aes(x=years,y=cVDPVness,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))

meta |> filter(classification=='iVDPV') |>
  ggplot() +
  geom_point(aes(x=years,y=d_n,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))


meta |> filter(classification=='VDPV') |>
  group_by(paper) |> 
  summarize(mean_cVDPVness = mean(cVDPVness)) |> 
  arrange(desc(mean_cVDPVness)) |> 
  View()

meta |> filter(classification=='VDPV') |>
  ggplot() +
  geom_point(aes(x=years,y=cVDPVness,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))

meta |> filter(classification=='cVDPV') |>
  ggplot() +
  geom_point(aes(x=years,y=d_n,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))


meta |> filter(classification=='wild') |>
  group_by(paper) |> 
  summarize(mean_cVDPVness = mean(cVDPVness)) |> 
  arrange(desc(mean_cVDPVness)) |> 
  View()

meta |> filter(classification=='wild') |>
  ggplot() +
  geom_point(aes(x=years,y=cVDPVness,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))

meta |> filter(classification=='wild') |>
  ggplot() +
  geom_point(aes(x=years,y=d_n,color=cVDPVness)) +
  facet_wrap('paper',labeller = label_wrap_gen(width = 40))

meta |> filter(classification=='wild') |>
  arrange(cVDPVness) |> select(cVDPVness, accession, paper) |> View()

names(meta)
p11=ggplot(meta,aes(x=adjusted_d_s,y=relative_fitness)) +
  geom_point(aes(color=cVDPVness,fill=cVDPVness),alpha=0.2) +
  geom_point(data=meta |> filter(accession %in% c('AF233112','JN812668')), color='red') + #MOL91-PV1 and BRA85-PV3
  geom_point(data=meta |> filter(accession %in% c('AF233098')), color='green') + # AZ59-PV1
  geom_point(data=meta |> filter(accession %in% c('FJ842177') | grepl('finland',paper,ignore.case = TRUE)), color='magenta') + #Finland
  geom_point(data=meta |> filter(grepl('the Republic of Congo in 2010',paper,ignore.case = TRUE)), color='orange') + #ROC2010
  geom_hline(data=sabin_df,aes(yintercept = relative_fitness,linetype=serotype)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('relative fitness') +
  xlab('adjusted synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_y_continuous(breaks=seq(0,2,by=0.2)) +
  guides(shape=guide_legend('serotype')) +
  facet_grid('.~classification', scale='free_x')
p11

p12=ggplot(meta,aes(x=adjusted_d_s,y=relative_fitness)) +
  geom_point(aes(color=serotype,fill=serotype),alpha=0.2) +
  geom_point(data=meta |> filter(accession %in% c('AF233112','JN812668')), color='red') + #MOL91-PV1 and BRA85-PV3
  geom_point(data=meta |> filter(accession %in% c('AF233098')), color='green') + # AZ59-PV1
  geom_point(data=meta |> filter(accession %in% c('FJ842177') | grepl('finland',paper,ignore.case = TRUE)), color='magenta') + #Finland
  geom_point(data=meta |> filter(grepl('the Republic of Congo in 2010',paper,ignore.case = TRUE)), color='orange') + #ROC2010
  geom_hline(data=sabin_df,aes(yintercept = relative_fitness,linetype=serotype)) +
  # geom_smooth(method='gam') +
  theme_bw() +
  ylab('relative fitness') +
  xlab('adjusted synonymous substitutions') +
  scale_x_continuous(breaks=seq(0,300,by=20)) +
  scale_y_continuous(breaks=seq(0,2,by=0.2)) +
  guides(shape=guide_legend('serotype')) +
  facet_grid('.~classification', scale='free_x')
p12

meta |> 
  # filter( (classification == 'wild' & !grepl('finland',paper, ignore.case = TRUE))) |>
  filter((classification == 'wild' & !grepl('finland',paper, ignore.case = TRUE)) | (classification == 'cVDPV' & years>=1)) |>
  filter(relative_fitness>0.5) |>
  group_by(serotype) |> 
  summarize(median_fitness = median(relative_fitness),
            n=n())

