# Himba extra-pair paternity estimation
# This script constructs a simple Bayesian compatibility interval for the rate of extra-pair paternity (EPP) in the sample.
# The model clusters by mother, to account for variation among women, and by mother-husband pairs, to account for variation among dyads.

library(rethinking) # need version 1.88 or higher

#############################
# load data

d <- read.csv( "HimbaPaternityResults_Nov2018.csv" )

# [1] "ParticipantID"              "SocFatherID"               
# [3] "MotherID"                   "BoyfriendID"               
# [5] "SF_Result"                  "BF_Result"                 
# [7] "Mother_Result"              "Actual.Father.ID"          
# [9] "ParentsMd"                  "Biofather.assertion.origin"
#[11] "Omoka.Dad"                  "Omoka.Mom"                 
#[13] "Self.Omoka.Status"   

# export version with only ID columns (first 4 cols)
# dout <- d[,1:4]
# write.csv( dout , file="Himba_id_columns.csv" , row.names=FALSE )

# read back in with corrections
din <- read.csv( "Himba_id_columns-final.csv" )
# replace in original
d[,1:4] <- din[,1:4]

#############################
# raw EPP rate

# include_idx <- which( d$ParentsMd=="Yes" & !is.na(d$SocFatherID) )
include_idx <- which( !is.na(d$SocFatherID) )

mean( 1 - d$SF_Result[ include_idx ] )
# [1] 0.4858757

#############################
# construct interval by estimation

# when MotherID missing, assign unique value
# mom_na <- which( is.na(d$MotherID) )
mom_id <- d$MotherID
# mom_id[mom_na] <- mom_na + max(d$MotherID,na.rm=TRUE)
mom_id <- mom_id[ include_idx ]

# model
dat_list <- list(
    y = 1L - d$SF_Result[ include_idx ],
    mom_id = coerce_index( as.factor( as.character(mom_id) ) )
)
check_index(dat_list$mom_id)

# model without any clustering
m1 <- ulam(
    alist(
        y ~ bernoulli( p ),
        p ~ beta(2,2) # 2,20 gives a prior skewed to small values
    ) , data=dat_list , chains=4 , log_lik=TRUE )

precis(m1,prob=0.95)

# incorporating 5% false positive paternity assignment
m1f <- ulam(
    alist(
        y|y==1 ~ custom( log_mix( 
            0.95 ,
            bernoulli_lpmf( 1 | p ),
            bernoulli_lpmf( 0 | p )
        )),
        y|y==0 ~ bernoulli( p ),
        p ~ beta(2,2) # 2,20 gives a prior skewed to small values
    ) , data=dat_list , chains=4 , sample=TRUE )

# model clustering on mom
# need to use non-centered parameterization due to sparseness for many mothers
m2 <- ulam(
    alist(
        y ~ bernoulli( p ),
        logit(p) <- a + z[mom_id]*sigma,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1)
    ) , data=dat_list , chains=4 , log_lik=TRUE , 
    constraints=list(sigma="lower=0") )

precis(m2)
post <- extract.samples(m2)
p <- inv_logit(post$a)
precis(list(p=p),prob=0.95)

# with false positive rate
m2f <- ulam(
    alist(
        y|y==1 ~ custom( log_mix( 
            0.95 ,
            bernoulli_lpmf( 1 | p ),
            bernoulli_lpmf( 0 | p )
        )),
        y|y==0 ~ bernoulli( p ),
        logit(p) <- a + z[mom_id]*sigma,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1)
    ) , data=dat_list , chains=4 , log_lik=FALSE , 
    constraints=list(sigma="lower=0") )

precis(m2f)
post <- extract.samples(m2f)
p <- inv_logit(post$a)
precis(list(p=p),prob=0.95)

# cluster on mother-husband dyads
# this is model to use

sf_id <- d$SocFatherID[ include_idx ]
sf_na <- which( is.na(sf_id) )
sf_id[sf_na] <- sf_na + max(sf_id,na.rm=TRUE)

dat_list$sf_id <- coerce_index( as.factor(sf_id) )

mf_dyad <- paste( dat_list$mom_id , dat_list$sf_id , sep="_" )
dat_list$dyad_id <- coerce_index( as.factor(mf_dyad) )
check_index(dat_list$dyad_id)

m3 <- ulam(
    alist(
        y ~ bernoulli( p ),
        logit(p) <- a + z[mom_id]*sigma + x[dyad_id]*tau,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat_list , chains=4 , iter=4000 , log_lik=TRUE ,
    constraints=list(sigma="lower=0",tau="lower=0") )

precis(m3)
dashboard(m3)

post3 <- extract.samples(m3)
p3 <- inv_logit( post3$a )
precis( list(p=p3) , prob=0.95 )

# now with false positive rate
m3f <- ulam(
    alist(
        y|y==1 ~ custom( log_mix( 
            0.95 ,
            bernoulli_lpmf( 1 | p ),
            bernoulli_lpmf( 0 | p )
        )),
        y|y==0 ~ bernoulli( p ),
        logit(p) <- a + z[mom_id]*sigma + x[dyad_id]*tau,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat_list , chains=4 , iter=4000 , log_lik=FALSE ,
    constraints=list(sigma="lower=0",tau="lower=0") )

precis(m3f)
dashboard(m3f)

post3f <- extract.samples(m3f)
p3f <- inv_logit( post3f$a )
precis( list(p=p3f) , prob=0.95 )

# dyads only
m4 <- ulam(
    alist(
        y ~ bernoulli( p ),
        logit(p) <- a + x[dyad_id]*tau,
        a ~ normal(0,1.5),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat_list , chains=4 , log_lik=TRUE ,
    constraints=list(tau="lower=0") )

precis(m4)

post4 <- extract.samples(m4)
p4 <- inv_logit( post4$a )
precis( list(p=p4))

##############
# comparison
compare( m1 , m2 , m3 , m4 )

#    WAIC pWAIC dWAIC weight   SE  dSE
#m3 196.5  31.4   0.0   0.68 9.96   NA
#m4 198.0  31.1   1.5   0.32 9.71 0.76
#m2 220.2  27.5  23.7   0.00 8.14 4.75
#m1 247.0   0.9  50.6   0.00 0.56 9.95

##############
# plotting

# plot post against prior
p_prior <- inv_logit( rnorm(1e4,0,1.5) )
post1 <- extract.samples(m1)
post1f <- extract.samples(m1f)

blank2(h=0.6)

dens( post1$p , lwd=1.5 , xlab="extra-pair paternity" , xlim=c(0,1) , col=col.alpha("black",0.4) )
dens( post1f$p , lwd=1.5 , add=TRUE )

#dens( p , lwd=1.5 , add=TRUE , col="red" )

dens( p_prior , add=TRUE , lty=2 )

text( 0.58 , 8 , "raw" , cex=0.8 )
text( 0.72 , 4 , "clustered" , cex=0.8 , col="red" )
text( 0.8 , 2 , "prior" , cex=0.8 )

dens( p3f , add=TRUE , col="red" , lwd=1.5 )
dens( p3 , add=TRUE , col=col.alpha("red",0.4) , lwd=1.5 )
#dens( p4 , add=TRUE , col="blue" , lwd=1.5 )



# variance among moms
blank2(h=0.6)

dens( post3$sigma , xlim=c(0,5) , xlab="standard deviation" , lwd=1.5 )
dens( post3$tau , add=TRUE , col="red ")

# plot every mother's posterior
dens( p3 , lwd=1.5 , xlab="extra-pair paternity" , xlim=c(0,1) , col="red" , ylim=c(0,7) )
for ( i in dat_list$mom_id ) {
    pm <- inv_logit( post3$a + post3$z[,i]*post3$sigma )
    dens( pm , add=TRUE , col=grau(0.1) )
}
mtext( "EPP by mother" )

# plot every dyad's posterior
dens( p3 , lwd=1.5 , xlab="extra-pair paternity" , xlim=c(0,1) , col="red" , ylim=c(0,7) )
for ( i in dat_list$dyad_id ) {
    pm <- inv_logit( post3$a + post3$x[,i]*post3$tau )
    dens( pm , add=TRUE , col=grau(0.1) )
}
mtext( "EPP by mother-husband dyad" )

# trace plots

traceplot( m3f , pars=c("a","sigma","tau","z[1]","z[2]","z[3]","x[1]","x[2]","x[3]") )

trankplot( m3f , pars=c("a","sigma","tau","z[1]","z[2]","z[3]","x[1]","x[2]","x[3]") )
