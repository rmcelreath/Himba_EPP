# Himba extra-pair paternity estimation
# This script constructs a simple Bayesian compatibility interval for the rate of extra-pair paternity (EPP) in the sample.
# The model clusters by mother, to account for variation among women, and by mother-husband pairs, to account for variation among dyads.

library(rethinking)

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
dout <- d[,1:4]
write.csv( dout , file="Himba_id_columns.csv" , row.names=FALSE )

#############################
# raw EPP rate

mean( d$SF_Result[d$ParentsMd=="Yes"] )
# [1] 0.4891304

#############################
# construct interval by estimation

# when MotherID missing, assign unique value
mom_na <- which( is.na(d$MotherID) )
mom_id <- d$MotherID
mom_id[mom_na] <- mom_na + max(d$MotherID,na.rm=TRUE)
mom_id <- mom_id[d$ParentsMd=="Yes"]

# model
dat_list <- list(
    y = d$SF_Result[d$ParentsMd=="Yes"],
    mom_id = coerce_index( as.factor(mom_id) )
)

# model without any clustering
m1 <- ulam(
    alist(
        y ~ bernoulli( p ),
        p ~ beta(2,2) # 2,20 gives a prior skewed to small values
    ) , data=dat_list , chains=4 )

precis(m1)

# model clustering on mom
# need to use non-centered parameterization due to sparseness for many mothers
m2 <- ulam(
    alist(
        y ~ bernoulli( p ),
        logit(p) <- a + zm[mom_id]*sigma_mom,
        a ~ normal(0,1.5),
        zm[mom_id] ~ normal(0,1),
        sigma_mom ~ exponential(1)
    ) , data=dat_list , chains=4 )

precis(m2)
post <- extract.samples(m2)
p <- inv_logit(post$a)
precis(list(p=p),prob=0.95)

# cluster on mother-husband dyads
sf_id <- d$SocFatherID[d$ParentsMd=="Yes"]
sf_na <- which( is.na(sf_id) )
sf_id[sf_na] <- sf_na + max(sf_id,na.rm=TRUE)

dat_list$sf_id <- coerce_index( as.factor(sf_id) )

mf_dyad <- paste( dat_list$mom_id , dat_list$sf_id , sep="_" )
dat_list$dyad_id <- coerce_index( as.factor(mf_dyad) )

m3 <- ulam(
    alist(
        y ~ bernoulli( p ),
        logit(p) <- a + zm[mom_id]*sigma_mom + zd[dyad_id]*sigma_dyad,
        a ~ normal(0,1.5),
        zm[mom_id] ~ normal(0,1),
        sigma_mom ~ exponential(1),
        zd[dyad_id] ~ normal(0,1),
        sigma_dyad ~ exponential(2)
    ) , data=dat_list , chains=4 )

precis(m3,2)

post3 <- extract.samples(m3)
p3 <- inv_logit( post3$a )
precis( list(p=p3))

##############
# plotting

# plot post against prior
p_prior <- inv_logit( rnorm(1e4,0,1.5) )
post1 <- extract.samples(m1)

blank2(h=0.6)

dens( post1$p , lwd=1.5 , xlab="extra-pair paternity" , xlim=c(0,1) )

dens( p , lwd=1.5 , add=TRUE , col="red" )

dens( p_prior , add=TRUE , lty=2 )

text( 0.58 , 8 , "raw" , cex=0.8 )
text( 0.23 , 4 , "clustered\non mother" , cex=0.8 , col="red" )
text( 0.8 , 2 , "prior" , cex=0.8 )

dens( p3 , add=TRUE , col=rangi2 , lwd=1.5 )

# variance among moms
blank2(h=0.6)

dens( post$sigma_mom , xlim=c(0,5) , xlab="standard deviation among mothers" , lwd=1.5 )

# plot every mother's posterior
dens( p , lwd=1.5 , xlab="extra-pair paternity" , xlim=c(0,1) , col="red" , ylim=c(0,7) )
for ( i in dat_list$mom_id ) {
    pm <- inv_logit( post$a + post$zm[,i]*post$sigma_mom )
    dens( pm , add=TRUE , col=grau(0.1) )
}
mtext( "EPP by mother" )

