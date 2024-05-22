library(tidyverse)

# fetch runtime data
###################################################################
###################################################################
###################################################################
domain_name <- "parking"
p.ffpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/ffpo-parking.csv"
p.cgpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cgpo-parking.csv"
p.cg <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cg-parking.csv"
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
x = c(9,10,11,12,13,14,16,17,18,19,20,22,23,24,25,26,28,29,30,31,32,34,35,36,37,39,40,41,42,43)
dat.ffpo$instance[1:8] <- x[1:8]
dat.ffpo$instance[9] <- x[10]
for(i in 1:length(dat.ffpo$instance)){
  rho1 <- dat.ffpo$instance[i]
  dat.ffpo$instance[i] <- 1+ sqrt( (9-rho1)^2 + (16-(2*rho1-2))^2 )
}
dat.cgpo$instance <- dat.ffpo$instance[1:8]
dat.cg$instance <- dat.ffpo$instance[1:8]
rm(p.ffpo, p.cgpo, p.cg,i,x,rho1)


domain_name <- "gripper"
p.ffpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/ffpo-gripper.csv"
p.cgpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cgpo-gripper.csv"
p.cg <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cg-gripper.csv"
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
dat.ffpo$instance <- 15 + 5*dat.ffpo$instance
dat.cgpo$instance <- 15 + 5*dat.cgpo$instance
dat.cg$instance <- 15 + 5*dat.cg$instance
rm(p.ffpo, p.cgpo, p.cg)


domain_name <- "scanalyzer"
p.ffpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/ffpo-scanalyzer.csv"
p.cgpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cgpo-scanalyzer.csv"
p.cg <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cg-scanalyzer.csv"
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
params.scanalyzer <- c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.scanalyzer[i]}
for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.scanalyzer[i]}
for(i in dat.cg$instance){dat.cg$instance[i] <- params.scanalyzer[i]}
rm(p.ffpo, p.cgpo, p.cg, i)

dat.ffpo.tibble <- tibble(dat.ffpo) %>% mutate(prob_size=instance, planner="ffpo") %>% select(planner,prob_size,runtime)
dat.cgpo.tibble <- tibble(dat.cgpo) %>% mutate(prob_size=instance, planner="cgpo") %>% select(planner,prob_size,runtime)
dat.cg.tibble <- tibble(dat.cg) %>% mutate(prob_size=instance, planner="cg") %>% select(planner,prob_size,runtime)
dat <- bind_rows(dat.ffpo.tibble, dat.cgpo.tibble, dat.cg.tibble,)



# perform BEE
stats <- SimulStats(dat.ffpo=dat.ffpo, dat.cgpo=dat.cgpo, dat.cg=dat.cg,
                    num.iters=1000, plot=FALSE, clean.thresh=0.25)


# plot runtime data
ggplot(
  data = dat,
  mapping = aes(x=prob_size, y=runtime)
) + 
  geom_point(mapping=aes(color=planner, shape=planner), size=1.5) +
  labs(
    x = "Problem size",y = "Runtime"
  ) + 
  theme(legend.position = c(0.12,.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)
        )


# k densities
dk1 <- tibble(stats$k.ffpo) %>% 
  mutate(planner="ffpo", k_val=stats$k.ffpo) %>% 
  select(planner, k_val)
dk2 <- tibble(stats$k.cgpo) %>% 
  mutate(planner="cgpo", k_val=stats$k.cgpo) %>% 
  select(planner, k_val)
dk3 <- tibble(stats$k.cg) %>% 
  mutate(planner="cg", k_val=stats$k.cg) %>% 
  select(planner, k_val)
dk <- bind_rows(dk1, dk2, dk3)

ggplot(
  data = dk,
  mapping = aes(x=k_val, color = planner)
) + geom_density(size=1.25, bw=0.0005) + 
  labs(x="k value", y="density") +
  theme(legend.position = c(0.12,.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)
  )



# m densities
dm1 <- tibble(stats$m.ffpo) %>% 
  mutate(planner="ffpo", m_val=stats$m.ffpo) %>% 
  select(planner, m_val)
dm2 <- tibble(stats$m.cgpo) %>% 
  mutate(planner="cgpo", m_val=stats$m.cgpo) %>% 
  select(planner, m_val)
dm3 <- tibble(stats$m.cg) %>% 
  mutate(planner="cg", m_val=stats$m.cg) %>% 
  select(planner, m_val)
dm <- bind_rows(dm1, dm2, dm3)

ggplot(
  data = dm,
  mapping = aes(x=m_val, color = planner)
) + labs(
  x = "m value",
  y = "density"
  ) + geom_density(size=1.25) + 
  theme(legend.position = c(0.88,.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)
  )




# pairwise compare parameters

# k
ParamCompare(stats$k.ffpo, stats$k.ffpo); ParamCompare(stats$k.ffpo, stats$k.cgpo); ParamCompare(stats$k.ffpo, stats$k.cg); 
ParamCompare(stats$k.cgpo, stats$k.ffpo); ParamCompare(stats$k.cgpo, stats$k.cgpo); ParamCompare(stats$k.cgpo, stats$k.cg); 
ParamCompare(stats$k.cg, stats$k.ffpo); ParamCompare(stats$k.cg, stats$k.cgpo); ParamCompare(stats$k.cg, stats$k.cg); 

# m
ParamCompare(stats$m.ffpo, stats$m.ffpo); ParamCompare(stats$m.ffpo, stats$m.cgpo); ParamCompare(stats$m.ffpo, stats$m.cg); 
ParamCompare(stats$m.cgpo, stats$m.ffpo); ParamCompare(stats$m.cgpo, stats$m.cgpo); ParamCompare(stats$m.cgpo, stats$m.cg); 
ParamCompare(stats$m.cg, stats$m.ffpo); ParamCompare(stats$m.cg, stats$m.cgpo); ParamCompare(stats$m.cg, stats$m.cg); 








