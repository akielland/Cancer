#

n1 <- df03 |> filter(timepoint == "SUR") |> 
  count(TrialArmNeo) |> 
  filter(TrialArmNeo == "Letro+Ribo")

n2 <- df03 |> filter(timepoint == "SUR") |> 
  count(TrialArmNeo) |> 
  filter(TrialArmNeo == "AC+pacli")

y1 <- df03 |> filter(timepoint == "SUR") |> 
  filter(TrialArmNeo == "Letro+Ribo") |> 
  filter(ROR.S.Group..Subtype.Only. == "low") |> 
  count(ROR.S.Group..Subtype.Only.)

y2 <- df03 |> filter(timepoint == "SUR") |> 
  filter(TrialArmNeo == "AC+pacli") |> 
  filter(ROR.S.Group..Subtype.Only. == "low") |> 
  count(ROR.S.Group..Subtype.Only.)

n <- c(48, 51)
n = c(as.numeric(n1[2]), as.numeric(n2[2]))
show(n)
y=c(39,37)
y = c(as.numeric(y1[2]), as.numeric(y2[2]))
show(y)

# uniform prior on theta
alpha=beta=1

# simulating posterior draws
S = 10000
p.sim = cbind(rbeta(S,alpha+y[1],beta+n[1]-y[1]),
              rbeta(S,alpha+y[2],beta+n[2]-y[2]))
hist(p.sim[,1]-p.sim[,2], breaks = 50)
show(mean(p.sim[,1]>p.sim[,2]))
show(quantile(p.sim[,1]-p.sim[,2],c(0.025,0.975)))


#b: 
p1 = seq(0,1,length=1000)
p2 = seq(0,1,length=1000)
p = expand.grid(p1,p2)
delta = diff(p1[1:2])*diff(p2[1:2])
pr.hat = delta*sum(dbeta(p[,1],alpha+y[1],beta+n[1]-y[1])*
                     dbeta(p[,2],alpha+y[2],beta+n[2]-y[2])*as.numeric(p[,1]>p[,2]))
show(pr.hat)
