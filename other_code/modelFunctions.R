# lognormal plot

lognormal <- function(Ndep, p4, p5){
  y = exp(-0.5 * (log(Ndep / p4)/ p5) ^ 2)
}

p4 = 2000
p5 = 3.8

curve(lognormal(x, p4 = p4, p5 = p5),from=0, to=4000000,ylab='growth',xlab = 'ndep',ylim = c(0,1))
abline(v = 3, col = "red")
abline(v = 35, col = "red")

growth <- function(tree_effect, agb, p2){
  y = tree_effect * agb ^ p2
}

tree_effect = -7
p2 = 0.475

curve(growth(x, tree_effect = tree_effect, p2 = p2),from=0, to=80000,ylab='growth',xlab = 'size')


comp <- function(ba_gt,p3){
  exp(-ba_gt*p3) 
}

p3 = 2.05
curve(comp(x, p3 = p3),from = 0, to = 10)

curve(dlnorm(x, meanlog = 2000, sdlog = 3.8), from= 0, to = 40000000)

sig <- function(ndep,p4){
  y = p4 * ndep
  }
p4 = 0.9
curve(sig(x, p4 = p4),from = 0, to = 10)

