# Example orthogonalization code
# Author: Mary Lofton
# Date: 09OCT25

current_dat <- read_csv("./data/example_ecoregion_data.csv") # EDIT THIS TO BE PATH TO THE FILE ON YOUR COMPUTER

r_s <- current_dat %>%
  pull(Dep_Sdiff)

a_s <- current_dat %>%
  pull(Dep_Shistoric)

r_n <- current_dat %>%
  pull(Dep_Ndiff)

a_n <- current_dat %>%
  pull(Dep_Nhistoric)

n_s = length(r_s)
n_n = length(r_n)

S = cbind(1, r_s, a_s) # design matrix needs to include intercept [1, S]
P_s = diag(n_s) - (S %*% solve(t(S) %*% S, t(S))) # projection onto span{1,S}
rtilde_n = as.numeric(P_s %*% r_n)
atilde_n = as.numeric(P_s %*% a_n)

rtilde_n_centered = scale(rtilde_n, center = TRUE)
atilde_n_centered = scale(atilde_n, center = TRUE)

P_n = diag(n_n) - (rtilde_n_centered %*% solve(t(rtilde_n_centered) %*% rtilde_n_centered, t(rtilde_n_centered)))

astar_n = as.numeric(P_n %*% atilde_n_centered)

c0 = (t(rtilde_n_centered) %*% atilde_n_centered) / (t(rtilde_n_centered) %*% rtilde_n_centered)
c = c0[1,]

current_dat$Dep_Nhistoric_ortho = astar_n
current_dat$Dep_Ndiff_ortho = rtilde_n_centered
current_dat$c = c
