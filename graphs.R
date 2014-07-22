## Random graphs. 
# Erdos-Renyi
gnp = function(n=30, p=0.05)  erdos.renyi.game(n=n, p.or.m=p)
sworld <- function(dim=1, size=10, nei=2, p=0.15) 
    simplify(watts.strogatz.game(dim=dim, size=size, nei=nei, p=p))

barabasi = function(n=20, power=0.8) barabasi.game(n, power=power,directed=F)