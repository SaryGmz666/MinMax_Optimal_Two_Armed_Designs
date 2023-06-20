# Comparacion de Frentes de Pareto obtenidos por los algoritmos PF-CE and EPCEA.

# Cargar los frentes de pareto del EPCEA.
n <- 100
p <- 6
G_SolND <- FALSE   #Si es TRUE cada iteracion se guardan las SND, FALSE para guardar la mejor solucion

load.file.pxe <- paste("Frentes_Pareto_EPCEA/nObs_", n, "_nCovar_", p, ".RData", sep = '')
load(load.file.pxe)

ii <- 5
load.file.pxe <- paste("PXE_nObs_", n, "_nCovar_", p, "_Inst_", ii, ".RData", sep = '')
load(load.file.pxe)
X.Sol.pxe <- X.Sol

load.file.pxe <- paste("AEC_nObs_", n, "_nCovar_", p, "_Inst_", ii,"_Opt_", 1,  ".RData", sep = '')
load(load.file.pxe)
X.Sol.G <- X.Sol
indicador <- find_nondominated_points(t(X.Sol.G[c(n + 1,n + 2),]))
X.Sol.G <- X.Sol.G[,indicador]

load.file.pxe <- paste("AEC_nObs_", n, "_nCovar_", p, "_Inst_", ii,"_Opt_", 2,  ".RData", sep = '')
load(load.file.pxe)
X.Sol.I <- X.Sol
indicador <- find_nondominated_points(t(X.Sol.I[c(n + 1,n + 2),]))
X.Sol.I <- X.Sol.I[,indicador]

frente.EPCEA <- cbind('Alg' = 1, 'I' = pareto.fronts[[ii]][,1],
                       'G' = pareto.fronts[[ii]][,2])
frente.PFCE <- cbind('Alg' = 2, 'I' = X.Sol.pxe[n + 2,],
                     'G' = X.Sol.pxe[n + 1,])
frente.AEC_G <- cbind('Alg' = 3, 'I' = X.Sol.G[n + 2,],
                     'G' = X.Sol.G[n + 1,])
frente.AEC_I <- cbind('Alg' = 4, 'I' = X.Sol.I[n + 2,],
                      'G' = X.Sol.I[n + 1,])

frentes <- rbind(frente.EPCEA, frente.PFCE, frente.AEC_G, frente.AEC_I)
G.lim <- c(min(frentes[,'G']),max(frentes[,'G']))
I.lim <- c(min(frentes[,'I']),max(frentes[,'I']))
plot(x=1, ylab = 'G-optimalidad', 
     xlab = 'I-optimalidad', pch = 3, type = 'n',
     xlim = I.lim, ylim = G.lim, 
     main = paste("Frentes de Pareto Instancia", ii))
points(frente.EPCEA[,2], frente.EPCEA[,3], pch = 3, col = 'red')
points(frente.PFCE[,2], frente.PFCE[,3], pch = 1, col = 'blue')
points(frente.AEC_G[,2], frente.AEC_G[,3], pch = 2, col = 'brown')
points(frente.AEC_I[,2], frente.AEC_I[,3], pch = 4, col = 'green')


# Crear grafica de soluciones no-dominadas
aux.vec <- find_nondominated_points(frentes[,c(2,3)])
non.dom <- frentes[aux.vec, ]

G.lim <- c(min(frentes[,'G']),max(frentes[,'G']))
I.lim <- c(min(frentes[,'I']),max(frentes[,'I']))
plot(x=1, ylab = 'G-optimalidad', 
     xlab = 'I-optimalidad', pch = 3, type = 'n',
     xlim = I.lim, ylim = G.lim, 
     main = paste("Soluciones No Dominadas Instancia", ii))
points(non.dom[non.dom[,1]==1,2], non.dom[non.dom[,1]==1,3], pch = 3, col = 'red')
points(non.dom[non.dom[,1]==2,2], non.dom[non.dom[,1]==2,3], pch = 1, col = 'blue')

