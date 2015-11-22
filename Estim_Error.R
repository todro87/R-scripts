#################################################################################
#################################################################################
# Ponizej sa funkcje do okreslania niepewnosci

# Funkcja zaleznosci j0' od j0 i n
fj <- function(res, barriers) {
        
        set.seed(100)
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        temp <- 300
        
        phi <- barriers[2]
        
        noise <- k*temp/q
        
        j0 <- res[3,1]
        
        n <- res[3,2]
        
        nn <- rnorm(1000, n, 0.1/3)
        
        alpha <- (nn-n)/n
        
        cbind(nn,j0*exp(alpha*q*phi/(nn*k*temp)))
        
}

#Funkcja wybiera dowolne parametry i wprowadza do nich szum, zwraca zaszumione
#paremetry i wspolczynnik jakosci dopasowania
noise <- function(pars, sel, data, temp=300, sd=0.1/3, nit=1000) {
        
        set.seed(100)
        
        
        j0 <- pars[1]
        n <- pars[2]
        rs <- pars[3]
        rsh <- pars[4]
        
        
        j0s <- rnorm(nit, j0, sd*j0)
        ns <- rnorm(nit, n, sd*n)
        rss <- rnorm(nit, rs, sd*rs)
        rshs <- rnorm(nit, rsh, sd*rsh)
        
        
        npars <- cbind(j0s, ns, rss, rshs)
        
        jv <- jvReal(pars, data=data, temp=temp)
        jsum <- sum((data[,2]-jv[,2])^2)
        
        tmppars <- pars
        
        store <- data.frame()
        
        for (i in seq(nit)) {
                
                tmppars[sel] <- npars[i,sel]
                
                print(c(i, tmppars))
                
                tmpjv <- jvReal(tmppars, data)
                
                tmpsum <- sum((data[,2]-tmpjv[,2])^2)
                
                store <- rbind(store, c(tmppars[sel], tmpsum))
                
        }
        
        names(store) <- c(seq(length(sel)+1))
        
        store
        
}

# Wyznacza j00 i phi z zaszumionych danych j i n
findBarrier <- function(data, res, temp) {
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        n <- data[,2]
        j <- data[,1]
        
        recipn <- 1/n/k/temp*q
        lnj <- log(j)
        
        #plot(cbind(recipn, lnj))
        
        res <- lm(lnj~recipn)
        
        j00 <- res$coefficients[[1]]
        
        phi <- res$coefficients[[2]]
        
        c(exp(j00),-phi)
        
}

# Wyznacza tylko roznice pomiedzy dopasowaniem a danymi dla zaszumionych parametrow
findErr <- function(pars, barriers, data, temp, sd=0.1/3, nit=1000) {
        
        j00 <- barriers[1]
        phi <- barriers[2]
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        j0 <- pars[1]
        n <- pars[2]
        rs <- pars[3]
        rsh <- pars[4]
        
        ns <- runif(nit, n-n*sd, n+n*sd)
        j0s <- j00*exp(-q*phi/(ns*k*temp))
        rss <- runif(nit, rs-rs*sd, rs+rs*sd)
        rshs <- runif(nit, rsh-rsh*sd, rsh+rsh*sd)
        
        #         ns <- rnorm(nit, n, sd*n)
        #         j0s <- j00*exp(-q*phi/(ns*k*temp))
        #         rss <- rnorm(nit, rs, sd*rs)
        #         rshs <- rnorm(nit, rsh, sd*rsh)
        
        npars <- cbind(j0s, ns, rss, rshs)
        
        store <- data.frame()
        
        for (i in seq(nit)) {
                
                tmpjv <- jvReal(npars[i,], data)
                
                tmpsum <- sum((data[,2]-tmpjv[,2])^2)
                
                store <- rbind(store, c(npars[i,], tmpsum))
                
                print(c(i, npars[i,]))
                
        }
        
        names(store) <- c("j0", "n", "rs", "rsh", "err")
        store
        
}