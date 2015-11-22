# Wczytywanie danych, tworzy liste z danymi dla kazdej temperatury
readData <- function(path, pattern) {
        
        #oldwd <- getwd()
        
        #setwd(path)
        
        files <- list.files(path, pattern, full.names = TRUE)
        
        
        # sortowanie po temepraturze
        indexsort <- order(as.numeric(gsub("(.*)/cf_|K_(.*)", "", files)))
        
        files <- files[indexsort]
        
        #nfiles <- length(files)
        
        datalist <- list()
        
        for (i in seq_along(files)) {
                
                data <- read.csv(files[i], header = FALSE, sep="\t")
                
                data <- data[c(1,2)]
                
                names(data) <- c("freq", "capacitance")
                
                #data <- data[data[,1]<2e5,]
                
                datalist[[i]] <- data
        }
        
        #setwd(oldwd)
        
        names(datalist) <- gsub("^(.*)/","",files)
        
        datalist
        
}

# Wyrysowanie wszystkich przebiegow C-f
plotAll <- function(data, ylim=c(4e-10, 8e-10)) {
        
        len <- length(data)
        
        plot(data[[1]], log="x", ylim=ylim, type="l", xaxp=c(1e2,1e6,3))
        
        for (i in seq(from=2, to=len, by=1)) {
                
                lines(data[[i]])
                
        }
        
        
}

# Oblicza pochodna dla jednego pliku
getDerivative <- function(data) {
        
        len <- length(data[,1])
        
        der <- data.frame()
        
        for (i in seq(len-2)) {
                
                tmp <- (data[i+2,2]-data[i,2])/(data[i+2,1]-data[i,1])
                
                der <- rbind(der, c(data[i+1,1], abs(data[i+1,1]*tmp)))
                
        }
        
        der <- as.data.frame(der)
        
        names(der) <- c("f", "absfdcdf")
        
        der
        
}

# Oblicza pochodne wszystkich zestawow, wykorzystuje funkcje dla pojednyczego zestawu
derAll <- function(data) {
        
        len <- length(data)
        
        derdata <- list()
        
        for (i in seq(len)) {
                
                der <- getDerivative(data[[i]])
                
                derdata[[i]] <- der
                
        }
        
        derdata
        
}

# poszukuje maximow w wykresach wdC/dw
getMax <- function(data, n=10) {
        
        
        
        data <- data[order(data[,2], decreasing = TRUE),]
        
        data2 <- data[1:n,]
        
        # interpolacja w celu dokladniejszego wyznaczenia maximum
        spl <- spline(data2)
        
        data3 <- cbind(spl$x, spl$y)
        data3 <- as.data.frame(data3)
        names(data3) <- names(data)
        
        len <- data3[,1]
        
        #plot(data, log="x", ylim=c(1e-11,8e-11))
        #lines(data2, type="p", col="green")
        lines(data3, col="green")
        
        # sortowanie po wdC/dw
        data3 <- data3[order(data3[,2], decreasing = TRUE),]
        
        # zwrocenie najwiekszej wartosci sygnalu wdC/dw i odpowiadajacej mu czestotliwosci
        points(data3[1,], pch=16)
        data3[1,]
        
}

# polaczenie temperatur, maximow sygnalu wdC/dw i odpowiadajacych im czestotliwosci
combine <- function(data) {
        
        
        len <- length(data)
        
        datacomb <- data.frame()
        
        for (i in seq(len)) {
                
                # wykorzystanie wyrazenia regularnego do wyznaczenia wartosci temperatur
                temperature <- as.numeric(gsub("cf_|K_(.*)", "", names(data[i])))
                # wyznaczenie pochodnych
                der <- getDerivative(data[[i]])
                # wyznaczenie maximow pochodnych
                maximum <- as.numeric(getMax(der))
                
                # polaczenie temperatur, czestotliwosci i maximow wdC/dw
                datacomb <- rbind(datacomb, c(temperature, maximum))
                
        }
        
        #names(datacomb) <- c("temp", "freq", "max")
        
        
        # dodanie odwrotnosci temperatury
        datacomb <- cbind(datacomb,1/datacomb[,1])
        # czestosc emisji
        datacomb <- cbind(datacomb,datacomb[,2]/2)
        # e/T^2
        datacomb <- cbind(datacomb,log(datacomb[,2]/2/datacomb[,1]^2))
        
        datacomb <- datacomb[order(datacomb[,1]),]
        
        names(datacomb) <- c("temp", "freq", "max", "recipT", "emissionrate", "logemission")
        
        datacomb
        
}

# funkcja wykonujw wszystkie powyzsze czynnosci i zwraca zestaw danych
doAll <- function() {
        
        # podanie sciezki do plikow C-f
        path <- "./CF REL"
        # wyrazenie regularne do wybrania odpowiednich plikow
        pattern <- "*.dat"
        # zakres osi y
        ylim <- c(1e-11, 5e-10)
        
        # zakres czestotliwosci
        fmin <- 5e2
        fmax <- 5e5
        
        # wczytanie danych
        data <- readData(path = path, pattern = pattern)
        
        # obliczenie pochodnych
        der <- derAll(data)
        
        # wyrysowanie pochodnych
        plotAll(data = der, ylim = ylim)
        
        # wektor do wybrania odpowiednich zestawow danych
        sel <- 1:8
        
        # wybranie odpowiednich zestawow danych
        data2 <- data[sel]
        der2 <- der[sel]
        
        # nalozenie ograniczenia czestotliwosci
        for (i in seq_along(sel)) {
                der2[[i]] <- der2[[i]][(der2[[i]][,1]>fmin)&(der2[[i]][,1]<fmax),]
                data2[[i]] <- data2[[i]][(data2[[i]][,1]>fmin)&(data2[[i]][,1]<fmax),]
                lines(der2[[i]], col="red", type="p")
        }
        
        # wyznaczenie parametrow
        datacomb <- combine(data = data2)
        
        datacomb
        
}

# wyznaczenie energi aktywacji z zaleznosci bez tunelowania
activationEnergy <- function(res) {
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        model <- lm(logemission~recipT, res[c(4,5)])
        
        Ea <- -model$coefficients[[2]]*k/q
        slope <- model$coefficients[[1]]
        
        res <- data.frame(Ea, slope)
        names(res) <- c("Ea", "sigma")
        
        res
        
        
}

# wyznaczenie sigma, E, i F w modelu z tunelowaniem
emissionTAT <- function(pars, temp) {
        
        # stale fizyczne
        q <- 1.602177e-19
        k <- 1.3806488e-23
        h <- 6.62607004e-34
        me <- 9.10938356e-31
        effmass <- 0.2 * me
        
        # parametry dopasowania
        sigma <- pars[1]
        Energy <- pars[2]
        Field <- pars[3]
        
        # funkcja podcalkowa
        integrand <- function(x) {
                
                exp(x - x^1.5 * (4*(2*effmass)^0.5*(k*temp)^1.5)/(3*q*h*Field/2/pi))
                
        }
        
        # gorna granica calkowania
        uplim <- q*Energy/(k*temp)
        
        # dolna granica calkowania
        lowlim <- 0
        
        #integral <- integrate(integrand, lower = lowlim, upper = uplim)
        
        # podzial zakresu energii na n odcinkow
        n <- 1000
        energies <- seq(from=0, to=uplim, length.out = n)
        values <- integrand(energies)
        
        # wektor do przechowywania calki
        areas <- vector()
        
        # calkowanie trapezami
        for (i in seq(n-1)) {
                
                tmp <- (values[i+1] + values[i])/2 * (energies[i+1] - energies[i])
                
                areas <- append(areas, tmp)
                
                
        }
        
        # wartosc calki
        integral <- sum(areas)
        
        # predkosc termiczna elektronow
        vth <- sqrt(3*k*temp/effmass)
        
        # efektywna gestosc stanow
        Nc <- 2*(2*pi*effmass*k*temp/h^2)^1.5
        
        # czestosc emisji w modelu bez tunelowania
        e0 <- sigma*Nc*vth*exp(-q*Energy/(k*temp))
        
        # czestosc emisji w modelu z tunelowaniem
        eTAT <- e0*(1+integral)
        
        #log(eTAT/temp^2)
        
        # zwrocenie wartosci czestosci emisji
        
        
        #log(sigma*vth*Nc/temp^2)
        #log(integral)
        
        eTAT
        
}

# wyznaczenie sigma, E, i F w modelu z tunelowaniem
emissionTAT2 <- function(pars, reciptemp) {
        
        # stale fizyczne
        q <- 1.602177e-19
        k <- 1.3806488e-23
        h <- 6.62607004e-34
        me <- 9.10938356e-31
        effmass <- 0.2 * me
        
        temp <- 1/reciptemp
        
        # parametry dopasowania
        c <- pars[1]
        Energy <- pars[2]
        Field <- pars[3]
        
        # funkcja podcalkowa
        integrand <- function(x) {
                
                exp(x - x^1.5 * (4*(2*effmass)^0.5*(k*temp)^1.5)/(3*q*h*Field/2/pi))
                
        }
        
        # gorna granica calkowania
        uplim <- q*Energy/(k*temp)
        
        # dolna granica calkowania
        lowlim <- 0
        
        #integral <- integrate(integrand, lower = lowlim, upper = uplim)
        
        # podzial zakresu energii na n odcinkow
        n <- 1000
        energies <- seq(from=0, to=uplim, length.out = n)
        values <- integrand(energies)
        
        # wektor do przechowywania calki
        areas <- vector()
        
        # calkowanie trapezami
        for (i in seq(n-1)) {
                
                tmp <- (values[i+1] + values[i])/2 * (energies[i+1] - energies[i])
                
                areas <- append(areas, tmp)
                
                
        }
        
        # wartosc calki
        integral <- sum(areas)
        
        c - q*Energy/(k*temp) + log(1+integral)
        
}

# wyliczenie czestosci emisji w zaleznosci od temperatury
emissionRates <- function(pars, temps) {
        
        rates <- sapply(temps, function(x) {emissionTAT(pars = pars, temp = x)})
        
        data <- cbind(temps, rates)
        
        data <- as.data.frame(data)
        
        
        
        #names(data) <- c("recipT", "logemissionrate")
        names(data) <- c("temp", "emissionrate")
        
        data
        
}

# wyliczenie czestosci emisji w zaleznosci od temperatury
emissionRates2 <- function(pars, reciptemps) {
        
        logrates <- sapply(reciptemps, function(x) {emissionTAT2(pars = pars, reciptemp = x)})
        
        data <- cbind(reciptemps, logrates)
        
        data <- as.data.frame(data)
        
        names(data) <- c("recipT", "logemissionrate")
        
        data
        
}

# funkcja do minimalizacji
diff <- function(pars, res) {
        
        temps <- res[,1]
        
        rates <- emissionRates(pars, temps)
        
        #sum((log(res[,5]/res[,1]^2) - log(rates[,2]/res[,1]^2))^2)
        #sum((res[,5] - rates[,2])^2)
        sum((res[,2] - rates[,2])^2/res[,2])
        
}

diff2 <- function(pars, res) {
        
        reciptemps <- res[,1]
        
        logrates <- emissionRates2(pars, reciptemps)
        
        #sum((log(res[,5]/res[,1]^2) - log(rates[,2]/res[,1]^2))^2)
        #sum((res[,5] - rates[,2])^2)
        sum((res[,2] - logrates[,2])^2)
        
}