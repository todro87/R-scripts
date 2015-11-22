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
        
        
        # posortowanie od najmniejszej do najwiekszej wartosci
        data <- data[order(data[,2], decreasing = TRUE),]
        
        # wybranie pierwszych n najwiekszych wartosci
        data2 <- data[1:n,]
        
        # interpolacja w celu dokladniejszego wyznaczenia maximum
        spl <- spline(data2)
        
        data3 <- cbind(spl$x, spl$y)
        data3 <- as.data.frame(data3)
        names(data3) <- names(data)
        
        len <- data3[,1]
        
        #plot(data, log="x", ylim=c(1e-11,8e-11))
        lines(data2, type="p", col="blue")
        lines(data3, col="green")
        
        # sortowanie po wdC/dw
        data3 <- data3[order(data3[,2], decreasing = TRUE),]
        
        points(data3[1,], pch=16)
        
        # zwrocenie najwiekszej wartosci sygnalu wdC/dw i odpowiadajacej mu czestotliwosci
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
        ylim <- c(1e-11, 7e-10)
        
        # zakres czestotliwosci
        fmin <- 2e3
        fmax <- 4e5
        
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

# wyznaczenie log(e/T^2) w funkcji c, E i F
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

# wyliczenie e/T^2 w zaleznosci od odwrotnosci temperatury
emissionRates2 <- function(pars, reciptemps) {
        
        logrates <- sapply(reciptemps, function(x) {emissionTAT2(pars = pars, reciptemp = x)})
        
        data <- cbind(reciptemps, logrates)
        
        data <- as.data.frame(data)
        
        names(data) <- c("recipT", "logemissionrate")
        
        data
        
}

# funkcja do minimalizacji
diff2 <- function(pars, res) {
        
        reciptemps <- res[,1]
        
        logrates <- emissionRates2(pars, reciptemps)
        
        #sum((log(res[,5]/res[,1]^2) - log(rates[,2]/res[,1]^2))^2)
        #sum((res[,5] - rates[,2])^2)
        sum((res[,2] - logrates[,2])^2)
        
}

fit <- function(pars, res) {
        
        res <- res[c(4,6)]
        
        res <- optim(pars, diff2, res=res)
        
        res
        
}

# Funkcja do oszacowania niepewnosci
estimErrDirection <- function(pars, res, sel, n, err=0.02) {
        
        set.seed(100)
        
        # podzbior zbioru z wynikami dopasowania
        subres <- res[c(4,6)]
        
        # wyliczenie 1/T i log e/T^2
        logrates <- emissionRates2(pars, subres[,1])
        
        # odchylenie standardowe, wykorzystane w chi^2
        sd <- sqrt(sum((logrates[,2] - subres[,2])^2)/length(subres[,1]))
        
        # n wartosci zaszumionych parametrow
        values <- sapply(pars[sel], function(x) {rnorm(n, x, err*x)})
        
        # przeksztalcenie do ramki danych
        values <- as.data.frame(values)
        
        # dodanie nazw kolumn w postaci P1, P2, P3
        names(values) <- sapply(sel, function(x) {paste(c("P", x), collapse = "")})
        
        # replikacja parametrow najlepszego dopasowania, powstala macierz jest transponowana
        reppars <- t(replicate(n, pars, simplify = TRUE))
        
        # zreplikowane parametry jako ramka danych
        reppars_df <- as.data.frame(reppars)
        
        # podstawienie zaszumionych paramerow w miejsce zreplikowanych
        reppars_df[,sel] <- values
        
        # przeksztalcenie ramki danych na liste, umozliwia to uzywanie operacji lapply
        #reppars_list <- as.list(as.data.frame(t(reppars_df)))
        
        # wyliczenie kwantyli rozkladu chi^2
        #chi_values <- lapply(reppars_list, function(x) {logrates <- emissionRates2(x, 1/subres[,1]); sum(((logrates[,2]-subres[,2])/sd)^2)})
        # w ten sposob nie trzeba tworzyc list, 1 mowi ze obliczamy dla wierszy
        chi_values <- apply(reppars_df, 1, function(x) {logrates <- emissionRates2(x, subres[,1]); sum(((logrates[,2]-subres[,2])/sd)^2)})
        
        # dla kazdego wiersza z zaszumionymi parametrami dodana jest wartosc kwantyli chi^2
        #pars_chi <- cbind(reppars_df, cbind(chi_values))
        
        # liczba stopni swobody
        df <- length(subres[,1]) - length(sel) - 1
        
        # wyliczenie wartosci krytycznej kwantyla chi^2
        chi_crit <- qchisq(0.95, df)
        
        # zaakceptowane wartosci
        #accepted <- pars_chi[pars_chi[,4]<chi_crit,]
        
        # indeksy zaakceptowanych wartosci
        indexes <- chi_values <= chi_crit
        
        # zrobienie wykresu punktow, ktore spelniaja warunek <= chi krytyczne
        plot(values[indexes,])
        
        # stworzenie formuly do fitowania
        formula <- as.formula(paste(rev(names(values)), collapse = "~"))
        
        # model dopasowania
        model <- lm(formula, values[indexes,])
        
        # model danych zbudowany na podstawie dopasowania
        modelled_data <- cbind(values[indexes,1],values[indexes,]*model$coefficients[2] + model$coefficients[1])
        
        # posortowanie w celu lepszego rysowania
        modelled_data <- modelled_data[order(modelled_data[,1]),]
        
        # dodanie linii do wykresu z punktami
        lines(modelled_data)
        
        # stworznie ramki danych dla parametru
        modelparams <- data.frame(intercept=model$coefficients[[1]], slope=model$coefficients[[2]])
        
        # zwrocenie ramki danych
        modelparams
        
}

estimErr <- function(pars, directions, sel, n, err=0.05, k, res) {
        
        # podzbior zbioru z wynikami dopasowania
        subres <- res[c(4,6)]
        
        # wyliczenie 1/T i log e/T^2
        logrates <- emissionRates2(pars, subres[,1])
        
        # odchylenie standardowe, wykorzystane w chi^2
        sd <- sqrt(sum((logrates[,2] - subres[,2])^2)/length(subres[,1]))
        
        # zakres skanowania
        min <- pars[sel[1]] - err * pars[sel[1]]
        max <- pars[sel[1]] + err * pars[sel[1]]
        
        # szerokosci przedzialow do skanowania
        delta <- (max-min)/n
        
        # wektor parametru x
        par1 <- seq(from=min, to=max, by=delta)
        
        # wektor parametru y stworzony z kierunku zmian bledow (?????)
        par2 <- par1*directions$slope + directions$intercept
        
        # stworzenie macierzy z x i y
        values <- cbind(par1, par2)
        
        # przeksztalcenie do ramki danych
        values <- as.data.frame(values)
        
        # replikacja parametrow najlepszego dopasowania
        reppars <- t(replicate(n+1, pars, simplify=TRUE))
        
        # przeksztalcenie do ramki danych
        reppars_df <- as.data.frame(reppars)
        
        # podstawienie parametrow skanowanych
        reppars_df[,sel] <- values
        
        # chi min
        chi_min <- sum((logrates[,2] - subres[,2])^2/sd^2)
        
        #print(chi_min)
        
        # wyliczenie wartosci chi^2
        chi_values <- apply(reppars_df, 1, function(x) {logrates <- emissionRates2(x, subres[,1]); sum(((logrates[,2]-subres[,2])/sd)^2)})
        
        # liczba stopni swobody
        df <- length(subres[,1]) - length(sel) - 1
        
        # wyliczenie wartosci krytycznej kwantyla chi^2
        chi_crit <- qchisq(0.95, df)
        
        #print(chi_crit)
        
        # indeksy parametrow spelniajacych warunek
        indexes <- chi_values <= (chi_min + k)
        
        # wybor zaakceptowanych parametrow
        accepted <- cbind(reppars_df[indexes,], chi_values[indexes])
        
        names(accepted) <- c("P1", "E", "F", "Chi2")
        
        upper_err <- abs(head(accepted[1:3], 1) - f$par)
        lower_err <- abs(tail(accepted[1:3], 1) - f$par)
        
        # wyliczenie niepwenosci z parametrow pierwszy i ostatnich zaakceptowanych
        
        
        # zwrocenie jako wektor z niepewnosciami
        error1 <- max(upper_err[1], lower_err[1])
        error2 <- max(upper_err[2], lower_err[2])
        error3 <- max(upper_err[3], lower_err[3])
        
        errors <- cbind(error1, error2, error3)
        
        errors <- as.data.frame(errors)
        
        names(errors) <- c("uP1", "uE", "uF")
        
        errors
        
}

# Funkcja do wyrysowania e/T^2(1/T), najlepszego dopasowania i zakresu + i - 3sigma
plotFit <- function(pars, res) {
        
        subres <- res[c(4,6)]
        
        plot(subres)
        
        logrates <- emissionRates2(pars, subres[,1])
        
        lines(logrates)
        
        sd <- sqrt(sum((logrates[,2] - subres[,2])^2)/length(subres[,1]))
        
        lines(cbind(logrates[,1], logrates[,2]+sd), col="red")
        
        lines(cbind(logrates[,1], logrates[,2]-sd), col="red")
        
        
}