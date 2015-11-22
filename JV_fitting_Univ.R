# Wyciaga temperatury z plikow, file - sciezka do pliku, pos - pozycja temperatury w nazwie
getTemperature <- function(file, pos) {
        
        #print(file)
        
        splitted <- strsplit(file, split = "_")
        #print(splitted)
        
        temp <- splitted[[1]][pos]
        #print(temp)
        
        temp1 <- strsplit(temp, split="T")
        #print(temp1)
        
        temp2 <- strsplit(temp1[[1]][2], split=" ")
        #print(temp2)
        
        len <- length(temp2[[1]])
        
        temp2[[1]][len]
}

# Przygotowanie danych, obcina ostatnie 10% punktow co do wartosci pradu
# ta wersja funkcji jest lepsza, k - poziom obciecia, cutlin - parametr odciecia nieliniowej czesci rs
preparedata2 <- function(data, Vmin=0.1, k=0.99, cutlin=0.00001) {
        
        # wczytanie danych
        data <- read.csv(data, sep="\t")
        
        names(data) <- c("v", "i", "j")
        
#         smoothedI <- smooth.spline(data[c(1,2)])
#         smoothedJ <- smooth.spline(data[c(1,3)])
        
#         datasmI <- as.data.frame(cbind(smoothedI$x, smoothedI$y))
#         datasmJ <- as.data.frame(cbind(smoothedJ$x, smoothedJ$y))
#         names(datasmI) <- c("v", "i")
#         names(datasmJ) <- c("v", "j")
        
        # minimalne wartosci I i J
        minI <- min(data[,2])
        minJ <- min(data[,3])
        
        # sprawdza czy IV jest jasny czy ciemny, jezeli jasny to przesuwa w gore o Isc, Jsc
        if (minI < 0) {
                
                #data[,2] <- data[,2] + abs(minI)
        
                #data[,3] <- data[,3] + abs(minJ)
                
#                 shiftI <- lm(i~v, datasmI[datasmI[,1]<=Vmin,c(1,2)])
#                 
#                 shiftJ <- lm(j~v, datasmJ[datasmJ[,1]<=Vmin,c(1,2)])
                
                shiftI <- lm(i~v, data[data[,1]<=Vmin,c(1,2)])
                
                shiftJ <- lm(j~v, data[data[,1]<=Vmin,c(1,3)])
                
                
                data[,2] <- data[,2] - shiftI$coefficients[[1]]
                
                data[,3] <- data[,3] - shiftJ$coefficients[[1]]
        }
        
        # usuwa prady i napiecia mniejsze badz rowne 0
        data <- data[data[,1]>0,]
        
        data <- data[data[,2]>0, ]
        
        # maksymalna wartosc pradu
        maxCurrent <- max(data[,2])
        
        # obcina ostatnie 1-k pradu
        data <- data[data[,2]<(k*maxCurrent),]
        
        #data <- abs(data)
        
        
        
        logdata <- cbind(data[,1], log(data[,2]))
        logdata <- as.data.frame(logdata)
        names(logdata) <- c("v", "logi")
        
        # obcina czesc liniowa w skali pollogarytmicznej
        for (i in 1:100) {
                
                res <- lm(logi~v, tail(logdata, i+2))
                qual <- sum((res$residuals)^2)
                
                index <- i+2
                
                if (qual>cutlin) {break}
                
        }
        
        #plot(data[1:2], log="y")
        #lines(cbind(data[,1],exp(data[,1]*res$coeff[2]+res$coeff[1])))
        #points(tail(data[1:2], index), col="blue")
        
        len <- length(data[,1])
        
        #lines(data[1:(len-index),],col="red", type="p")
        
        data[1:(len-index),]
        
}

# JV tylko z j0, n, rs jako zmienne
jv1 <- function(vec, rsh, rsh2, alpha, data, temp) {
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        j0 <- vec[1]
        n <- vec[2]
        rs <- vec[3]
        #rsh <- vec[4]
        
        v <- data[,1]
        j <- data[,2]
        
        jv <- data.frame()
        
        for (i in seq_along(v)) {
                jReal <- j0 * (exp(q * (v[i] - j[i] * rs) / (n * k * temp)) - 1) + (v[i] - j[i] * rs)/rsh + v[i]^alpha/rsh2
                jv <- rbind(jv, c(v[i], jReal))
        }
        names(jv) <- c("v", "j")
        jv
}

# JV tylko z rsh, rsh2, alpha jako zmienne
jv2 <- function(vec, j0, n, rs, data, temp) {
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        #j0 <- vec[1]
        #n <- vec[2]
        #rs <- vec[1]
        rsh <- vec[1]
        rsh2 <- vec[2]
        alpha <- vec[3]
        
        v <- data[,1]
        j <- data[,2]
        
        jv <- data.frame()
        
        for (i in seq_along(v)) {
                jReal <- j0 * (exp(q * (v[i] - j[i] * rs) / (n * k * temp)) - 1) + (v[i] - j[i] * rs)/rsh + v[i]^alpha/rsh2
                jv <- rbind(jv, c(v[i], jReal))
        }
        names(jv) <- c("v", "j")
        jv
}

# Optymalizacha jv1
diff1 <- function(vec, rsh, rsh2, alpha, data, temp) {
        
        jv <- jv1(vec, rsh, rsh2, alpha, data, temp = temp)
        
        sum((log(data[,2]) - log(jv[,2]))^ 2)
        #sum(abs((data[,2]) - (jv[,2]) )/jv[,2])
        
}

# optymalizacja jv2
diff2 <- function(vec, j0, n, rs, data, temp) {
        
        jv <- jv2(vec, j0, n, rs, data, temp = temp)
        
        sum((log(data[,2]) - log(jv[,2]) )^ 2)
        #sum(abs((data[,2]) - (jv[,2]) )/jv[,2])
        
}

# JV z wszystkimi parametrami
jvReal <- function(vec, data, temp) {
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        j0 <- vec[1]
        n <- vec[2]
        rs <- vec[3]
        rsh <- vec[4]
        rsh2 <- vec[5]
        alpha <- vec[6]
        
        v <- data[,1]
        j <- data[,2]
        
        jv <- data.frame()
        
        for (i in seq_along(v)) {
                jReal <- j0 * (exp(q * (v[i] - j[i] * rs) / (n * k * temp)) - 1) + (v[i] - j[i] * rs)/rsh + v[i]^alpha/rsh2
                jv <- rbind(jv, c(v[i], jReal))
        }
        names(jv) <- c("v", "j")
        jv
}

# Funkcja okresla poczatkowe rs i rsh z dopasowania liniowego, rsh i rs z zaleznosci liniowej
# Pobiera dane, dolny zakres napiec do okreslenia rsh, gorny zakres napiec do okreslenia rs
findres <- function(data, lowerVlim=0.1, upperVlim=0.1) {
        
        v <- data[,1]
        #j <- data[,2]
        
        data <- data[1:2]
        
        #logv <- log(v)
        #logj <- log(j)
        
        vmax <- max(v)
        vmin <- min(v)
        
        rshuntData <- data[data[,1]<=(vmin+lowerVlim),]
        rseriesData <- data[data[,1]>=(vmax-upperVlim),]
        
        #rshuntDataLog <- log(rshuntData)
        
        
        names(rseriesData) <- c("vseries", "jseries")
        
        names(rshuntData) <- c("vshunt", "jshunt")
        #names(rshuntDataLog) <- c("vshuntLog", "jshuntLog")
        
        #plot(rshuntData)
        
        resShunt <- lm(jshunt~vshunt-1, rshuntData)
        #resShuntLog <- lm(jshuntLog~vshuntLog-1, rshuntDataLog)
        resSeries <- lm(jseries~vseries, rseriesData)
        
        
        rs <- 1/resSeries$coefficients[[2]]
        
        #rsh <- 1/exp(resShuntLog$coefficients[[1]])
        rsh <- 1/resShunt$coefficients[[1]]
        
        #c(rs, rsh, resShunt$coefficients[[1]])
        c(rs, rsh)
        
}

# usuwa z danych punkty o ujemnych wartoscia pradu oraz punkty z nimi sasiadujace
cleanup <- function(data) {
        
        len <- length(data[,1])
        
        # wartosci mniejsze zero
        isNegative <- (data[,2] <= 0)
        
        # sprawdza czy pierwszy i ostatni punkt sa mniejsze badz rowne zero
        if (data[1,2] <= 0) {isNegative[1] <- TRUE}
        if (data[len,2] <= 0) {isNegative[len] <- TRUE}
        
        # skanuje punkty i zamienia je i sasiadow na mniejsze rowne 0
        for (i in seq(from=2, to=(len-1))) {
                
                if (data[i,2]<=0) {
                        
                        isNegative[i-1] <- TRUE
                        isNegative[i] <- TRUE
                        isNegative[i+1] <- TRUE
                        
                }
                
        }
        
        isNegative
        
        data[!isNegative,]
        
}

# Wyznaczanie poczatkowego J0 i n
# Pobiera dane, temperature, dolny i gorny zakres napiec do okreslenia rsh i rs
# zakres napiec do okreslenia rsh2 i alpha
initialparams <- function(data, temp, lowerVlim=0.1, upperVlim=0.1, vRange=c(0.2, 0.4)) {
        
#         v <- data[,1]
#         j <- data[,2]
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        # wyznacza poczatkowe rs i rsh
        res1 <- findres(data, lowerVlim, upperVlim)
        
        rsh <- res1[2]
        rs <- res1[1]
        
        data <- cleanup(data[1:2])
        
        plot(data[1:2], log="y")
        
        # rysuje jshunt(rsh) w funkcji v
        lines(cbind(data[,1],data[,1]/rsh), col="red", lty="dashed")
        
        # odejmuje u*rsh, usuwa punkty mniejsze rowne zero
        index2 <- data[,1] > lowerVlim
        data2 <- cbind(data[index2,1],data[index2,2]-data[index2,1]/rsh)
        data2 <- cleanup(data2)
#         data2 <- smooth.spline(data2)
#         data2 <- cbind(data2$x, data2$y)
        data2 <- as.data.frame(data2)
        names(data2) <- c("v", "i")
        
        
        # rysuje jv po odjeciu rsh
        lines(data2, col="red")
        
        
        # wyznacza rsh2 i alpha
        res2 <- lm(i~v, log(data2[(data2[,1]>=vRange[1]) & (data2[,1]<=vRange[2]),]))
        #res2 <- lm(i~v, log(data2[data2[,1]<=vRange[2],]))
        rsh2 <- exp(-res2$coeff[[1]])
        alpha <- res2$coeff[[2]]
        
        # rysuje jshunt(rsh, rsh2, alpha) w funkcji napiecia
        lines(cbind(data[,1], data[,1]/rsh+data[,1]^alpha/rsh2), col="green", lty="dashed")
        
        # odejmuje u^alpha*rsh2, usuwa punkty mniejsze rowne zero
        index3 <- data[,1] > vRange[2]
        data3 <- cbind(data[index3,1],data[index3,2] - data[index3,1]/rsh - data[index3,1]^alpha/rsh2)
        data3 <- cleanup(data3)
#         data3 <- smooth.spline(data3)
#         data3 <- cbind(data3$x, data3$y)
        data3 <- as.data.frame(data3)
        names(data3) <- c("v", "i")
        
        
        # rysuje jv po odjeciu wplywu rsh i rsh2 z alpha
        lines(data3, col="green")
        
        #print(findres(data3, 0.1, 0.1))

        # dane w skali pollogarytmicznej do celow okreslenia j0 i n
        logdata3 <- cbind(data3[,1],log(data3[,2]))
        logdata3 <- as.data.frame(logdata3)

        names(logdata3) <- c("v", "logi")
        
        med <- median(logdata3[,1])
        
        #print(med)

        #res3 <- lm(logi~v, logdata3)
        res3 <- lm(logi~v, logdata3[logdata3[,1]<med,])
        #res3 <- lm(logi~v, logdata3[logdata3[,1]<(max(data[,1])-0.2),])
        
        # wartosci z dopasowania liniowego
        j0 <- exp(res3$coefficients[[1]])
        
        n <- 1/(res3$coefficients[[2]]*k*temp/q)
        
        # v w funkcji j
        data4 <- cbind(data3[,2], n*k*temp/q*log(data3[,2]/j0))
        
        # zamiana v i j miejscami
        data4 <- cbind(data4[,2], data4[,1])
        
        lines(data4, col="blue")
#         
#         data5 <- cbind(data3[,1]-data4[,1], data4[,2])
#         
#         lines(data5)
#         
#         res2 <- findres(as.data.frame(data5), 0.1, 0.1)
#         
#         print(res2[1])
#         
#         rs <- res2[1]
        
        as.numeric(c(j0, n, rs, rsh, rsh2, alpha))
        
}

# Funkcja poszukuje dokladne wartosci parametrow
# Pobiera wektor parametrow poczatkowych, dane, wartosc temperatury oraz liczbe iteracji
fit <- function(vec, data, temp, nit) {
        
        data <- data[1:2]
        
        params1 <- vec[1:3]
        params2 <- vec[4:6]
        
        #alpha <- vec[5]
        
        # okreslenie poczatkowej jakosci dopasowania z prametrami inicjalizujacymi
        jv <- jvReal(vec, data = data, temp = temp)
        jvsum <- sum((log(data[,2])-log(jv[,2]))^2)
        
        # alpha - reflection, beta - contraction, gamma - expansion
        control <- list(alpha=1.0, beta=0.75, gamma=1.5)
        
        for (i in seq(nit)) {
                
                # wykorzystanie parametrow poczatkowych tylko w pierwszej iteracji
                if (i==1) {
                        
                        # wypisanie parametrow do konsoli
                        print(c(i-1, vec[1], vec[2], vec[3], vec[4], vec[5], vec[6], jvsum))
                        
                        res1 <- optim(params1, diff1, rsh=vec[4], rsh2=vec[5], alpha=vec[6], data=data, temp=temp, control=control)
                        j0_tmp <- res1$par[1]
                        n_tmp <- res1$par[2]
                        rs_tmp <- res1$par[3]
                        res2 <- optim(params2, diff2, j0=j0_tmp, n=n_tmp, rs=rs_tmp, data=data, temp=temp, control=control)
                        rsh_tmp <- res2$par[1]
                        rsh2_tmp <- res2$par[2]
                        alpha_tmp <- res2$par[3]
                        
                } else {
                        res1 <- optim(c(j0,n,rs), diff1, rsh=rsh_tmp, rsh2=rsh2_tmp, alpha=alpha_tmp, data=data, temp=temp, control=control)
                        j0_tmp <- res1$par[1]
                        n_tmp <- res1$par[2]
                        rs_tmp <- res1$par[3]
                        res2 <- optim(c(rsh,rsh2, alpha), diff2, j0=j0_tmp, n=n_tmp, rs=rs_tmp, data=data, temp=temp, control=control)
                        rsh_tmp <- res2$par[1]
                        rsh2_tmp <- res2$par[2]
                        alpha_tmp <- res2$par[3]
                        
                }
                
                
                
                jv_tmp <- jvReal(c(j0_tmp,n_tmp,rs_tmp,rsh_tmp, rsh2_tmp, alpha_tmp), data = data, temp = temp)
                jvsum_tmp <- sum((log(data[,2])-log(jv_tmp[,2]))^2)
                
                # okresla ile wzglednie nowe dopasowanie jest lepsze od poprzedniego
                distance <- abs(jvsum_tmp - jvsum)
                
                # jezeli nowe dopasowanie jest lepisze o 1e-4 to nastepuje podstawienie nowych parametrow
                if ((jvsum_tmp < jvsum) & (distance > 1e-6)) {
                        
                        j0 <- j0_tmp
                        n <- n_tmp
                        rs <- rs_tmp
                        rsh <- rsh_tmp
                        rsh2 <- rsh2_tmp
                        alpha <- alpha_tmp
                        jvsum <- jvsum_tmp
                        
                        # wypisanie parametrow do konsoli
                        print(c(i,j0,n,rs,rsh, rsh2, alpha, jvsum, distance))
                        
#                         plot(data, log="y")
#                         lines(jvReal(c(j0,n,rs,rsh, rsh2, alpha), data = data, temp = temp))
#                         lines(cbind(data[,1], data[,1]/rsh), col="red", lty="dashed")
#                         lines(cbind(data[,1], data[,1]/rsh+data[,1]^alpha/rsh2), col="green", lty="dashed")
                        
                } else {
                        
                        # w przeciwnym razie nastepuje zatrzymanie petli
                        if (i>1) 
                        {
                                
                                print("CONVERGENCE!!")
                                
                        } else {
                                
                                print("FAIL!!")
                                
                        }
                        break
                        
                }
                
                
                
        }
        
        # zwrocenie parametrow
        c(j0,n,rs,rsh, rsh2, alpha, jvsum)
        
}

# Wykonuje wszystkie kroki
doAll <- function() {
        
        if (sum(search()=="package:beepr")==0) {
                
                library(beepr)
                
        }
        
        #############################################
        # SET THIS PARAMETERS FIRST
        
        # pattern for selecting files
        pattern <- "*L40.dat"
        path <- "./IV WHITE REL/"
        
        # number of files to fit
        n <- 11
        
        # position of temperature value in filename
        pos <- 6
        
        # firs values tells to cut last 10% of current
        # second value is used to shift JV by average Jsc at low voltages < Vmin
        # third value is used to remove non-linear part of JV at high voltages
        k <- 0.99
        Vmin <- 0.1
        cutlin <- 0.00001
        
        # voltage limits to estimate initial rsh and rs
        lowerVlim <- 0.1
        upperVlim <- 0.1
        vRange=c(0.2, 0.45)
        
        # number of iterations
        nit <- 30
        
        #############################################
        
        # All files in directory
        files1 <- list.files(path = path, pattern = pattern, full.names = TRUE)
        # Only first n files
        len <- length(files1)
        first <- len-n+1
        files2 <- files1[first:len]
        
        # data frame for storing results
        resAll <- data.frame()
        
        
        for (i in seq(n)) {
                filename <- files2[i]
                print("File name:")
                print(filename)
                
                # gets temperature
                temp <- as.numeric(getTemperature(filename, pos))
                
                # prepares data
                data <- preparedata2(filename, Vmin, k, cutlin)
                #data <- rawdata[[1]]
                area <- mean(data[,2]/data[,3]*1000)
                
                # finds initial parameters
                p <- try(initialparams(data=data, temp = temp, lowerVlim = lowerVlim, upperVlim = upperVlim, vRange = vRange))
                if (class(p)=="try-error") {p <- c(NA,NA,NA,NA,NA,NA)}
                
                print("Initial values")
                print(p)
                print("Fitting results:")
                
                # fits, returns j0, n, rs, rsh and guality factor, simple error handling
                res <- try(fit(vec = p, data = data, temp = temp, nit = nit))
                if (class(res)=="try-error") {res <- c(NA,NA,NA,NA,NA,NA,NA)}
                print("")
                
                v <- c(temp, res, area, 1/temp, log(res[1]*1000/area), res[2]*log(res[1]*1000/area), 1/res[2])
                
                resAll <- rbind(resAll, v)
                
        }

        #Sygnal dzwiekowy na zakonczenie skryptu
        beep(4)
        Sys.sleep(0.5)
        
        names(resAll) <- c("Temperature", "i0", "n", "rs", "rsh", "rsh2", "alpha","Quality", "Area", "recipT", "lnj0", "Alnj0", "recipA")
        resAll
        
}

# Rysowanie wszystkich plikow z katalogu
plotAll <- function(path, pattern) {
        
        files <- list.files(path = path, pattern = pattern, full.names = TRUE)
        
        n <- length(files)
        
        for (i in n:1) {
                
                data <- preparedata2(files[i])
                
                if (i==n) {
                        
                        plot(data[1:2], log="y", type="l", xlim=c(0,1.5), ylim=c(1e-6, 2e-2), xaxp=c(0,1.5,15))
                } else {
                        
                        lines(data[1:2])
                        
                }
                
        }
        
}

# Rysowanie konkretnego zestawu danychi jego dopasowania
plotSpecific <- function(path, res, i) {
        
        data <- preparedata2(path)
        
        plot(data[1:2], log="y",  main=path)
        
        temp <- res[i,1]
        
        p <- as.numeric(res[i,2:7])
        
        lines(jvReal(p, data, temp))
        
        
}

# Wyeksportowanie dopasowanych parametrow do pliku, dorobic sprawdzanie czy plik istnieje
exportData <- function(res, path) {
        
        write.table(x = res, file = path, append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
        
}

# Import dopasowanych parametrow
importData <- function(path) {
        
        
        read.table(path, header = TRUE)
        
}

# Funkcja do oszacowania niepewnosci
estimErrDirection <- function(pars, data, temp, sel, n, err=0.02) {
        
        set.seed(100)
        
        # jv wyliczony z dopasowanych parametrow
        jv <- jvReal(as.numeric(pars), data, temp)
        
        # odchylenie standardowe, wykorzystane w chi^2
        sd <- sqrt(sum((log(data[,2]) - log(jv[,2]))^2)/length(data[,1]))
        
        # n-krotne zaszumienie wybranych parametrow najlepszego dopasowania
        values <- sapply(pars[sel], function(x) {rnorm(n, x, err*x)})
        
        # przeksztalcenie wyrbanych zaszumionych parametrow do ramki danych
        values <- as.data.frame(values)
        
        # dodanie nazw kolumn w postaci P1, P2, P3, ...
        #names(values) <- sapply(sel, function(x) {paste(c("P", x), collapse = "")})
        
        # replikacja wszystkich parametrow najlepszego dopasowania
        # powstala macierz musi zostac transponowana bo replicate automatycznie transponuje
        #reppars <- t(replicate(n, pars, simplify = TRUE))
        reppars <- apply(pars,2,function(x) rep(as.numeric(x),n))
        
        # zreplikowane parametry jako ramka danych
        reppars_df <- as.data.frame(reppars)
        
        # podstawienie zaszumionych paramerow w miejsce zreplikowanych
        # od tego moment dysponujemy zestawem i0, n, rs, rsh, rsh2, alpha gdzie wyrbane parametry sa zaszumione
        reppars_df[,sel] <- values
        
        # obliczenie wartosci chi^2 dla kazdego zestawu zaszumionych parametrow
        chi_values <- apply(reppars_df, 1, function(x) {jv <- jvReal(as.numeric(x), data, temp); sum(((log(jv[,2])-log(data[,2]))/sd)^2)})
        
        # liczba stopni swobody
        df <- length(data[,1]) - length(sel) - 1
        
        # wyliczenie wartosci krytycznej kwantyla chi^2
        chi_crit <- qchisq(0.95, df)
        
        # indeksy zaakceptowanych wartosci
        indexes <- chi_values <= chi_crit
        
        # zrobienie wykresu zaakceptowanych wartosci
        reppars_df_acc <- reppars_df[indexes,sel]
        plot(reppars_df_acc)
        
        # punkty zaakceptowane w skali n~log(i0), ktore spelniaja warunek <= chi krytyczne
        # logarytmowanie bo i0~exp(n)
        logvalues_acc <- cbind(log(reppars_df_acc[,1]), reppars_df_acc[,2])
        logvalues_acc <- as.data.frame(logvalues_acc)
        names(logvalues_acc) <- names(reppars_df_acc[,sel])
        
        # tylko zaakceptowane parametry, zamienione do ramki danycyh z wartosciami numerycznymi
        #reppars_df_acc <- as.data.frame(t(apply(reppars_df[indexes,], 1, as.numeric)))
        
        # stworzenie formuly do fitowania
        formula <- as.formula(paste(rev(names(logvalues_acc)), collapse = "~"))
        
        # model dopasowania
        model <- lm(formula, logvalues_acc)
        intercept <- model$coefficients[[1]]
        slope <- model$coefficients[[2]]
        
        # model danych zbudowany na podstawie dopasowania
        modelled_data <- cbind(exp(logvalues_acc[,1]),logvalues_acc[,1]*slope + intercept)
        
        # posortowanie w celu lepszego rysowania
        modelled_data <- modelled_data[order(modelled_data[,1]),]
        
        # dodanie linii do wykresu z punktami
        lines(modelled_data)
        
        # stworznie ramki danych dla parametru
        modelparams <- data.frame(intercept=intercept, slope=slope)
        
        # zwrocenie ramki danych
        modelparams
        
}

estimErr <- function(pars, data, temp, directions, sel, n, err=0.1, k) {
        
        
        # wyliczenie jv
        jv <- jvReal(as.numeric(pars), data, temp)
        
        # odchylenie standardowe, wykorzystane w chi^2
        sd <- sqrt(sum((log(jv[,2]) - log(data[,2]))^2)/length(data[,1]))
        
        # zakres skanowania
        min <- as.numeric(pars[sel[1]] - err * pars[sel[1]])
        max <- as.numeric(pars[sel[1]] + err * pars[sel[1]])
        
        # szerokosci przedzialow do skanowania
        delta <- as.numeric((max-min)/n)
        
        # wektor parametru x
        par1 <- seq(from=min, to=max, by=delta)
        
        
        # wektor parametru y
        par2 <- log(par1)*directions$slope + directions$intercept
        
        
        
        # stworzenie macierzy z x i y
        values <- cbind(par1, par2)
        
        # przeksztalcenie do ramki danych
        values <- as.data.frame(values)
        
        # replikacja parametrow najlepszego dopasowania
        #reppars <- t(replicate(n+1, pars, simplify=TRUE))
        reppars <- apply(pars,2,function(x) rep(as.numeric(x),n+1))
        
        # przeksztalcenie do ramki danych
        reppars_df <- as.data.frame(reppars)
        
        # podstawienie parametrow skanowanych
        reppars_df[,sel] <- values
        
        # chi min
        chi_min <- sum((log(jv[,2]) - log(data[,2]))^2/sd^2)
        
        # wyliczenie wartosci chi^2
        chi_values <- apply(reppars_df, 1, function(x) {jv <- jvReal(as.numeric(x), data, temp); sum(((log(jv[,2])-log(data[,2]))/sd)^2)})
        print(pars)
        print(chi_min)
        print(reppars_df[chi_values==min(chi_values),])
        print(min(chi_values))
        
        plot(cbind(par2,chi_values), col="red", log="y")
        
        # liczba stopni swobody
        df <- length(data[,1]) - length(sel) - 1
        chi_crit <- qchisq(0.95, df)
        print(chi_crit)
        
        
        # indeksy parametrow spelniajacych warunek
        indexes <- chi_values <= chi_crit
        
        # wybor zaakceptowanych parametrow
        accepted <- as.data.frame(cbind(reppars_df[indexes,], chi_values[indexes]))
        
        print(accepted)
        
        # przeksztalcenie wszystkich wartosci do postaci numerycznej
        #accepted <- as.data.frame(t(apply(accepted, 1, as.numeric)))
        
        #accepted <- apply(accepted, 2, as.numeric)
        
        # wyliczenie niepwenosci z parametrow pierwszy i ostatnich zaakceptowanych
        upper_err <- abs(head(accepted[sel], 1) - pars[sel])
        lower_err <- abs(tail(accepted[sel], 1) - pars[sel])
        
        # zwrocenie jako wektor z niepewnosciami
        error1 <- max(upper_err[1], lower_err[1])
        error2 <- max(upper_err[2], lower_err[2])
        
        errors <- cbind(error1, error2)
        
        errors <- as.data.frame(errors)
        
        errors
        
}