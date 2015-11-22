prepareData <- function(path, pattern="*.dat") {
        
        files <- list.files(path=path, pattern=pattern, full.names = TRUE)
        
        parameters <- data.frame()
        
        for (i in 1:length(files)) {
                
                data <- read.table(files[i], header=FALSE, col.names = c("v", "i", "j"))
                
                # data in the 1st quadrant of coordinates system
                data1q <- data[data$i>0,]
                
                # data in the 4th quadrant of coordinates system
                data4q <- data[data$v>=0 & data$i<=0,]
                
                # wyznacza voc z dopasowania liniowego do 4 punktow w poblizu 0
                datavoc <- data.frame(rbind(tail(data4q, 2), head(data1q, 2)))
                lmfit <- lm(i~v, datavoc)
                intercept <- lmfit$coefficients[[1]]
                slope <- lmfit$coefficients[[2]]

                voc <- -round(intercept/slope, 5)
                isc <- abs(data4q[1,2])
                jsc <- abs(data4q[1,3])
                logisc <- log(isc)
                logjsc <- log(jsc)
                
                # usuwa wszystkie spacje ze sciezki pliku
                filename <- gsub(" ", "", files[i])
                #print(filename)
                
                # pozycja temperatury w sciezke do pliku
                tpos <- regexpr("_T[0-9]", filename)
                
                # pozycja mocy swiatle w sciezce do pliku
                lpos <- regexpr("_L[0-9]", filename)
                
                # pozycja rozszerzenie nazwy pliku
                epos <- regexpr(".dat", filename)
                
                temperature <- as.numeric(substr(filename, tpos[1]+2, lpos[1]-1))
                light <- as.numeric(substr(filename, lpos[1]+2, epos[1]-1))
                
                print(c(voc, isc, jsc, temperature, light))
                parameters <- rbind(parameters, c(voc, isc, logisc, jsc, logjsc, temperature, light))
                
        }
        
        #parameters <- as.data.frame(parameters)
        names(parameters) <- c("voc", "isc", "logisc","jsc", "logjsc", "temperature", "light")
        parameters <- parameters[parameters$light!=0,]
        parameters <- parameters[order(parameters$temperature, parameters$light),]
        parameters
        
}

diff <- function(rsh, data, delta) {
        
        voc <- data[,1]
        isc <- data[,2]+delta
        
        #x <- voc
        y <- isc - voc/rsh
        
        #print(y)
        
        ispositive <- y > 0
        
        #print(sum(ispositive))
        
        logy <- log(y[ispositive])
        
        lmfit <- lm(logy~x, data.frame(logy=logy, x=voc[ispositive]))
        
        qual <- lmfit$residuals
        
        sum(abs(qual))/length(qual)
        
}

getParams <- function(data, tlim1=250, tlim2=200, skip=2, drop=3, vlim=c(0.4,0.9)) {
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        light <- unique(data$light)
        temperature <- unique(data$temperature)
        
        params <- data.frame()
        
        len <- length(light)
        
        for (i in light[(skip+1):(len-drop)]) {
                
                #wyznaczenie energii aktywacji dla kazdego z natezen
                
                tmp <- data[(data$temperature>=tlim1)&(data$light==i),]
                
                lmfit <- lm(voc~temperature, tmp)
                
                Ea <- lmfit$coefficients[[1]]
                
                #print(c(i, intercept))
                
                #dataEa <- rbind(dataEa, c(i, Ea))
        
        
        # wyznaczenie j0, n dla kazdej z temperatur
        for (t in temperature[temperature >= tlim2]) {
                
                tmp2 <- data[data$temperature==t,c(5,1)]
                
                #len2 <- length(tmp2[,1])
                
                #print(len2)
                
                if (t==temperature[temperature >= tlim2][1]) {
                        plot(tmp2, ylim=vlim)
                } else {
                        points(tmp2)
                }
                
                tmp3 <- tmp2[(skip+1):(len-drop),]
                
                points(tmp3, col="green")
                
                lmfit2 <- lm(voc~logjsc, tmp3)
                
                abline(lmfit2)
                
                intercept <- lmfit2$coefficients[[1]]
                slope <- lmfit2$coefficients[[2]]
                
                A <- (Ea-intercept)
                n <- slope*q/k/t
                logj00 <- A*q/k/t/n
                j00 <- exp(logj00)
                j0 <- j00*exp(-q*Ea/(n*k*t))
                
                #print(c(t, i, j00, n, logj00, Ea))
                
                params <- rbind(params, c(t, i, j00, logj00, n, j0, Ea))
                
        }
                
                
                
        }
        names(params) <- c("temperature", "light", "j00", "logj00", "n", "j0", "Ea")
        params
        
}