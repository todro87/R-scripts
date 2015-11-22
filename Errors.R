# Funkcja do oszacowania niepewnosci
estimErrDirectionDT <- function(pars, data, temp, sel, n, err=0.02) {
        
        set.seed(100)
        
        q <- 1.602177e-19
        k <- 1.3806488e-23
        
        # jv wyliczony z dopasowanych parametrow
        jv <- jvReal(as.numeric(pars), data, temp)
        
        # odchylenie standardowe, wykorzystane w chi^2
        sd <- sqrt(sum((data[,2] - jv[,2])^2)/length(data[,1]))
        
        # n-krotne zaszumienie wybranych parametrow najlepszego dopasowania
        values <- data.table(sapply(pars[sel], function(x) {rnorm(n, x, err*x)}))
        
        # przeksztalcenie wyrbanych zaszumionych parametrow do ramki danych
        #values <- as.data.frame(values)
        
        # dodanie nazw kolumn w postaci P1, P2, P3, ...
        #names(values) <- sapply(sel, function(x) {paste(c("P", x), collapse = "")})
        
        # replikacja wszystkich parametrow najlepszego dopasowania
        # powstala macierz musi zostac transponowana bo replicate automatycznie transponuje
        reppars <- data.table(apply(pars,2,function(x) rep(as.numeric(x),n)))
        
        # zreplikowane parametry jako ramka danych
        #reppars_df <- as.data.frame(reppars)
        
        # podstawienie zaszumionych paramerow w miejsce zreplikowanych
        # od tego moment dysponujemy zestawem i0, n, rs, rsh, rsh2, alpha gdzie wyrbane parametry sa zaszumione
        rndpars <- cbind(values[,.(i0,n)], reppars[,.(rs,rsh,rsh2,alpha)])

        # obliczenie wartosci chi^2 dla kazdego zestawu zaszumionych parametrow
        #chi_values <- apply(rndpars, 1, function(x) {jv <- jvReal(as.numeric(x), data, temp); sum(((jv[,2]-data[,2])/sd)^2)})
         
        #rndpars[2,chivals:=sum(((i0*exp(q*(data[,1]-data[,2]*rs)/(n*k*temp))+rsh/data[,1]+data[,1]^alpha/rsh2-data[,2])/sd)^2)]
        #for(i in 1:n) {
        #rndpars[i,chivals:=sum(((jvReal(c(i0,n,rs,rsh,rsh2,alpha),data,temp)[,2] - data[,2])/sd)^2)]
        #}
        
        rndpars[,chivals:=apply(rndpars,1,function(x) sum(((jvReal(x,data,temp)[,2] - data[,2])/sd)^2))]
        
        print(rndpars)
        
        # liczba stopni swobody
        df <- length(data[,1]) - length(sel) - 1
        
        # wyliczenie wartosci krytycznej kwantyla chi^2
        chi_crit <- qchisq(0.95, df)
        
        plot(rndpars[chivals<=chi_crit, .(i0,n)])
        
        # indeksy zaakceptowanych wartosci
        #indexes <- chi_values <= chi_crit
#         
#         # zrobienie wykresu zaakceptowanych wartosci
#         reppars_df_acc <- reppars_df[indexes,sel]
#         plot(reppars_df_acc)
#         
#         # punkty zaakceptowane w skali n~log(i0), ktore spelniaja warunek <= chi krytyczne
#         # logarytmowanie bo i0~exp(n)
#         logvalues_acc <- cbind(log(reppars_df_acc[,1]), reppars_df_acc[,2])
#         logvalues_acc <- as.data.frame(logvalues_acc)
#         names(logvalues_acc) <- names(reppars_df_acc[,sel])
#         
#         # tylko zaakceptowane parametry, zamienione do ramki danycyh z wartosciami numerycznymi
#         #reppars_df_acc <- as.data.frame(t(apply(reppars_df[indexes,], 1, as.numeric)))
#         
#         # stworzenie formuly do fitowania
#         formula <- as.formula(paste(rev(names(logvalues_acc)), collapse = "~"))
#         
#         # model dopasowania
#         model <- lm(formula, logvalues_acc)
#         intercept <- model$coefficients[[1]]
#         slope <- model$coefficients[[2]]
#         
#         # model danych zbudowany na podstawie dopasowania
#         modelled_data <- cbind(exp(logvalues_acc[,1]),logvalues_acc[,1]*slope + intercept)
#         
#         # posortowanie w celu lepszego rysowania
#         modelled_data <- modelled_data[order(modelled_data[,1]),]
#         
#         # dodanie linii do wykresu z punktami
#         lines(modelled_data)
#         
#         # stworznie ramki danych dla parametru
#         modelparams <- data.frame(intercept=intercept, slope=slope)
#         
#         # zwrocenie ramki danych
#         modelparams
        
}