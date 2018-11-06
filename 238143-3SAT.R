# Ładowanie biblioteki do algorytmów genetycznych
library(genalg)

# Funkcje pomocnicze
## Funkcja, która zamienia bity
swap <- function(a) {
    ifelse(a==0, return(1), return(0))
}

## Funkcja, która ładuje listę klauzul z pliku do wektora 
load_clauses_to_a_vector <- function(file_location) {
    connection <- file(file_location) 
    len <- length(readLines(connection)) - 11 # eliminacja 8 linii z góry, 3 z dołu
    data <- read.table(connection , skip=8, nrows=len)[,1:3] # kolumny od 1 do 3
    vector <- as.vector(unname(unlist(t(data)))) # niejawna konwersja data.frame do macierzy przez t() oraz konwersja macierzy do wektora
    return(vector)
}

## Zmienna a'la globalna w C
clauses_vectorized <- NULL

## Ustalenie ścieżki do katalogu z klauzulami
path <- "./samples/"

## Funkcja fitness
fitness3SAT <- function(chromosome) {
    score <- 0 
    for(i in seq(1, length(clauses_vectorized), by = 3)) { # pętla przechodząca co 3 elementy (klauzulę) w wektorze klauzul
        tmp <- c() # tymczasowy wektor dla wartości logicznych
        for(j in 0:2) { # dla każdego elementu spośród tych 3 wcześniej wybranych, czyli dla każdego elementu z klauzuli
            if (clauses_vectorized[i+j] < 0) { # jeśli element wskazuje na miejsce w chromosomie liczbą ujemną, to oznacza negację, więc trzeba 
                tmp <- c(tmp,swap(chromosome[abs(clauses_vectorized[i+j])]))    # znaleźć wartość bezwzględną l.p. miejsca, sprawdzić, co się 
                                                                                # kryje w chromosomie na tym miejscu i zanegować, czyli zeswapować
            } else {
                tmp <- c(tmp,chromosome[abs(clauses_vectorized[i+j])])  # "zwyczajne" (w porównaniu do poprzedniego) sprawdzenie wartości 
                                                                        # w chromosomie
            } 
        } # tym sposobem w wktorze tmp mamy wartości logiczne, które po podstawieniu do klauzuli czynią ją spełnialną
        
        if(1 %in% tmp) { # jeśli w wektorze wartości logicznych jest 1, klauzula jest spełniona, więc
            score <- score + 1 # należy dodać punkt dla tej klauzuli
        }
    } # i przejść do następnej 3 
    return(-score)
}


# CZĘŚĆ PIERWSZA - badanie instancji problemu - tu 50 
## Wektory pomocnicze
v_seq <- c()
v_spopSize <- c()
v_found_in_iteration <- c()
v_how_many_found <- c()
v_time <- c()
v_is_elitism <- c()

## Analiza
clauses_vectorized <- load_clauses_to_a_vector(paste(path, "/50/2.cnf", sep=""))
for(i in seq(0.00, 0.5, by=0.05)) {
    for(j in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) {
        for(k in c(T, F)) {
            v_seq <- c(v_seq, i)
            v_spopSize <- c(v_spopSize, j)
            
            algorithm <- NULL
            algorithm_running_time <- system.time(algorithm <- rbga.bin(size = length(unique(abs(clauses_vectorized))), popSize = j, iters = 100, 
                                                                        mutationChance = i, elitism = k, evalFunc = fitness3SAT))[1]
            
            v_found_in_iteration <- c(v_found_in_iteration, which.min(algorithm$best)) #https://stat.ethz.ch/pipermail/r-help/2007-January/123788.html
            v_how_many_found <- c(v_how_many_found, max(abs(algorithm$best)))
            v_time <- c(v_time, algorithm_running_time)
            v_is_elitism <- c(v_is_elitism, k)
        }
    }
}

comparison_frame <- data.frame(v_seq, v_spopSize, v_found_in_iteration, v_how_many_found, v_time, v_is_elitism)

## Dla jakich parametrów algorytm działa najszybciej?
fastest_order <- comparison_frame[order(v_time),]
write.csv(fastest_order,'fastest_order.csv')

## Dla jakich parametrów algorytm działa najefektywniej?
most_effective_order <- comparison_frame[order(-v_how_many_found),]
write.csv(most_effective_order,'most_effective_order.csv')


# CZĘŚĆ DRUGA
clause_sizes <- c(10,20,30,40,50)
times <- c()

for(i in clause_sizes) {
    running_time_summarized <- 0
    for(j in c(1,2,3)) {
        clauses_vectorized <- load_clauses_to_a_vector(paste(path, paste(i, paste("/", paste(j, ".cnf", sep=""), sep=""), sep=""), sep=""))
        algorithm_running_time <- system.time(rbga.bin(size = length(unique(abs(clauses_vectorized))), popSize = 200, iters = 100,
                                    mutationChance = 0.05, elitism = T, evalFunc = fitness3SAT))[1]
        running_time_summarized <- running_time_summarized + algorithm_running_time # "system.time" jest wektorem[?], który zawiera w sobie 3 wyniki, 1 pozycja to v_time uzytkownika
    }
    average_time <- running_time_summarized / 3.0 # dzielenie przez liczbę iteracji (przykładów dla danego rozmiaru formuły)
    times <- c(times, average_time)
}

## Tabela
times_frame <- data.frame(clause_sizes, times)
write.csv(times_frame,'times_frame.csv')

## Wykres
plot(clause_sizes, times, type="o", col="red", xlab="Ilość klauzul", ylab="Średni czas wykonania na podstawie 3 próbek", main="Czas wykonywania alg. genetycznego w zależności od długości formuły")
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
