# Pakiet obliczeniowych dla kwantowej metody trajektorii / Quantum Trajectory Method Package

Pakiet Quantum Trajectory Method (QTM), to implementacja metody kwantowych trajektorii z wykorzystaniem protokołu MPI. Pozwala to na wykorzystanie MPI do równoległego przetwarzania co przekłada się na możliwość uzyskanie większej wydajności w przypadku wiekszych problemów dla dużej liczby obliczanych trajektorii.

Pakiet QTM został opracowany w języku C++ z wykorzystaniem metod numerycznych rozwiązujących układy równań różniczkowych z pakietu ZVODE.


# Instalacja/ Installation

Instalacja pakiety ze źródeł wymaga obecności kompilatora jezyk C++, Fortran bibliotek blas (i ewentualnie Lapack) a także wybranej implementacji protokou MPI, może to być np. pakiet MPICH lub OpenMPI.

QTM był rozwijany w ramach systemu Ubuntu, toteż dalszy opis jego  kompilacji będzie realizowany na przykładzie systemu Ubuntu.

Pierwszym krokiem jest instalacja podstawowych pakietów związanych z narzędziami dla programistów:

sudo apt-get install build-essential

Następnie można doinstalować biblioteki BLAS i LAPACK:

sudo apt-get install libblas-dev liblapack-dev

Należy także wgrać pakiet do obsługi MPI, wybierzemy odmianę OpenMPI:

sudo apt-get install libopenmpi-dev

Aktualną wersję pakietu QTM można ściągnąć bezpośrednio z repozytorium:

git clone https://github.com/qMSUZ/QTM

Możemy przejść do nowo utworzonego katalogu:

cd QTM

Właściwą kompilację QTM musimy poprzedzić ściągnięcem plików z pakiey ZVODE, nie są on bowiem umieszczone w repozytorium QTM:

make zvode_download

Kompilację wybranego przykładu np. hamiltonianu trojliniowego, wykonamy jednym poleceniem:

make triham

Uruchomiemie aplikacji nie wymaga posiadania  klastra obliczeniowego opartego o MPI. Można korzystać z pakietu QTM również loklanie, wykorzystując dostępne rdzenie obliczeniowe. Uruchomienie jednak pozostaje typowe dla aplikacji MPI np.:

mpirun -n 4 triham-ex

Dla aplikacji dotyczącej symulacji hamiltonianiu trójliniowego.

# Przykład / Basic example

Pakiet QTM wymaga podania postaci hamiltonianu oraz operatorów collapsu oraz wartości oczekiwanej oraz  funkcji stosowanej przez pakiet ZVODE podczas rozwiązywania układu równań różniczkowych.

W pakiecie znadują się dodatkowe przykłady poniżej prezentujemy tylko przykład odnoszący się do tzw. hamiltonianu unitarnego.

