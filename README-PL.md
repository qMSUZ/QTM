# Pakiet obliczeniowych dla kwantowej metody trajektorii / Quantum Trajectory Method Package

Pakiet Quantum Trajectory Method (QTM), to implementacja metody kwantowych trajektorii z wykorzystaniem protokołu MPI. Pozwala to na wykorzystanie MPI do równoległego przetwarzania co przekłada się na możliwość uzyskanie większej wydajności w przypadku większych problemów dla dużej liczby obliczanych trajektorii.

Pakiet QTM został opracowany w języku C++ z wykorzystaniem metod numerycznych rozwiązujących układy równań różniczkowych z pakietu ZVODE.


# Instalacja/ Installation

Instalacja pakiety ze źródeł wymaga obecności kompilatora jezyk C++ oraz Fortran. Niezbędne są takżę biblioteki Blas (i Lapack), a także wybranej implementacji protokou MPI, może to być np. pakiet MPICH lub OpenMPI.

QTM był rozwijany w ramach systemu Ubuntu (http://www.ubuntu.com), toteż dalszy opis jego  kompilacji będzie realizowany na przykładzie systemu Ubuntu.

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

Kompilację wybranego przykładu np. hamiltonianu trojliniowego, wykonamy jednym poleceniem np.:

make triham-ex

Uruchomiemie aplikacji nie wymaga posiadania  klastra obliczeniowego opartego o MPI. Można korzystać z pakietu QTM również w ramach lokalnej maszyny, wykorzystując dostępne rdzenie obliczeniowe. Uruchomienie jednak pozostaje typowe dla aplikacji MPI np.:

mpirun -n 5 triham-ex

Co oznacza uruchomienie pięciu procesów MPI, jeden z nich będzie procesem zarządzającym, natomiast cztery będą pełnić rolę węzłów obliczeniowych. 

# Przykład / Basic example

Pakiet QTM  nie oferuje pracy interktywnej i wymaga podania postaci hamiltonianu oraz operatorów collapsu oraz wartości oczekiwanej oraz  funkcji stosowanej przez pakiet ZVODE podczas rozwiązywania układu równań różniczkowych. Nie mniej, są to jedyne elementy jakie trzeba samodzielnie określić. Pozostałe detale implementacji nie są eksponowane dla użytkownika, ktory chce rozwiązać dany problem obliczeniowy za pomocą metody kwantowych trajektorii.

W pakiecie znadują się dodatkowe przykłady poniżej prezentujemy tylko przykład odnoszący się do tzw. hamiltonianu unitarnego. Niezbędne definicje podstawowych struktur danych to:

#include "complexnum.h"
#include "qtm.h" 

const size_t COLLAPSE_OPERATORS = 1;
const size_t Ntrj = 1;
const size_t N = 100;
const size_t WAVEVECTOR_LEAD_DIM = 2;
const size_t WAVEVECTOR_LEAD_DIM_SQR = 4;

uMatrix< simpleComplex<double>, WAVEVECTOR_LEAD_DIM > c_ops[ COLLAPSE_OPERATORS ];
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > collapse_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > expect_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > H;
simpleComplex<double> alpha[WAVEVECTOR_LEAD_DIM];

#include "qtm.cc"

extra_options opt;

Ponieważ QTM to szkielet, to należy na początku pliku źródłowego włączyć dwa pliki nagłówkowe, obsługujące podstawowe typy danych oraz samą metodę QTM, wymaga ona także dołączenia treści funkcji mpi_main, która znajduje się  w pliku qtm.cc.

W treści funkcji main, możemy zatem podać następujący kod (msc to skrót funkcji make_simpleComplex tworzącej liczbę zespoloną):

int r = 0;
co[0] = msc( 0.0, 0.0 );  co[1] = msc( 0.05, 0.0 );
co[2] = msc( 0.05, 0.0 ); co[3] = msc( 0.0, 0.0 );

eo[0] = msc( 1.0, 0.0); eo[1] = msc( 0.0, 0.0);
eo[2] = msc( 0.0, 0.0); eo[3] = msc(-1.0, 0.0);

alpha[0] = msc( 1.0, 0.0);
alpha[1] = msc( 0.0, 0.0);
	
H[0] = msc( -0.00125, 0.0); 
H[1] = msc( 0.0, -0.62831853);
H[2] = msc( 0.0, -0.62831853); 
H[3] = msc( -0.00125, 0.0);
	
c_ops[0].rows=2; c_ops[0].cols=2;
c_ops[0].m = co; 

opt.type_output = OUTPUT_FILE;
opt.only_final_trj = 1;
opt.ode_method = METADAMS;
opt.tolerance = 1e-7;
opt.file_name = strdup("output-data.txt");
opt.fnc = &myfex_fnc_f1;
	
r = mpi_main<N, Ntrj, WV_LD, WV_LD_SQR, 1>(argc, argv, 1, 0, 10, 1, 1, opt);
      
Po kompilacji pliku źródłowego np. w taki bezpośredni sposób:

g++ unitary-mc-ex.cc rgen_lfsr113.cc zvode.o zgesl.o zgefa.o zgbsl.o zgbfa.o -o unitary-ex -lblas -lgfortran -lopenmpi

Możemy urchomić przykład w środowsku klasatra MPI, nawet na maszynie lokalnej np.

mpirun -n 9 unitary-ex

Uruchamiamy 9 procesów MPI, z jednym węzłem głównym i ośmioma roboczymi.
