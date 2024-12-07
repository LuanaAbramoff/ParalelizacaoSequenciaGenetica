#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cctype> 

#define NUM_BASES 4 // A, T, C, G

using namespace std;

void contarBases(const string &sequencia, int *contagemLocal) {
    #pragma omp parallel for reduction(+:contagemLocal[:NUM_BASES])
    for (size_t i = 0; i < sequencia.size(); i++) {
        char base = toupper(sequencia[i]); // Converter para maiúsculo
        switch (base) {
            case 'A': contagemLocal[0]++; break;
            case 'T': contagemLocal[1]++; break;
            case 'C': contagemLocal[2]++; break;
            case 'G': contagemLocal[3]++; break;
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    std::vector<char> sequencia;
    int tamPorProcesso;
    int contagemLocal[NUM_BASES] = {0}; // A, T, C, G
    int contagemGlobal[NUM_BASES] = {0};

    if (rank == 0) {
        // Leitura do arquivo FASTA
        std::ifstream arquivo("chr1.subst.fa");
        if (!arquivo.is_open()) {
            std::cerr << "Erro ao abrir o arquivo genoma.fa" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::string linha;
        while (std::getline(arquivo, linha)) {
            if (linha[0] != '>') { // Ignorar linha de descrição
                for (char base : linha) {
                    if (base != 'N' && base != 'n') { // Ignorar caracteres 'N'
                        sequencia.push_back(base);
                    }
                }
            }
        }
        arquivo.close();

        tamPorProcesso = sequencia.size() / numProcs;
    }

    // Manda o tamanho para todos os processos
    MPI_Bcast(&tamPorProcesso, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Distribuição dos dados para cada processo
    vector<char> parteSequencia(tamPorProcesso);
    MPI_Scatter(sequencia.data(), tamPorProcesso, MPI_CHAR, parteSequencia.data(), tamPorProcesso, MPI_CHAR, 0, MPI_COMM_WORLD);

    contarBases(parteSequencia.data(), contagemLocal);

    // Redução dos resultados
    MPI_Reduce(contagemLocal, contagemGlobal, NUM_BASES, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Contagem total de bases:" << std::endl;
        std::cout << "A: " << contagemGlobal[0] << std::endl;
        std::cout << "T: " << contagemGlobal[1] << std::endl;
        std::cout << "C: " << contagemGlobal[2] << std::endl;
        std::cout << "G: " << contagemGlobal[3] << std::endl;
    }

    MPI_Finalize();
    return 0;
}
