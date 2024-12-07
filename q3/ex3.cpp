#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cctype> 

using namespace std;

// Função para contar o número de códons "AUG"
int contarCodonAUG(const std::vector<char> &sequencia) {
    int count = 0;
    #pragma omp parallel for reduction(+:count)
    for (size_t i = 0; i < sequencia.size() - 2; i++) { // Verificar limite para evitar out-of-bound
        if (sequencia[i] == 'A' && sequencia[i + 1] == 'U' && sequencia[i + 2] == 'G') {
            count++;
        }
    }
    return count;
}

int main(int argc, char *argv[]) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    vector<char> sequencia;
    int tamPorProcesso;

    // Tempo de execução total
    double tempoInicioTotal = MPI_Wtime();

    double tempoInicioLeitura = 0.0, tempoFimLeitura = 0.0;
    if (rank == 0) {
        // Leitura do arquivo FASTA
        tempoInicioLeitura = MPI_Wtime();
        std::ifstream arquivo("rna_transcrito.txt");
        if (!arquivo.is_open()) {
            std::cerr << "Erro ao abrir o arquivo rna_transcrito.txt" << std::endl;
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
        tempoFimLeitura = MPI_Wtime();

        tamPorProcesso = sequencia.size() / numProcs;
    }

    // Manda o tamanho para todos os processos
    MPI_Bcast(&tamPorProcesso, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Distribuição dos dados para cada processo
    vector<char> parteSequencia(tamPorProcesso);
    MPI_Scatter(sequencia.data(), tamPorProcesso, MPI_CHAR, parteSequencia.data(), tamPorProcesso, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Tempo para contar os códons
    double tempoInicioContagem = MPI_Wtime();
    int contagemLocal = contarCodonAUG(parteSequencia);
    double tempoFimContagem = MPI_Wtime();

    // Redução dos resultados
    int contagemTotal = 0;
    MPI_Reduce(&contagemLocal, &contagemTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    double tempoFimTotal = MPI_Wtime();
    double tempoTotal = tempoFimTotal - tempoInicioTotal;

    // Processo mestre exibe os resultados e tempos
    if (rank == 0) {
        cout << "Número total de proteínas inicializadas (AUG): " << contagemTotal << endl;
        cout << "Tempo de leitura do arquivo: " << (tempoFimLeitura - tempoInicioLeitura) << " segundos" << endl;
        cout << "Tempo de contagem dos códons: " << (tempoFimContagem - tempoInicioContagem) << " segundos" << endl;
        cout << "Tempo total de execução: " << tempoTotal << " segundos" << endl;
    }

    MPI_Finalize();
    return 0;
}
