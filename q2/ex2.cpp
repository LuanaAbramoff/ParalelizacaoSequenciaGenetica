#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cctype> 

using namespace std;

void transcreveDNAparaRNA(vector<char> &sequencia) {
    #pragma omp parallel for
    for (size_t i = 0; i < sequencia.size(); i++) {
        if(sequencia[i] == 'T') {
            sequencia[i] = 'U';
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    vector<char> sequencia;
    int tamPorProcesso;

    // Variáveis para tempo
    double tempoInicioTotal = MPI_Wtime();
    double tempoInicioLeitura = 0.0, tempoFimLeitura = 0.0;
    double tempoInicioDistribuicao = 0.0, tempoFimDistribuicao = 0.0;
    double tempoInicioTranscricao = 0.0, tempoFimTranscricao = 0.0;
    
    if (rank == 0) {
        // Leitura do arquivo FASTA
        tempoInicioLeitura = MPI_Wtime();
        std::ifstream arquivo("arquivos_fasta/chr2.subst.fa");
        if (!arquivo.is_open()) {
            std::cerr << "Erro ao abrir o arquivo genoma.fa" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::string linha;
        while (std::getline(arquivo, linha)) {
            if (linha[0] != '>') { // Ignorar linha de descrição
                for (char base : linha) {
                    if (base != 'N' && base != 'n') { // Ignorar caracteres 'N'
                        sequencia.push_back(toupper(base));
                    }
                }
            }
        }
        arquivo.close();

        tamPorProcesso = sequencia.size() / numProcs;
        tempoFimLeitura = MPI_Wtime();
    }

    // Manda o tamanho para todos os processos
    MPI_Bcast(&tamPorProcesso, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Distribuição dos dados para cada processo
    tempoInicioDistribuicao = MPI_Wtime();
    vector<char> parteSequencia(tamPorProcesso);
    MPI_Scatter(sequencia.data(), tamPorProcesso, MPI_CHAR, parteSequencia.data(), tamPorProcesso, MPI_CHAR, 0, MPI_COMM_WORLD);
    tempoFimDistribuicao = MPI_Wtime();

    // Transcrição de DNA para RNA
    tempoInicioTranscricao = MPI_Wtime();
    transcreveDNAparaRNA(parteSequencia);
    tempoFimTranscricao = MPI_Wtime();

    // Redução dos resultados
    MPI_Gather(parteSequencia.data(), tamPorProcesso, MPI_CHAR, sequencia.data(), tamPorProcesso, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        ofstream arquivoSaida("rna_transcrito.txt");
        if(arquivoSaida.is_open()) {
            cout << sequencia.size() << endl;
            for (char base : sequencia) {
                arquivoSaida << base;
            }
        }
        arquivoSaida.close();
        cout << "Transcrição concluída e salva em rna_transcrito.txt" << endl;
        
        // Exibindo os tempos
        double tempoTotal = MPI_Wtime() - tempoInicioTotal;
        cout << "Tempo de leitura do arquivo: " << (tempoFimLeitura - tempoInicioLeitura) << " segundos" << endl;
        cout << "Tempo de distribuição dos dados: " << (tempoFimDistribuicao - tempoInicioDistribuicao) << " segundos" << endl;
        cout << "Tempo de transcrição de DNA para RNA: " << (tempoFimTranscricao - tempoInicioTranscricao) << " segundos" << endl;
        cout << "Tempo total de execução: " << tempoTotal << " segundos" << endl;
    }

    MPI_Finalize();
    return 0;
}
