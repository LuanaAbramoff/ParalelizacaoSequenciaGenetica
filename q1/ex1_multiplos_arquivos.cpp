#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <cctype>

#define NUM_BASES 4 // A, T, C, G

using namespace std;

// Função para contar as bases de uma sequência
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

// Função para processar um arquivo e retornar a sequência
string lerArquivoFASTA(const string &nomeArquivo) {
    ifstream arquivo(nomeArquivo);
    if (!arquivo.is_open()) {
        cerr << "Erro ao abrir o arquivo " << nomeArquivo << endl;
        return "";
    }

    string linha, sequencia;
    while (getline(arquivo, linha)) {
        if (linha[0] != '>') { // Ignorar linha de descrição
            for (char base : linha) {
                if (base != 'N' && base != 'n') { // Ignorar caracteres 'N'
                    sequencia += base;
                }
            }
        }
    }
    arquivo.close();
    return sequencia;
}

int main(int argc, char *argv[]) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // Vetor de contagem local e global
    int contagemLocal[NUM_BASES] = {0};
    int contagemGlobal[NUM_BASES] = {0};

    // Lista de arquivos disponíveis
    vector<string> arquivos;
    for (int i = 1; i <= 22; i++) {
        arquivos.push_back("arquivos_fasta/chr" + to_string(i) + ".subst.fa");
    }

    // Divisão dos arquivos entre os processos
    vector<string> arquivosLocais;
    for (size_t i = rank; i < arquivos.size(); i += numProcs) {
        arquivosLocais.push_back(arquivos[i]);
    }

    double tempo_inicio = MPI_Wtime();

    // Processar os arquivos atribuídos a este processo
    for (const auto &arquivo : arquivosLocais) {
        string sequencia = lerArquivoFASTA(arquivo);
        contarBases(sequencia, contagemLocal);
    }

    double tempo_final = MPI_Wtime();
    double tempo_local = tempo_final - tempo_inicio;
    double tempo_maximo;

    // Reduzir os tempos para obter o tempo máximo
    MPI_Reduce(&tempo_local, &tempo_maximo, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Redução das contagens
    MPI_Reduce(contagemLocal, contagemGlobal, NUM_BASES, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Contagem total de bases:" << endl;
        cout << "A: " << contagemGlobal[0] << endl;
        cout << "T: " << contagemGlobal[1] << endl;
        cout << "C: " << contagemGlobal[2] << endl;
        cout << "G: " << contagemGlobal[3] << endl;

        cout << "O tempo de execução do programa foi de " << tempo_maximo << " segundos." << endl;
    }

    MPI_Finalize();
    return 0;
}
