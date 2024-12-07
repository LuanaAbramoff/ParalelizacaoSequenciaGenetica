#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <cctype>
#include <climits>

using namespace std;

// Função para transcrever DNA para RNA
void transcreveDNAparaRNA(vector<char> &sequencia) {
    #pragma omp parallel for
    for (size_t i = 0; i < sequencia.size(); i++) {
        if (sequencia[i] == 'T') {
            sequencia[i] = 'U';
        }
    }
}

// Função para ler um arquivo FASTA e retornar a sequência
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
                    sequencia += toupper(base);
                }
            }
        }
    }
    arquivo.close();
    return sequencia;
}

// Função para salvar uma sequência em um arquivo
void salvarSequencia(const string &nomeArquivo, const string &sequencia) {
    ofstream arquivo(nomeArquivo, ios::app); // Adiciona ao arquivo existente
    if (arquivo.is_open()) {
        arquivo << sequencia << endl;
    } else {
        cerr << "Erro ao salvar no arquivo " << nomeArquivo << endl;
    }
    arquivo.close();
}

int main(int argc, char *argv[]) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // Lista de arquivos disponíveis
    vector<string> arquivos;
    for (int i = 1; i <= 22; i++) {
        arquivos.push_back("arquivos_fasta/chr" + to_string(i) + ".subst.fa");
    }

    // Divisão de arquivos entre os processos
    vector<string> arquivosLocais;
    for (size_t i = rank; i < arquivos.size(); i += numProcs) {
        arquivosLocais.push_back(arquivos[i]);
    }

    double tempoInicioTotal = MPI_Wtime();
    double tempoLeitura = 0.0, tempoTranscricao = 0.0;

    // Transcrição de cada arquivo localmente
    string rnaTranscritoTotal;
    for (const auto &arquivo : arquivosLocais) {
        double inicioLeitura = MPI_Wtime();
        string sequencia = lerArquivoFASTA(arquivo);
        double fimLeitura = MPI_Wtime();
        tempoLeitura += fimLeitura - inicioLeitura;

        vector<char> sequenciaVec(sequencia.begin(), sequencia.end());

        double inicioTranscricao = MPI_Wtime();
        transcreveDNAparaRNA(sequenciaVec);
        double fimTranscricao = MPI_Wtime();
        tempoTranscricao += fimTranscricao - inicioTranscricao;

        rnaTranscritoTotal += string(sequenciaVec.begin(), sequenciaVec.end());
    }

    // Preparar para MPI_Gatherv
    long long localSize = rnaTranscritoTotal.size();
    vector<long long> tamanhos(numProcs, 0);

    MPI_Gather(&localSize, 1, MPI_LONG_LONG, tamanhos.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    vector<int> tamanhosInt(numProcs, 0); // Vetor de tamanhos para MPI_Gatherv
    vector<int> deslocamentos(numProcs, 0);
    long long totalSize = 0;

    if (rank == 0) {
        for (int i = 0; i < numProcs; i++) {
            if (tamanhos[i] > INT_MAX) {
                cerr << "Erro: Tamanho de segmento excede o limite de int!" << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            tamanhosInt[i] = static_cast<int>(tamanhos[i]);
            deslocamentos[i] = static_cast<int>(totalSize);
            totalSize += tamanhos[i];
        }
    }

    string rnaFinal;
    if (rank == 0) {
        rnaFinal.resize(totalSize);
    }

    MPI_Gatherv(
        rnaTranscritoTotal.data(), static_cast<int>(localSize), MPI_CHAR,
        rank == 0 ? rnaFinal.data() : nullptr, tamanhosInt.data(),
        deslocamentos.data(), MPI_CHAR, 0, MPI_COMM_WORLD);

    double tempoFimTotal = MPI_Wtime();
    double tempoTotal = tempoFimTotal - tempoInicioTotal;

    // O processo mestre salva o resultado no arquivo final
    if (rank == 0) {
        salvarSequencia("rna_transcrito.txt", rnaFinal);
        cout << "Transcrição concluída e salva em rna_transcrito.txt" << endl;
        cout << "Tempo total de leitura dos arquivos: " << tempoLeitura << " segundos" << endl;
        cout << "Tempo total de transcrição: " << tempoTranscricao << " segundos" << endl;
        cout << "Tempo total de execução: " << tempoTotal << " segundos" << endl;
    }

    MPI_Finalize();
    return 0;
}
