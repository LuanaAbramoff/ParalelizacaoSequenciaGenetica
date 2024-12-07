#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <mpi.h>
#include <omp.h>

// Função para inicializar o mapa de tradução de códons para números (aminoácidos)
std::map<std::string, int> criarMapaCodons() {
    std::map<std::string, int> mapa;
    mapa["AUG"] = 1; // Metionina (início)
    mapa["CCA"] = 2; mapa["CCG"] = 2; mapa["CCU"] = 2; mapa["CCC"] = 2; // Prolina
    mapa["UCU"] = 3; mapa["UCA"] = 3; mapa["UCG"] = 3; mapa["UCC"] = 3; // Serina
    mapa["CAG"] = 4; mapa["CAA"] = 4; // Glutamina
    mapa["ACA"] = 5; mapa["ACC"] = 5; mapa["ACU"] = 5; mapa["ACG"] = 5; // Treonina
    mapa["UGC"] = 6; mapa["UGU"] = 6; // Cisteína
    mapa["GUG"] = 7; mapa["GUA"] = 7; mapa["GUC"] = 7; mapa["GUU"] = 7; // Valina
    mapa["UGA"] = 0; // Códons de parada
    return mapa;
}

// Função para ler o arquivo RNA e retornar a sequência como um vetor de caracteres
std::vector<char> lerArquivoRNA(const std::string &nomeArquivo) {
    std::ifstream arquivo(nomeArquivo);
    std::vector<char> sequencia;
    if (!arquivo.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << nomeArquivo << std::endl;
        return sequencia;
    }

    std::string linha;
    while (std::getline(arquivo, linha)) {
        for (char base : linha) {
            sequencia.push_back(base);
        }
    }
    arquivo.close();
    return sequencia;
}

std::vector<int> traduzirRNA(const std::vector<char> &sequencia, const std::map<std::string, int> &mapaCodons) {
    std::vector<int> aminoacidos;

    #pragma omp parallel
    {
        std::vector<int> localAminoacidos;
        bool localTraduzindo = false;

        #pragma omp for
        for (size_t i = 0; i < sequencia.size() - 2; i += 3) {
            std::string codon = {sequencia[i], sequencia[i + 1], sequencia[i + 2]};
            if (mapaCodons.count(codon)) {
                int aminoacido = mapaCodons.at(codon);
                if (!localTraduzindo && aminoacido == 1) { // Códon de início
                    localTraduzindo = true;
                }
                if (localTraduzindo) {
                    if (aminoacido == 0) { // Códon de parada
                        localTraduzindo = false; // Parar de traduzir nesta thread
                        continue;
                    }
                    localAminoacidos.push_back(aminoacido);
                }
            }
        }

        // Região crítica para combinar resultados
        #pragma omp critical
        aminoacidos.insert(aminoacidos.end(), localAminoacidos.begin(), localAminoacidos.end());
    }

    return aminoacidos;
}

int main(int argc, char *argv[]) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    std::vector<char> sequenciaTotal;
    std::map<std::string, int> mapaCodons = criarMapaCodons();

    double tempoInicioTotal = MPI_Wtime();
    double tempoLeitura = 0, tempoDistribuicao = 0, tempoTraducao = 0, tempoFimTotal = 0;

    if (rank == 0) {
        double tempoInicioLeitura = MPI_Wtime();
        if (argc < 2) {
            std::cerr << "Uso: mpirun -np <num_procs> ./traducao_proteica <arquivo_rna>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        sequenciaTotal = lerArquivoRNA(argv[1]);
        if (sequenciaTotal.empty()) {
            std::cerr << "Erro: Arquivo vazio ou não pôde ser lido." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        double tempoFimLeitura = MPI_Wtime();
        tempoLeitura = tempoFimLeitura - tempoInicioLeitura;
    }

    long long tamanhoTotal = sequenciaTotal.size();
    MPI_Bcast(&tamanhoTotal, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (tamanhoTotal == 0) {
        if (rank == 0) std::cerr << "Erro: Sequência de RNA está vazia." << std::endl;
        MPI_Finalize();
        return 1;
    }

    std::vector<int> sendCounts(numProcs, tamanhoTotal / numProcs);
    std::vector<int> displs(numProcs, 0);
    int resto = tamanhoTotal % numProcs;

    for (int i = 0; i < resto; ++i) {
        sendCounts[i]++;
    }
    for (int i = 1; i < numProcs; ++i) {
        displs[i] = displs[i - 1] + sendCounts[i - 1];
    }

    std::vector<char> sequenciaLocal(sendCounts[rank]);
    double tempoInicioDistribuicao = MPI_Wtime();
    MPI_Scatterv(sequenciaTotal.data(), sendCounts.data(), displs.data(), MPI_CHAR,
                sequenciaLocal.data(), sendCounts[rank], MPI_CHAR, 0, MPI_COMM_WORLD);
    double tempoFimDistribuicao = MPI_Wtime();
    tempoDistribuicao = tempoFimDistribuicao - tempoInicioDistribuicao;

    double tempoInicioTraducao = MPI_Wtime();
    std::vector<int> aminoacidosLocais = traduzirRNA(sequenciaLocal, mapaCodons);
    double tempoFimTraducao = MPI_Wtime();
    tempoTraducao = tempoFimTraducao - tempoInicioTraducao;

    int localSize = aminoacidosLocais.size();
    std::vector<int> recvCounts(numProcs, 0);
    MPI_Gather(&localSize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> displsAminoacidos(numProcs, 0);
    if (rank == 0) {
        for (int i = 1; i < numProcs; ++i) {
            displsAminoacidos[i] = displsAminoacidos[i - 1] + recvCounts[i - 1];
        }
    }

    std::vector<int> aminoacidosFinais;
    if (rank == 0) {
        int totalSize = displsAminoacidos[numProcs - 1] + recvCounts[numProcs - 1];
        aminoacidosFinais.resize(totalSize);
    }
    MPI_Gatherv(aminoacidosLocais.data(), localSize, MPI_INT, aminoacidosFinais.data(),
                recvCounts.data(), displsAminoacidos.data(), MPI_INT, 0, MPI_COMM_WORLD);

    tempoFimTotal = MPI_Wtime();

    if (rank == 0) {
        std::ofstream arquivoSaida("aminoacidos.txt");
        for (int aminoacido : aminoacidosFinais) {
            if (aminoacido != 0) {
                arquivoSaida << aminoacido << " ";
            }
        }
        std::cout << "Tempo de Leitura: " << tempoLeitura << " segundos" << std::endl;
        std::cout << "Tempo de Distribuição: " << tempoDistribuicao << " segundos" << std::endl;
        std::cout << "Tempo de Tradução: " << tempoTraducao << " segundos" << std::endl;
        std::cout << "Tempo Total: " << (tempoFimTotal - tempoInicioTotal) << " segundos" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
