#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define ignite_prob 0  //A tree ignites with this probability (f), 0 <= ignite_prob <= 100
#define grow_prob 0    //A tree grows on an empty cell with this probability (p), 0 <= grow_prob <= 100
#define gen_prob 80     //Probability that a tree grows, used for the matrix generation, 0 <= gen_prob <= 100
#define S 25             //Number of steps

//MACRO per grafica
#define TREE "🌲"
#define BURN "🔥"
#define EMPTY "❌"

//inizializza una matrice m x n con alberi o celle vuote
void matrix_init(char *matrix, int m, int n){
    int random;
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            random = rand() % 101;
            if (random < 0) 
                random *= -1;

            if (random <= gen_prob){
                matrix[(i*n) + j] = 'T';
            }else{
                matrix[(i*n) + j] = 'E';
            }

            if(i == j){
                matrix[(i*n) + j] = 'B';
            }
        }
    }
}

//stampa una matrice m x n
void print_matrix(char *matrix, int m, int n){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            printf("[%c]", matrix[(i*n) + j]);
        }
        printf("\n");
    }
    printf("\n");
}

//stampa una matrice m x n utilizzando le emoji delle macro
void print_graphic_matrix(char *matrix, int m, int n, FILE *file){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            if(matrix[(i*n) + j] == 'T'){
                fprintf(file, "[%s]", TREE);
            } else if (matrix[(i*n) + j] == 'B'){
                fprintf(file, "[%s]", BURN);
            } else if (matrix[(i*n) + j] == 'E'){
                fprintf(file, "[%s]", EMPTY);
            }
            fprintf(file, " ");
            
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

//controlla i 4 vicini della cella recv[i][j] e se la cella già sta bruciando oppure dev'essere bruciata in base a f
void check_neighbors_ignite(char *recv, char *next, int i, int j, int n, int nrows){
    if (((i*n) + j - n) >= 0 && recv[((i*n) + j - n)] == 'B')    
        next[(i*n) + j] = 'B';
    if (((i*n) + j + n) < nrows*n && recv[((i*n) + j + n)] == 'B')    
        next[(i*n) + j] = 'B';
    if (((i*n) + j - 1)%n != n-1 && recv[((i*n) + j - 1)] == 'B')    
        next[(i*n) + j] = 'B';
    if (((i*n) + j + 1)%n != 0 && recv[((i*n) + j + 1)] == 'B')    
        next[(i*n) + j] = 'B';

    int random = rand() % 101;
    if (random < 0) 
        random *= -1;
    if (next[(i*n) + j] == 'B' || (random <= ignite_prob && random > 0))
        next[(i*n) + j] = 'B';
    else
        next[(i*n) + j] = 'T';
}

//controlla la probabilità p di crescita di un albero in una cella vuota, n è il numero di colonne
void check_tree_grows(char *next, int i, int j, int n){
    int random = rand() % 100;
    if (random < 0) 
        random *= -1;

    if (random <= grow_prob && random > 0)
        next[(i*n) + j] = 'T';
    else
        next[(i*n) + j] = 'E';
}

//creazione griglia cartesiana per MPI, con nproc righe e 1 colonna
void create_topology(int nproc, MPI_Comm *new_comm){
    int dim = 2, ndim[2] = {nproc, 1}, period[2] = {0, 0}, reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, dim, ndim, period, reorder, new_comm);
}

int main(int argc, char **argv){
    int rank, P;
    MPI_Comm comm_topology;
    MPI_Status status;
    MPI_Request request[2];
    FILE *file;

    //variabili per sottomatrici
    char *matrix, *recv, *next, *temp;
    int start_index, work_rows, nrows, R, k, M = atoi(argv[1]), N = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    //creazione griglia cartesiana processi
    create_topology(P, &comm_topology);
    MPI_Comm_rank(comm_topology, &rank);

    //variabili controllo matrice vuota
    int empty_counter = 0, empty_matrix = 0; 

    //variabili gatherv
    int gat_counts[P], gat_displ[P];

    //partizionamento della matrice e allocazione sottomatrici
    R =  M % P;
    //calcolo numero di righe di ogni sottomatrice (nrows), primo ed ultimo processo ricevono una riga di confine in meno
    if (rank == 0 || rank == P-1){
        nrows = M/P + 1;
    } else {
        nrows = M/P + 2;
    }
    // i primi R processi ricevono una riga in più, quando N non è multiplo di P
    if (rank < R){
        nrows += 1;
    }
    recv = (char *) malloc(sizeof(char) * nrows * N);
    next = (char *) malloc(sizeof(char) * nrows * N);

    //assegnazione indici riga iniziale e finale delle sottomatrici di lavoro (righe di scrittura work_rows)
    if (rank == 0) {
        start_index = 0;
        work_rows = nrows - 1; 
    } else if(rank == P-1) {
        start_index = 1;
        work_rows = nrows - 1;
    } else {
        start_index = 1;
        work_rows = nrows - 2;
    }

    srand(rank*rank);
    if (rank == 0) {
        file = fopen("correttezza_par.txt", "w");
        //ricalcolo per ogni processo con quante righe deve lavorare, necessario per la scatterv
        //non conviene comunicare i valori precedentemente computati (solo P elementi)
        int a_nrows[P], a_work_rows[P];
        for(int i=0; i<P; i++){
            if (i == 0 || i == P-1){
                a_nrows[i] = M/P + 1;
            } else {
                a_nrows[i] = M/P + 2;
            }
            if (i < R){
                a_nrows[i] += 1;
            }

            if (i == 0 || i == P-1) {
                a_work_rows[i] = a_nrows[i] - 1; 
            } else {
                a_work_rows[i] = a_nrows[i] - 2;
            }
        }

        //generazione matrice iniziale
        matrix = (char *) malloc(sizeof(char) * M * N);
        
        matrix_init(matrix, M, N);
        fprintf(file, "\nMatrice iniziale\n");
        print_graphic_matrix(matrix, M, N, file);

        //generazione array per scatterv e gatherv
        int counts[P];
        for(int i=0; i<P; i++){
            counts[i] = a_nrows[i] * N;
            gat_counts[i] = a_work_rows[i] * N;
        }

        int displacements[P];
        displacements[0] = 0;
        gat_displ[0] = 0;
        for(int i=1; i<P; i++){
            //posizionamento sulla riga precedente (di confine) a quella di lavoro del processo
            displacements[i] = displacements[i-1] + counts[i-1] - 2*N;
            gat_displ[i] = gat_displ[i-1] + gat_counts[i-1];
        }

        MPI_Scatterv(matrix, counts, displacements, MPI_CHAR, recv, nrows*N, MPI_CHAR, 0, comm_topology);
    } else {
        MPI_Scatterv(NULL, NULL, NULL, NULL, recv, nrows*N, MPI_CHAR, 0, comm_topology);
    }

    //simulazione su S passi discreti
    for(k=0; k<S; k++){
        empty_counter = 0;
        //comunicazione con altri processi
        if (rank == 0){
            //invio ultima riga al processo successivo
            MPI_Isend(&recv[work_rows*N - N], N, MPI_CHAR, rank+1, 0, comm_topology, &request[0]);
            //ricezione prima riga dal processo successivo
            MPI_Irecv(&recv[work_rows*N], N, MPI_CHAR, rank+1, 1, comm_topology, &request[1]);
        } else if (rank == P-1) {
            //invio prima riga al processo precedente
            MPI_Isend(&recv[N], N, MPI_CHAR, rank-1, 1, comm_topology, &request[1]);
            //ricezione ultima riga dal processo precedente
            MPI_Irecv(recv, N, MPI_CHAR, rank-1, 0, comm_topology, &request[0]);
        } else {
            //invio prima riga al precedente e ultima riga al successivo
            MPI_Isend(&recv[N], N, MPI_CHAR, rank-1, 1, comm_topology, &request[1]);
            MPI_Isend(&recv[work_rows*N], N, MPI_CHAR, rank+1, 0, comm_topology, &request[0]);
            //ricezione prima riga dal precedente e ultima riga dal successivo
            MPI_Irecv(recv, N, MPI_CHAR, rank-1, 0, comm_topology, &request[0]);
            MPI_Irecv(&recv[work_rows*N + N], N, MPI_CHAR, rank+1, 1, comm_topology, &request[1]);
        }

        //lavoro su righe non di confine, non dipendono da altri processi
        for(int i=start_index + 1; i<start_index + work_rows - 1; i++){
            for(int j=0; j<N; j++){
                if (recv[(i*N) + j] == 'E'){
                    check_tree_grows(next, i, j, N);
                } else if (recv[(i*N) + j] == 'T') {
                    check_neighbors_ignite(recv, next, i, j, N, nrows);
                } else if (recv[(i*N) + j] == 'B') {
                    next[(i*N) + j] = 'E';
                }
                
                if (next[(i*N) + j] == 'E')  empty_counter += 1;
            }
        }

        MPI_Waitall(2, request, MPI_STATUSES_IGNORE);

        //lavora esclusivamente sulle righe di bordo ricevute dagli altri processi
        int i = start_index;
        for(int j=0; j<N; j++){
            if (recv[(i*N) + j] == 'E'){
                check_tree_grows(next, i, j, N);
            } else if (recv[(i*N) + j] == 'T') {
                check_neighbors_ignite(recv, next, i, j, N, nrows);
            } else if (recv[(i*N) + j] == 'B') {
                next[(i*N) + j] = 'E';
            }

            if (next[(i*N) + j] == 'E')  empty_counter += 1;
        }
        
        i = start_index + work_rows - 1;
        if(start_index != start_index + work_rows - 1){//quando un processo ha solo una riga da calcolare, N = P
            for(int j=0; j<N; j++){
                if (recv[(i*N) + j] == 'E'){
                    check_tree_grows(next, i, j, N);
                } else if (recv[(i*N) + j] == 'T') {
                    check_neighbors_ignite(recv, next, i, j, N, nrows);
                } else if (recv[(i*N) + j] == 'B') {
                    next[(i*N) + j] = 'E';
                }

                if (next[(i*N) + j] == 'E')  empty_counter += 1;
            }
        }

        //swap matrice di scrittura con quella di lettura
        temp = recv;
        recv = next;
        next = temp;

        MPI_Gatherv(&recv[start_index*N], work_rows*N, MPI_CHAR, matrix, gat_counts, gat_displ, MPI_CHAR, 0, comm_topology);
    
        if(rank == 0){
            fprintf(file, "\n\nMatrice al %d passo discreto\n\n", k+1);
            print_graphic_matrix(matrix, M, N, file);
        }

        //controllo se tutte le celle sono vuote, ha senso solo quando grow prob -> 0
        MPI_Reduce(&empty_counter, &empty_matrix, 1, MPI_INT, MPI_SUM, 0, comm_topology);
        if (rank == 0){
            if (empty_matrix == M*N){
                empty_matrix = 1;
            } else {
                empty_matrix = 0;
            }
        }
        MPI_Bcast(&empty_matrix, 1, MPI_INT, 0, comm_topology);
        if (empty_matrix == 1)
            break;
        
    }
    if (rank == 0){
        fprintf(file, "\nEsperimento finito\n");
        free(matrix);
        fclose(file);
    }  

    free(recv);
    free(next);
    //MPI_Request_free(&request[0]);
    //MPI_Request_free(&request[1]);
    MPI_Finalize();
    return 0;
}