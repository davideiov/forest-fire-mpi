#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define ignite_prob 15  //A tree ignites with this probability (f), 0 <= ignite_prob <= 100
#define grow_prob 5     //A tree grows on an empty cell with this probability (p), 0 <= grow_prob <= 100
#define gen_prob 80     //Probability that a tree grows, used for the matrix generation, 0 <= gen_prob <= 100
#define S 250           //Number of steps

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
    if (((i*n) + j + n) <= nrows*n && recv[((i*n) + j + n)] == 'B')    
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

//controlla la probabilità p di crescita di un albero in una cella vuota
void check_tree_grows(char *next, int i, int j, int n){
    int random = rand() % 100;
    if (random < 0) 
        random *= -1;

    if (random <= grow_prob && random > 0)
        next[(i*n) + j] = 'T';
    else
        next[(i*n) + j] = 'E';
}

int main(int argc, char **argv){
    
    char *matrix, *next, *temp;
    int k, empty_counter = 0, empty_matrix = 0, N = atoi(argv[1]);

    //variabili per prestazioni
    double start, end;

    //generazione matrice iniziale
    matrix = (char *) malloc(sizeof(char) * N * N);
    next = (char *) malloc(sizeof(char) * N * N);
        
    matrix_init(matrix, N, N);

    start = MPI_Wtime();

    //simulazione su S passi discreti
    for(k=0; k<S; k++){
        empty_counter = 0;

        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                if (matrix[(i*N) + j] == 'E'){
                    check_tree_grows(next, i, j, N);
                } else if (matrix[(i*N) + j] == 'T') {
                    check_neighbors_ignite(matrix, next, i, j, N, N);
                } else if (matrix[(i*N) + j] == 'B') {
                    next[(i*N) + j] = 'E';
                }
                
                if (next[(i*N) + j] == 'E')  empty_counter += 1;
            }
        }

        //swap matrice di scrittura con quella di lettura
        temp = matrix;
        matrix = next;
        next = temp;

        if (empty_counter == N*N){
            empty_matrix = 1;
        } else {
            empty_matrix = 0;
        }
    
        if (empty_matrix == 1)
            break;
    }

    end = MPI_Wtime();
    printf("Esperimento completato in %d passi discreti in %.4fs.\n", k, end-start);
    
    free(matrix);
    free(next);
    
    return 0;
}