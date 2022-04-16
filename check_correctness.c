#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
    
    FILE *file_par, *file_seq;

    file_par = fopen("correttezza_par.txt", "r");
    file_seq = fopen("correttezza_seq.txt", "r");

    int nrow = 0;

    char seq_row[4096];
    char par_row[4096];
    
    if ((file_par = fopen("correttezza_par.txt", "r")) != NULL && 
    (file_seq = fopen("correttezza_seq.txt", "r")) != NULL){
        
        while((fgets(par_row, 4096, file_par) != NULL) &&
        (fgets(seq_row, 4096, file_seq) != NULL)){
            if(strcmp(par_row, seq_row) != 0){
                printf("Differenza nella riga %d!\n", nrow);
                fclose(file_par);
                fclose(file_seq);
                return 0;
            }
            nrow++;
        }

    } else {
        printf("FILE OPENING ERROR!\n");
    }
    printf("Tutto ok!\n");
    
    fclose(file_par);
    fclose(file_seq);

    return 0;
}