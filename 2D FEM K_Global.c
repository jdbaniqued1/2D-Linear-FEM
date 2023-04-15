#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main()
{
    char t[10];
    int i, j, k, tmp, gather, n, elementnum;
    double temp;
    FILE *klocal, *input;
    klocal = fopen("K_Local.txt", "r");
    fscanf(klocal, "%s %d", &t, &n);

    Matrix K_Local, K_Global, K_Gtemp, L, LT;

    K_Local = make_mat(8,8);
    K_Global = make_mat(n,n);
    K_Gtemp = make_mat(n,n);
    L = make_mat(8,n);
    LT = make_mat(n,8);

    input = fopen("Input.txt", "r");
    for (i = 0; i < 5; i++)
    {
        fscanf(input, "%s", &t);
    }
    for (i = 0; i < n*3*0.5; i++)
    {
        fscanf(input, "%s", &t);
    }

    fscanf(input, "%s %d", &t, &elementnum);
    fscanf(input, "%s %s %s %s %s", &t, &t, &t, &t, &t);

    for(i = 0; i < elementnum; i++)
    {
        fscanf(input, "%s", &t);
        for(j = 0; j < 8; j++)
        {
            fscanf(input, "%d", &gather);
            L.mat[j][gather*2-2] = 1;
            L.mat[j][gather*2-1] = 1;
        }
        LT = Matrix_Trans(L);

        fscanf(klocal, "%s", &t);
        fscanf(klocal, "%d", &tmp);
        for(j = 0; j < 8; j++)
        {
            for(k = 0; k < 8; k++)
            {
                fscanf(klocal, "%lf", &K_Local.mat[j][k]);
            }
        }

        K_Gtemp = Matrix_Mult(Matrix_Mult(LT,K_Local),L);

        K_Global = Matrix_Add(K_Global, K_Gtemp);

        zero_mat(L);
        zero_mat(LT);
        zero_mat(K_Gtemp);

    }

    FILE *out;

    out = fopen("K_Global.txt", "w");

    fprintf(out, "K_Global Matrix\n");
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k++)
        {
            fprintf(out, "%lf\t", K_Global.mat[j][k]);
        }
        fprintf(out, "\n");
    }

    printf("The output file is K_Global.txt. Press any key to end.\n");
    getch();

    return 0;
}
