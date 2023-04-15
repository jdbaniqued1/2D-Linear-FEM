#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "Conjugate_Gradient.h"

int main()
{
    int g, h, i, j, k, nodenum, bcnum, a, b, elementnum;
    char t[20];

    FILE *input;
    input = fopen("Input.txt", "r");
    fscanf(input, "%s %d", &t, &nodenum);
    for (i = 0; i < 3; i++)
    {
        fscanf(input, "%s", &t);
    }
    for (i = 0; i < nodenum*3; i++)
    {
        fscanf(input, "%s", &t);
    }
    fscanf(input, "%s %n", &t, &elementnum);
    for (i = 0; i < 5; i++)
    {
        fscanf(input, "%s", &t);
    }
    for (i = 0; i < elementnum*5; i++)
    {
        fscanf(input, "%s", &t);
    }
    fscanf(input, "%s %s %d", &t, &t, &bcnum);
    for (i = 0; i < 3; i++)
    {
        fscanf(input, "%s", &t);
    }
    h = bcnum*2;
    g = nodenum*2 - h;

    Matrix Kgg, Kgh, Khg, Khh, fg, uh, Kggi, I, ug, B, K;
    Kgg = make_mat(g,g);
    Kgh = make_mat(g,h);
    Khg = make_mat(h,g);
    Khh = make_mat(h,h);
    /*fg = make_mat(g,1);
    uh = make_mat(h,1);
    ug = make_mat(g,1);
    B = make_mat(g,1);
    K = make_mat(nodenum*2, nodenum*2);

    int c[h];

    for (i = 0; i < h; i++)
    {
        fscanf(input, "%d %d %lf", &a, &b, &uh.mat[i][0]);
        if (b = 1)
            c[i] = a * 2 - 2;
        else
            c[i] = a * 2 - 1;
    }

    FILE *Kin;

    Kin = fopen("K_Global.txt", "r");
    fscanf(Kin, "%s %s", &t, &t);

    for (i = 0; i < nodenum * 2; i++)
    {
        for (j = 0; j < nodenum * 2; j++)
        {
            fscanf(Kin, "%lf", &K.mat[i][j]);
        }
    }

    print_mat(K);

    FILE *cond;
    cond = fopen("Conditions.txt", "r");
    fscanf(cond, "%s %s %s %s", &t, &t, &t, &t);

    for (i = 0; i < h; i++)
    {
        fscanf(cond, "%s %lf", &t, &uh.mat[i][0]);
    }

    fscanf(cond, "%s %s %s", &t, &t, &t);
    for (i = 0; i < g; i++)
    {
        fscanf(cond, "%s %lf", &t, &fg.mat[i][0]);
    }

    B = Matrix_Sub(fg,Matrix_Mult(Kgh,uh));

    ug = Conjugate_Gradient(Kgg,B);

    FILE * output = fopen("Solution.csv", "w");

    fprintf(output, "Solution\nNode,Value\n");
    for (i=0; i<ug.row; i++)
    {
        fprintf(output, "%d, %lf\n", i+1, ug.mat[i][0]);
    }

    printf("The solution can be accessed in Solution.csv.\nPress any key to continue.");
    getch();*/

}
