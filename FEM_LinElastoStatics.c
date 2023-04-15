#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

int main()
{
    FILE *in;
    in = fopen("Input.txt", "r");
    char tem[20];
    int nodenum, elementnum, temp, i, j, e, k, ngp = 3, r, m;

    fscanf(in, "%s %d", &tem, &nodenum);

    double node[nodenum][2];
    fscanf(in, "%s %s %s", &tem, &tem, &tem);
    for (i = 0; i < nodenum; i++)
    {
        fscanf(in, "%d", &temp);
        fscanf(in, "%lf %lf", &node[temp-1][0], &node[temp-1][1]);
    }

    fscanf(in, "%s %d", &tem, &elementnum);
    printf("%d\n", elementnum);

    int element[elementnum][4];
    fscanf(in, "%s %s %s %s %s", &tem, &tem, &tem, &tem, &tem);
    for (i = 0; i < elementnum; i++)
    {
        fscanf(in, "%d", &temp);
        fscanf(in, "%d %d %d %d", &element[temp-1][0], &element[temp-1][1], &element[temp-1][2], &element[temp-1][3]);
    }

    ngp = 3;
    double s, t, p, q, X, Y, scalar;
    double Kloc[elementnum][8][8], xi[ngp], phi[4], dphids[4], dphidt[4], x[elementnum][4], y[elementnum][4], W[ngp];

    Matrix J, Ji, dphi, B, M, I, Bnew, BTnew, D, T;
    J = make_mat(2,2);
    Ji = make_mat(2,2);
    dphi = make_mat(2,4);
    B = make_mat(2,4);
    M = make_mat(4,4);
    I = make_mat(2,2);
    Bnew = make_mat(3,8);
    BTnew = make_mat(8,3);
    D = make_mat(3,3);
    T = make_mat(8,3);

    I.mat[0][0] = 1;
    I.mat[0][1] = 0;
    I.mat[1][0] = 0;
    I.mat[1][1] = 1;

    double E, v, factor;
    E = 200;
    v = 0.3;
    factor = E*pow(1+v,-1);
    D.mat[0][0] = factor;
    D.mat[0][1] = factor * v;
    D.mat[1][0] = factor * v;
    D.mat[1][1] = factor;
    D.mat[2][2] = factor * (1-v) * 0.5;

    p = 0;
    q = 0;
    r = 0;

    for (i = 0; i < elementnum; i++)
    {
        for (j = 0; j < 8; j++)
        {
            for (k = 0; k < 8; k++)
            {
                Kloc[i][j][k] = 0;
            }
        }
    }

    xi[0] = 0.7745966692;
    xi[1] = -0.7745966692;
    xi[2] = 0;

    W[0] = 0.5555555556;
    W[1] = 0.5555555556;
    W[2] = 0.8888888889;

    for (e = 0; e < elementnum; e++)
    {
        for (i = 0; i < ngp; i++)
        {
            for (j = 0; j < ngp; j++)
            {
                s = xi[i];
                t = xi[j];

                phi[0] = (1-s)*(1-t)*0.25;
                phi[1] = (1+s)*(1-t)*0.25;
                phi[2] = (1+s)*(1+t)*0.25;
                phi[3] = (1-s)*(1+t)*0.25;

                dphi.mat[0][0] = 0.25*(-1)*(1-t);
                dphi.mat[0][1] = 0.25*(1-t);
                dphi.mat[0][2] = 0.25*(1+t);
                dphi.mat[0][3] = 0.25*(-1)*(1+t);

                dphi.mat[1][0] = 0.25*(-1)*(1-s);
                dphi.mat[1][1] = 0.25*(-1)*(1+s);
                dphi.mat[1][2] = 0.25*(1+s);
                dphi.mat[1][3] = 0.25*(1-s);

                X = 0;
                Y = 0;
                for (r = 0; r < 4; r++)
                {
                    X += x[e][r] * phi[r];
                    Y += y[e][r] * phi[r];
                }

                for (m = 0; m < 2; m++)
                {
                    for (r = 0; r < 2; r++)
                    {
                        J.mat[m][r] = 0;
                    }
                }

                x[e][0] = node[element[e][0]-1][0];
                y[e][0] = node[element[e][0]-1][1];
                x[e][1] = node[element[e][1]-1][0];
                y[e][1] = node[element[e][1]-1][1];
                x[e][2] = node[element[e][2]-1][0];
                y[e][2] = node[element[e][2]-1][1];
                x[e][3] = node[element[e][3]-1][0];
                y[e][3] = node[element[e][3]-1][1];

                for (r = 0; r < 4; r++)
                {
                    J.mat[0][0] += x[e][r] * dphi.mat[0][r];
                    J.mat[0][1] += y[e][r] * dphi.mat[0][r];
                    J.mat[1][0] += x[e][r] * dphi.mat[1][r];
                    J.mat[1][1] += y[e][r] * dphi.mat[1][r];
                }

                Ji = Gauss_Jordan(J,I);

                B = Matrix_Mult(Ji,dphi);

                Bnew.mat[0][0] = B.mat[0][0];
                Bnew.mat[0][2] = B.mat[0][1];
                Bnew.mat[0][4] = B.mat[0][2];
                Bnew.mat[0][6] = B.mat[0][3];
                Bnew.mat[1][1] = B.mat[1][0];
                Bnew.mat[1][3] = B.mat[1][1];
                Bnew.mat[1][5] = B.mat[1][2];
                Bnew.mat[1][7] = B.mat[1][3];
                Bnew.mat[2][0] = B.mat[1][0];
                Bnew.mat[2][2] = B.mat[1][1];
                Bnew.mat[2][4] = B.mat[1][2];
                Bnew.mat[2][6] = B.mat[1][3];
                Bnew.mat[2][1] = B.mat[0][0];
                Bnew.mat[2][3] = B.mat[0][1];
                Bnew.mat[2][5] = B.mat[0][2];
                Bnew.mat[2][7] = B.mat[0][3];
                print_mat(Bnew);
                printf("\n\n");

                scalar = W[i]*W[j]*det(J);

                BTnew = Matrix_Trans(Bnew);
                T = Matrix_Mult(BTnew,D);
                M = Matrix_Scalar(Matrix_Mult(T,Bnew),scalar);

                //print_mat(M);


                for (m = 0; m < 8; m++)
                {
                    for (r = 0; r < 8; r++)
                    {
                        Kloc[e][m][r] += M.mat[m][r];
                    }
                }

                printf("e = %d\n", e);
            }
        }
    }

    FILE *out;

    out = fopen("K_Local.txt", "w");

    fprintf(out, "n\t%d\n", nodenum * 2);

    for (i = 0; i < elementnum; i++)
    {
        fprintf(out, "Element %d\n", i+1);
        for (j = 0; j < 8; j++)
        {
            for (k = 0; k < 8; k++)
            {
                fprintf(out, "%lf\t", Kloc[i][j][k]);
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }

    printf("The output file is K_Local.txt. Press any key to end.\n");
    getch();

    return 0;
}
