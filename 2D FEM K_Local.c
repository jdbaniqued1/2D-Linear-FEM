#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

int main()
{
    FILE *geom;
    geom = fopen("Geometry.txt", "r");
    char tem[20];
    int n,m, nx, ny, Lx, Ly, i, j, k, ngp, r, e;
    n = 4;
    ngp = 3;

    fscanf(geom, "%s %s", &tem, &tem);
    fscanf(geom, "%d", &Lx);
    fscanf(geom, "%s %s", &tem, &tem);
    fscanf(geom, "%d", &Ly);
    fscanf(geom, "%s %s %s %s", &tem, &tem, &tem, &tem);
    fscanf(geom, "%d", &nx);
    fscanf(geom, "%s %s %s %s", &tem, &tem, &tem, &tem);
    fscanf(geom, "%d", &ny);

    double s, t, p, q, X, Y, scalar;
    double Kloc[nx*ny][n][n], xi[ngp], phi[4], dphids[4], dphidt[4], x[nx*ny][4], y[nx*ny][4], W[ngp];

    Matrix J, Ji, dphi, B, BT, T, M, D, Bnew, Dnew;
    J = make_mat(2,2);
    Ji = make_mat(2,2);
    dphi = make_mat(2,4);
    B = make_mat(2,4);
    BT = make_mat(4,2);
    T = make_mat(4,2);
    M = make_mat(4,4);
    D = make_mat(2,2);
    Bnew = make_mat(3,8);
    Dnew = make_mat(3,3);

    D.mat[0][0] = 1;
    D.mat[0][1] = 0;
    D.mat[1][0] = 0;
    D.mat[1][1] = 1;

    double E, v, factor;
    E = 200;
    v = 0.3;
    factor = E*pow(1+v,-1);
    Dnew.mat[0][0] = factor;
    Dnew.mat[0][1] = factor * v;
    Dnew.mat[1][0] = factor * v;
    Dnew.mat[1][1] = factor;
    Dnew.mat[2][2] = factor * (1-v) * 0.5;

    print_mat(Dnew);
    printf("\n\n");

    p = 0;
    q = 0;
    r = 0;
    for (i = 0; i < ny; i++)
    {
        for (j = 0; j < nx; j++)
        {
            x[r][0] = p;
            x[r][1] = p + Lx * pow(nx,-1);
            x[r][2] = p + Lx * pow(nx,-1);
            x[r][3] = p;

            y[r][0] = q;
            y[r][1] = q;
            y[r][2] = q + Ly * pow(ny,-1);
            y[r][3] = q + Ly * pow(ny,-1);

            p = p + Lx * pow(nx,-1);
            r = r+1;
        }
        p = 0;
        q = q + Ly * pow(ny,-1);
    }

    for (i = 0; i < nx*ny; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
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

    for (e = 0; e < nx*ny; e++)
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

                for (r = 0; r < 4; r++)
                {
                    J.mat[0][0] += x[e][r] * dphi.mat[0][r];
                    J.mat[0][1] += y[e][r] * dphi.mat[0][r];
                    J.mat[1][0] += x[e][r] * dphi.mat[1][r];
                    J.mat[1][1] += y[e][r] * dphi.mat[1][r];
                }

                Ji = Gauss_Jordan(J,D);

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

                BT = Matrix_Trans(B);
                T = Matrix_Mult(BT,D);
                M = Matrix_Scalar(Matrix_Mult(T,B),scalar);

                for (m = 0; m < n; m++)
                {
                    for (r = 0; r < n; r++)
                    {
                        Kloc[e][m][r] += M.mat[m][r];
                    }
                }
            }
        }
    }

    FILE *out;

    out = fopen("K_Local.txt", "w");

    fprintf(out, "n\t%d\n", (nx+1)*(ny+1));

    for (i = 0; i < nx*ny; i++)
    {
        fprintf(out, "Element %d\n", i+1);
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
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
