#include <stdio.h>
#include <stdlib.h>

main()
{
    int tmp, nx, ny, intcount, extcount, i, j;
    char t[20];
    FILE *geom;
    geom = fopen("Geometry.txt", "r");
    fscanf(geom, "%s %s %d", &t, &t, &tmp);
    fscanf(geom, "%s %s %d", &t, &t, &tmp);
    fscanf(geom, "%s %s %s %s %d", &t, &t, &t, &t, &nx);
    fscanf(geom, "%s %s %s %s %d", &t, &t, &t, &t, &ny);
    int p[nx+1][ny+1];
    intcount = 1;
    extcount = nx + 2;
    for (i = 0; i < nx+1; i++)
    {
        for (j = 0; j < ny+1; j++)
        {
            if(i == 0 || j == 0 || i == nx || j == ny)
            {
                p[i][j] = extcount;
                extcount += 1;
            }
            else
            {
                p[i][j] = intcount;
                intcount += 1;
            }
        }
    }

    FILE *elnode = fopen("Element_Node.txt", "w");
    fprintf(elnode, "Element\tNodeA\tNodeB\tNodeC\tNodeD\n");

    int h = 0, v = 0;
    for (i = 0; i < nx*ny; i++)
    {
        if (h < nx)
        {
            fprintf(elnode, "%d\t", i+1);
            fprintf(elnode, "%d\t", p[v][h]);
            fprintf(elnode, "%d\t", p[v][h+1]);
            fprintf(elnode, "%d\t", p[v+1][h+1]);
            fprintf(elnode, "%d\n", p[v+1][h]);
            h += 1;
        }
        else
        {
            v += 1;
            h = 0;
            i = i - 1;
        }

    }

}
