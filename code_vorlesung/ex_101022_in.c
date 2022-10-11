#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
    printf("%d\n", numerow(0., 0.1, 1, 1000, 0., 10.));
}

double numerov(float u_0, float u_1, int K, unsigned int N, float x_min, float x_max)
{
    double *u = malloc(N * sizeof(double) + 2 * sizeof(double));
    u[0] = u_0;
    u[1] = u_1;
    double h = (x_min - x_max)  / N;
    for (int i = 2; i < N; i++)
    {
        u[i] = (2 * u[i - 1] * (1 - (5 / 12) * pow(h, 2) * K) - u[i - 2] * (1 + (pow(h, 2) / 12) * K)) / (1 + 1 / 12 * pow(h, 2) * K);
    }
    return *u;
}