#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}

void fft(complex s[], complex t[], int n, int sign) {
    // 1. Construir um sinal sp de tamanho n/2 cujos elementos são
    //    aqueles nos índices pares de s
    complex sp[n/2];
    for (int i = 0; i < n; i++){
        if (i % 2 == 0){
            sp[i/2] = s[i];
        }
    }

    // 2. Construir um sinal si de tamanho n/2 cujos elementos são
    //    aqueles nos índices ímpares de s
    complex si[n/2];
    for (int i = 0; i < n; i++){
        if (i % 2 != 0){
            si[(i-1)/2] = s[i]
        }
    }

    complex tp[n/2];
    complex ti[n/2];
    // 3. Obter recursivamente a transformada de Fourier tp do sinal sp
    fft(sp, tp, n/2, sign);

    // 4. Obter recursivamente a transformada de Fourier ti do sinal ti
    fft(si, ti, n/2, sign);

    // 5. Olhar no enunciado (Construir um sinal t de tamanho n)


    for (int k = 0; k < n/2; k++){
        double x = sign * 2 * PI * k / n;

        double cosx = cos(x);
        double sinx = sin(x);

        t[k].a = tp[k].a + (ti[k].a * cosx - ti[k].b * sinx);
        t[k].b = tp[k].b + (ti[k].a * sinx + ti[k].b * cosx);

        t[k+(n/2)].a = tp[k].a - (ti[k].a * cosx - ti[k].b * sinx);
        t[k+(n/2)].b = tp[k].a - (ti[k].a * sinx + ti[k].b * cosx);
    }
}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
