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
    if (n == 1){ // Condicao de parada da recursao (quando o vetor s possuir tamanho 1)
        t[0] = s[0]; 
        return ;
    }

    // Constroi o sinal sp de tamanho n/2 (apenas com os indices pares de s)
    complex sp[n/2];
    for (int i = 0; i < n; i++){
        if (i % 2 == 0){
            sp[i/2] = s[i];
        }
    }

    // Constroi o sinal si de tamanho n/2 (apenas com os indices impares de s)
    complex si[n/2];
    for (int i = 0; i < n; i++){
        if (i % 2 != 0){
            si[(i-1)/2] = s[i];
        }
    }

    // Declara tp e ti
    complex tp[n/2];
    complex ti[n/2];
    // Calcula tp e ti com a recursividade
    fft(sp, tp, n/2, sign); // (Observe que diminuimos o tamanho do vetor pela metade)
    fft(si, ti, n/2, sign); 


    // Aqui calculamos o sinal t de tamanho n
    for (int k = 0; k < n/2; k++){
        double x = sign * 2 * PI * k / n; // Declara x

        double cosx = cos(x);
        double sinx = sin(x);

        t[k].a = tp[k].a + (ti[k].a * cosx - ti[k].b * sinx); // Parte real de t (do indice 0 ao (n/2)-1)
        t[k].b = tp[k].b + (ti[k].a * sinx + ti[k].b * cosx); // Parte imaginaria de t (do indice 0 ao (n/2)-1)

        t[k+(n/2)].a = tp[k].a - (ti[k].a * cosx - ti[k].b * sinx); // Parte real de t (do indice n/2 ao n-1)
        t[k+(n/2)].b = tp[k].b - (ti[k].a * sinx + ti[k].b * cosx); // Parte imaginaria de t (do indice n/2 ao n-1)
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
    complex l[width];

    // Para cada linha da matriz
    for (int j = 0; j < height; j++){ // percorro linha por linha
        // Cria sinal com dados da linha
        for (int i = 0; i < width; i++){ 
            l[i] = matrix[j][i];
        }

        // Calcula transformada
        complex t[width];
        fft_forward(l, t, width);
        
        // Atualiza a propria matriz com a transformada
        for (int i = 0; i < width; i++){ 
            matrix[j][i] = t[i]; 
        }
    }

    complex c[height];

    // Para cada coluna da matriz
    for (int i = 0; i < width; i++){
        // Cria sinal com dados da coluna
        for (int j = 0; j < height; j++){
            c[j] = matrix[j][i];
        }

        // Calcula transformada
        complex t[height];
        fft_forward(c, t, height);

        // Atualiza a propria matriz com a transformada
        for (int j = 0; j < height; j++){
            matrix[j][i] = t[j];
        }
    } 
    
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    complex c[height];

    // Para cada coluna da matriz
    for (int i = 0; i < width; i++){
        // Cria sinal com dados da coluna
        for (int j = 0; j < height; j++){
            c[j] = matrix[j][i];
        }

        // Calcula transformada
        complex t[height];
        fft_inverse(c, t, height);

        // Atualiza a propria matriz com a transformada
        for (int j = 0; j < height; j++){
            matrix[j][i] = t[j];
        }
    }     

    complex l[width];

    // Para cada linha da matriz
    for (int j = 0; j < height; j++){ // percorro linha por linha
        // Cria sinal com dados da linha
        for (int i = 0; i < width; i++){ 
            l[i] = matrix[j][i]; 
        }

        // Calcula transformada
        complex t[width];
        fft_inverse(l, t, width);
        
        // Atualiza a propria matriz com a transformada
        for (int i = 0; i < width; i++){ 
            matrix[j][i] = t[i]; 
        }
    }

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
