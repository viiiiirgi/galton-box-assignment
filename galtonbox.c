#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void galton_board(int n, int N, int *bins) {
    
    // reset the bins
    for (int i = 0; i <= n; i++) {
        bins[i] = 0;
    }

    // simulate N balls dropping through n levels
    for (int i = 0; i < N; i++) {
        int left_moves = 0;

        for (int level = 0; level < n; level++) {

            // move left with probability 0.5, or move right with probability 0.5
            if ((double)rand() / RAND_MAX < 0.5) {
                left_moves++; // move left = stay in the same column
            }
        }

        //the final position is (i, n - i), where i is the number of left moves
        bins[left_moves]++;
    }
}


void binomial_distribution(int n, int N, double *binomial_probability) {

    for (int i = 0; i <= n; i++) {
        
        double binomial_coefficient = 1;
        for (int j = 0; j < i; j++) {
            binomial_coefficient *= (double)(n - j) / (j + 1);
        }
        
        binomial_probability[i] = binomial_coefficient * pow(0.5, n);
    }
}


void normal_distribution(int n, int N, double *normal_pdf) {
    double mu = n / 2.0;
    double sigma = sqrt(n / 4.0);
    
    for (int i = 0; i <= n; i++) {
        double x = i;
        double coefficient = 1.0 / (sqrt(2 * M_PI) * sigma);
        normal_pdf[i] = coefficient * exp(-0.5 * pow((x - mu) / sigma, 2));

        normal_pdf[i] *= N;
    }
}


double compute_mse(int n, double *binomial_probability, double *normal_pdf) {
    double mse = 0.0;
    
    for (int i = 0; i <= n; i++) {
        double error = binomial_probability[i] - normal_pdf[i];
        mse += error * error; 
    }
    
    return mse / (n + 1); 
}

void compute_mse_per_bin(int N, int n, int *bins, double *binomial_probability, double *normal_pdf) {
    printf("Position\tEmpirical\tBinomial\tNormal\tMSE(Binomial)\tMSE(Normal)\n");
    
    for (int i = 0; i <= n; i++) {
        int empirical_count = bins[i];
        
        double theoretical_binomial = binomial_probability[i] * N;  
        double theoretical_normal = normal_pdf[i];

        double mse_binomial = (empirical_count - theoretical_binomial) * (empirical_count - theoretical_binomial);
        double mse_normal = (empirical_count - theoretical_normal) * (empirical_count - theoretical_normal);

        printf("%d\t\t%d\t\t%.6f\t%.6f\t%.6f\t\t%.6f\n", i, empirical_count, theoretical_binomial, theoretical_normal, mse_binomial, mse_normal);
    }
}

void save_data_to_file(int N, int n, int *bins, double *binomial_probability, double *normal_pdf) {
    FILE *file = fopen("galton_data.txt", "w");
    if (file == NULL) {
        printf("Error opening file\n");
        return;
    }

    fprintf(file, "#Position\tEmpirical\tBinomial\tNormal\n");
    for (int i = 0; i <= n; i++) {
        fprintf(file, "%d\t%d\t%.6f\t%.6f\n", i, bins[i], binomial_probability[i] * N, normal_pdf[i]);
    }

    fclose(file);
}


void generate_gnuplot_script(int n) {
    FILE *gnuplot_file = fopen("galton_plot.gnu", "w");
    if (gnuplot_file == NULL) {
        printf("Error opening Gnuplot script file\n");
        return;
    }

    fprintf(gnuplot_file, "set title 'Galton Board Simulation (n=%d)'\n", n);
    fprintf(gnuplot_file, "set xlabel 'Position on Board (bin)'\n");
    fprintf(gnuplot_file, "set ylabel 'Number of Balls'\n");
    fprintf(gnuplot_file, "set grid\n");
    fprintf(gnuplot_file, "set style data histograms\n");
    fprintf(gnuplot_file, "set style histogram clustered gap 1\n");
    fprintf(gnuplot_file, "set boxwidth 0.9\n");
    fprintf(gnuplot_file, "set style fill solid\n");
    fprintf(gnuplot_file, "plot 'galton_data.txt' using 2:xtic(1) title 'Empirical' with boxes, \\\n");
    fprintf(gnuplot_file, "     'galton_data.txt' using 3:xtic(1) title 'Binomial' with linespoints lt 1 lw 2 lc rgb 'red', \\\n");
    fprintf(gnuplot_file, "     'galton_data.txt' using 4:xtic(1) title 'Normal' with linespoints lt 2 lw 2\n");

    fclose(gnuplot_file);
}


void run_gnuplot() {
    system("gnuplot -p galton_plot.gnu");
}

int main(int argc, char *argv[]) {

    // check if there is correct number of arguments
    if (argc != 3) {
        printf("There should be  %s arguments, add the number of balls and of levels \n", argv[0]);
        return 1;
    }

    // the number of balls and levels should be positive
    int N = atoi(argv[1]);
    int n = atoi(argv[2]);

    if (N <= 0 || n <= 0) {
        printf("Error: Number of balls and levels must be positive integers.\n");
        return 1;
    }

    // seed random number generator
    srand(time(NULL));


    // arrays to store results
    int *bins = (int *)malloc((n + 1) * sizeof(int));
    double *binomial_probability = (double *)malloc((n + 1) * sizeof(double)); 
    double *normal_pdf = (double *)malloc((n + 1) * sizeof(double)); 

    // simulate the galton board
    galton_board(n, N, bins);

    // calculate the binomial distribution (PMF)
    binomial_distribution(n, N, binomial_probability);

    // calculate the normal distribution (PDF)
    normal_distribution(n, N, normal_pdf);

    // compute the mean quadratic error (MSE)
    double mse = compute_mse(n, binomial_probability, normal_pdf);
    printf("Mean Squared Error between the Binomial PMF and Normal PDF: %.6f\n", mse);

    compute_mse_per_bin(N, n, bins, binomial_probability, normal_pdf);

    save_data_to_file(N, n, bins, binomial_probability, normal_pdf);
    generate_gnuplot_script(n);
    run_gnuplot();
  
    free(bins);
    free(binomial_probability);
    free(normal_pdf);

    return 0;
}