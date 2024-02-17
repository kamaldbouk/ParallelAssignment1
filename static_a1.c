#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define WIDTH 1200
#define HEIGHT 1200
#define max_iter 250

 
int main(int argc, char** argv) {

	double ctime; //communivation time
	double start_time = clock();
	int rank, size;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double real_min = -2.0;
	double real_max = 1.0;
	double imag_min = -1.0;
	double imag_max = 1.0;

	int n = HEIGHT / size;
	int start_row = n * rank;
	int end_row = start_row + n;

	if (rank == size - 1)
		end_row = HEIGHT;

	int* row = (int*)malloc(sizeof(int) * WIDTH);
	int* data = (int*)malloc(sizeof(int) * WIDTH * n);

	for (int i = start_row; i < end_row; i++) {
		for (int j = 0; j < WIDTH; j++) {
			double c_real = real_min + (real_max-real_min) * j/WIDTH;
			double c_imag = imag_min + (imag_max-imag_min) * i / HEIGHT;
			double z_real = 0.0;
			double z_imag = 0.0;
			int iter = 0;

			while (z_real * z_real + z_imag * z_imag < 4.0 && iter < max_iter) {
				double next_z_real = z_real * z_real - z_imag * z_imag + c_real;
				double next_z_imag = 2.0 * z_real * z_imag + c_imag;

				z_real = next_z_real;
				z_imag = next_z_imag;
				iter++;
			}

			if (iter == max_iter)
				row[j] = 0;
			else
				row[j] = iter % 256;
		
		}
	
	int row_index = (i - start_row) * WIDTH;
        for (int j = 0; j < WIDTH; j++)
            data[row_index + j] = row[j];
    	}

	free(row);

	int* final_data;
	if (rank == 0)
	final_data = (int*)malloc(sizeof(int) * WIDTH * HEIGHT);


	double cstart = clock();

	MPI_Gather(data, WIDTH * n, MPI_INT, final_data, WIDTH * n, MPI_INT, 0, MPI_COMM_WORLD);

	double cend = clock();

	ctime = (cend-cstart)/CLOCKS_PER_SEC;

	free(data);

	if (rank == 0) {
		FILE* fp = fopen("mandelbrot.pgm", "wb");
		fprintf(fp, "P5\n%d %d\n255\n", WIDTH, HEIGHT);
		
		for (int i = 0; i < WIDTH * HEIGHT; i++)
	    		fputc(final_data[i], fp);

		fclose(fp);
		free(final_data);
	}

	double end_time = clock();
	double elapsed_time = (end_time-start_time) / CLOCKS_PER_SEC;
	
	printf("Elapsed Time: %f seconds.\n", elapsed_time);
	printf("Communication Time: %f seconds.\n", ctime);

	MPI_Finalize();

	return 0;
}
