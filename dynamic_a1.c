#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define WIDTH 1000
#define HEIGHT 1000
#define MAX_ITER 250


int mandelbrot(double x, double y) {

	int iter = 0;
	double zr = 0.0, zi = 0.0, zr_new, zi_new;

	while (zr * zr + zi * zi < 4.0 && iter < MAX_ITER) {
		zr_new = zr * zr - zi * zi + x;
		zi_new = 2.0 * zr * zi + y;
		zr = zr_new;
		zi = zi_new;
		iter++;
	}
    return iter;
}

int main(int argc, char **argv) {

	double ctime;
	double start_time = clock();
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double real_min = -2.0;
	double real_max = 2.0;
	double imag_min = -2.0;
	double imag_max = 2.0;
	
	double dx = (real_max - real_min) / WIDTH;
	double dy = (imag_max - imag_min) / HEIGHT;

	int n = HEIGHT / size;
	int start_row = rank * n;
	int end_row = start_row + n;

	if (rank == size - 1)
		end_row = HEIGHT;

	int *buffer = (int *) malloc(WIDTH * n * sizeof(int));

	for (int i = start_row; i < end_row; i++) {
		for (int j = 0; j < WIDTH; j++) {
		    double x = real_min + j * dx;
		    double y = imag_min + i * dy;
		    buffer[(i - start_row) * WIDTH + j] = mandelbrot(x, y);
		}
	}
 

	MPI_Request request;
	
	if (rank == 0) {
		int *output = (int *) malloc(WIDTH * HEIGHT * sizeof(int));

		for (int i = 0; i < n * WIDTH; i++)
			output[i] = buffer[i];

		double cstart = clock();

		for (int i = 1; i < size; i++) {
			MPI_Irecv(output + i * n * WIDTH, n * WIDTH, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, MPI_STATUS_IGNORE);
		}

		double cend = clock();
		ctime = cend-cstart;

		FILE *fp = fopen("mandelbrot2.pgm", "wb");
		fprintf(fp, "P2\n%d %d\n%d\n", WIDTH, HEIGHT, MAX_ITER);
		for (int row = 0; row < HEIGHT; row++) {
			for (int col = 0; col < WIDTH; col++)
				fprintf(fp, "%d ", output[row * WIDTH + col]);
		    fprintf(fp, "\n");
		}

		fclose(fp);
		free(output);

	}
	else {
		double cstart2 = MPI_Wtime();
		
		MPI_Isend(buffer, n * WIDTH, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
		
		MPI_Wait(&request, MPI_STATUS_IGNORE);
		
		double cend2 = MPI_Wtime();
		ctime += cend2-cstart2;
	}

	double end_time = clock();
	double elapsed_time = (end_time-start_time) / CLOCKS_PER_SEC;
	printf("Elapsed time: %f seconds\n", elapsed_time);
	free(buffer);
	
	ctime = ctime/CLOCKS_PER_SEC;
	printf("Communication time: %f seconds\n", ctime);

	MPI_Finalize();

	return 0;

}
