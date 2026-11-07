#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>

#include <src/util/output.h>


int write_proc_data(const char *directory, const char *filename, int n, int d, double **proc_data) {

    /* Create directory if it doesn't exist */
    if (mkdir(directory, 0777) != 0 && errno != EEXIST) {
        perror("Error creating directory");
        return -1;
    }

    /* Construct output path */
    char filepath[512];

    snprintf(filepath, sizeof(filepath), "%s/%s", directory, filename);

    FILE *fp = fopen(filepath, "wb");

    if (fp == NULL) {
        perror("Error opening output file");
        return -1;
    }

    int64_t d_file = d;
    int64_t n_file = n;

    /* Write header: D first, N second */
    if (fwrite(&d_file, sizeof(int64_t), 1, fp) != 1 ||
        fwrite(&n_file, sizeof(int64_t), 1, fp) != 1) {

        perror("Error writing file header");
        fclose(fp);
        return -1;
    }

    /* Write points */
    for (int i = 0; i < n; ++i) {

        if (fwrite(proc_data[i], sizeof(double), d, fp) != (size_t)d) {

            perror("Error writing point data");
            fclose(fp);
            return -1;
        }
    }

    if (fclose(fp) != 0) {
        perror("Error closing output file");
        return -1;
    }

    return 0;
}


int write_pivot(const char *filename, int d, double *pivot) {

    /* Construct output path */
    char filepath[256];

    snprintf(filepath, sizeof(filepath), "proc_data/%s", filename);

    FILE *fp = fopen(filepath, "wb");

    if (fp == NULL) {
        perror("Error opening pivot file");
        return -1;
    }

    int64_t d_file = d;
    int64_t n_file = 1;

    /* Write header: D first, N second */
    if (fwrite(&d_file, sizeof(int64_t), 1, fp) != 1 ||
        fwrite(&n_file, sizeof(int64_t), 1, fp) != 1) {

        perror("Error writing pivot header");
        fclose(fp);
        return -1;
    }

    /* Write pivot */
    if (fwrite(pivot, sizeof(double), d, fp) != (size_t)d) {

        perror("Error writing pivot data");
        fclose(fp);
        return -1;
    }

    if (fclose(fp) != 0) {
        perror("Error closing pivot file");
        return -1;
    }

    return 0;
}