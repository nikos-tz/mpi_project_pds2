#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <stdint.h>

int delete_directory(const char *directory);

FILE* open_file(char* filename);

void file_read_int64(int64_t *value, size_t count, FILE *fp);

void scan_one_int(int* value);

#endif