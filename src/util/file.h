#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <stdint.h>

int delete_directory(const char *directory);

FILE* open_file(char* filename);

#endif