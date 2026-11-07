#ifndef FILE_H
#define FILE_H

#include <stdio.h>

int delete_directory(const char *directory);

FILE* open_file(char* filename);

#endif