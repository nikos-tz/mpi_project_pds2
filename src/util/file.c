#define _XOPEN_SOURCE 700

#include <src/util/file.h>
#include <stdlib.h>
#include <ftw.h>
#include <stdio.h>

static int remove_entry(const char *path, const struct stat *sb, int type, struct FTW *ftwbuf){
  (void)sb;
  (void)type;
  (void)ftwbuf;

  return remove(path);
}

int delete_directory(const char *directory){
  if(nftw(directory, remove_entry, 64, FTW_DEPTH | FTW_PHYS) == -1){
    perror("Error deleting directory");
    return -1;
  }

  return 0;
}

FILE* open_file(char* filename){
  FILE* fp = NULL;

  fp = fopen(filename, "rb");

  if(fp == NULL){
    perror("Error opening input file");
    exit(EXIT_FAILURE);
  }

  return fp;
}