#define _XOPEN_SOURCE 700

#include "delete_directory.h"

#include <ftw.h>
#include <stdio.h>

static int remove_entry(const char *path,
                        const struct stat *sb,
                        int type,
                        struct FTW *ftwbuf)
{
    (void)sb;
    (void)type;
    (void)ftwbuf;

    return remove(path);
}

int delete_directory(const char *directory)
{
    if (nftw(directory, remove_entry, 64, FTW_DEPTH | FTW_PHYS) == -1) {
        perror("Error deleting directory");
        return -1;
    }

    return 0;
}