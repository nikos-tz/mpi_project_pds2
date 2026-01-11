CC = mpicc
CFLAGS = -O3
RM = rm -f

EXECUTABLES_DATASETS = build/main_datasets.out

EXECUTABLES_RAND = build/main_rand.out

DEPENDENCIES = src/distribute_by_median.c src/functions.c

datasets: $(EXECUTABLES_DATASETS)

build/main_datasets.out:  src/main_datasets.c $(DEPENDENCIES)
	$(CC) $(CFLAGS) $^ -o $@ -lm

rand: $(EXECUTABLES_RAND)

build/main_rand.out:  src/main_rand.c $(DEPENDENCIES)
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	$(RM) $(EXECUTABLES_RAND) $(EXECUTABLES_DATASETS)