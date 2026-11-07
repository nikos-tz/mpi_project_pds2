CC = mpicc
CFLAGS = -O3 -I.
RM = rm -f

BUILD_DIR = build

EXECUTABLES_DATASETS = $(BUILD_DIR)/main_datasets.out
EXECUTABLES_RAND     = $(BUILD_DIR)/main_rand.out

DEPENDENCIES = src/distribute_by_median.c src/aux/math.c src/aux/memory.c src/aux/sort.c src/util/output.c src/util/file.c

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.PHONY: all clean datasets rand

all: datasets rand

datasets: $(EXECUTABLES_DATASETS)

$(BUILD_DIR)/main_datasets.out: src/main_datasets.c $(DEPENDENCIES) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $^ -o $@ -lm

rand: $(EXECUTABLES_RAND)

$(BUILD_DIR)/main_rand.out: src/main_rand.c $(DEPENDENCIES) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	$(RM) $(EXECUTABLES_RAND) $(EXECUTABLES_DATASETS)
