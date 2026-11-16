CC = mpicc
CFLAGS = -O3 -I.
RM = rm -f

BUILD_DIR = build

EXECUTABLES_DATASETS = $(BUILD_DIR)/main.out

DEPENDENCIES = src/distribute_by_median.c src/aux/math.c src/aux/memory.c src/aux/sort.c src/util/output.c src/util/file.c

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.PHONY: all clean datasets

all: datasets

datasets: $(EXECUTABLES_DATASETS)

$(BUILD_DIR)/main.out: src/main.c $(DEPENDENCIES) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	$(RM) $(EXECUTABLES_DATASETS)
