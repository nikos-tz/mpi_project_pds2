CC = mpicc
CFLAGS = -O3
RM = rm -f

BUILD_DIR = build

EXECUTABLES_DATASETS = $(BUILD_DIR)/main_datasets.out
EXECUTABLES_RAND     = $(BUILD_DIR)/main_rand.out

DEPENDENCIES = src/distribute_by_median.c src/functions.c

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

datasets: $(EXECUTABLES_DATASETS)

$(BUILD_DIR)/main_datasets.out: $(BUILD_DIR) src/main_datasets.c $(DEPENDENCIES)
	$(CC) $(CFLAGS) $^ -o $@ -lm

rand: $(EXECUTABLES_RAND)

$(BUILD_DIR)/main_rand.out: $(BUILD_DIR) src/main_rand.c $(DEPENDENCIES)
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	$(RM) $(EXECUTABLES_RAND) $(EXECUTABLES_DATASETS)
