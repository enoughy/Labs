#include <cstdlib>
#include <fstream>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

constexpr int K = 100;
constexpr int N_MIN = 256;             // 1 KB
constexpr int N_MAX = 8 * 1024 * 1024; // 16 MB
constexpr int COUNT = 1;

void initialize_forward(std::vector<int> &arr) {
  int size = arr.size();
  for (int i = 0; i < size - 1; i++) {
    arr[i] = i + 1;
  }
  arr[size - 1] = 0;
}

void initialize_reverse(std::vector<int> &arr) {
  int size = arr.size();
  for (int i = size - 1; i > 0; i--) {
    arr[i] = i - 1;
  }
  arr[0] = size - 1;
}

void initialize_random(std::vector<int> &arr) {
  int size = arr.size();
  for (int i = 0; i < size; i++) {
    arr[i] = i;
  }

  for (int i = size - 1; i > 0; i--) {
    std::swap(arr[i], arr[rand() % (i + 1)]);
  }
}

long double measure_latency(std::vector<int> &arr, int size,
                            std::ofstream &log_file) {
  long double avg_latency;
  long double min_latency = std::numeric_limits<long double>::max();
  for (int i = 0; i < COUNT; i++) {
    volatile int x = 0;

    // Прогрев кэш
    for (int i = 0; i < size; i++) {
      x = arr[x];
    }

    auto start_time = __rdtsc();
    for (int i = 0; i < size * K; i++) {
      x = arr[x];
    }
    x = 0;
    auto end_time = __rdtsc();

    avg_latency = (end_time - start_time) / (long double)(size * K);

    min_latency = std::min(min_latency, avg_latency);
  }
  return min_latency;
}

void tieme_array_traversal(std::ofstream &output) {

  for (int size = N_MIN; size <= N_MAX; size *= 2) {

    std::vector<int> arr(size);

    initialize_forward(arr);
    long double forward_time = measure_latency(arr, size, output);

    initialize_reverse(arr);
    long double reverse_time = measure_latency(arr, size, output);

    initialize_random(arr);
    long double random_time = measure_latency(arr, size, output);

    output << size * static_cast<int>(sizeof(int)) << ";" << std::fixed
           << std::setprecision(4) << forward_time << ";" << reverse_time << ";"
           << random_time << "\n";
  }
}

int main() {

  std::ofstream out_file("out.csv");

  out_file << "Size;Forward(tact);Reverse(tact);Random(tact)\n";

  tieme_array_traversal(out_file);

  std::cout << "DONE\n";
  return 0;
}
