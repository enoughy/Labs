#include <fstream>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <locale>
#include <vector>

int const offset = 1024 * 8 * 4; // смещение 128кб для int(4)

std::vector<int> create_array(int size_cache, int fragments_count) {
  int size = offset * fragments_count;
  std::vector<int> array(size);
  for (int j = 0; j < size; j += offset) {

    for (int i = 0; i < size_cache / fragments_count / sizeof(int); i++) {
      array[i + j] = (i + j + offset) % (size);
    }
  }
  return array;
}

const int K = 100;
long double travaling_time(int fragments_count, int size_cache) {
  std::vector<int> arr = create_array(size_cache, fragments_count);
  int size = arr.size();
  volatile int x = 0;
  auto start = __rdtsc();
  for (int i = 0; i < size * K; i++) {
    x = arr[x];
  }
  auto end = __rdtsc();
  long double avg_latency = (end - start) / (double)(size * K);

  return avg_latency;
}

void benchmark_to_csv(int size_cache, int min_fragments, int max_fragments,
                      const std::string &filename) {
  std::ofstream csv_file(filename);
  if (!csv_file.is_open()) {
    std::cerr << "Не удалось открыть файл для записи!" << std::endl;
    return;
  }
  std::locale::global(
      std::locale("ru_RU.utf8")); // Локаль для русского// языка
                                  // std::cout.imbue(std::locale("ru_RU.utf8"));
  //  Записываем заголовок CSV
  csv_file << "Fragments Count,Elapsed Time (ticks)\n";

  // Проводим измерения
  for (int fragments = min_fragments; fragments <= max_fragments; ++fragments) {
    long double time = travaling_time(fragments, size_cache);
    csv_file << fragments << "," << time << "\n";
    // std::cout << std::setprecision(6) << "Fragments: " << fragments
    //          << ", Time: " << time << " ticks" << std::endl;
  }

  csv_file.close();
  std::cout << "Результаты записаны в файл: " << filename << std::endl;
}

int main() {
  int size_cache = 64 * 1024; // Размер кэша (например, 64 КБ)
  int min_fragments = 1;      // Минимальное количество фрагментов
  int max_fragments = 32;     // Максимальное количество фрагментов
  std::string filename = "cache_benchmark.csv";

  benchmark_to_csv(size_cache, min_fragments, max_fragments, filename);

  return 0;
}
