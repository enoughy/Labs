#include <chrono>
#include <cstdlib>
#include <immintrin.h>
#include <iostream>
#include <stdlib.h>

float *allocate_matrix(int N) {
  float *matrix;
  size_t size = N * N * sizeof(float);

  if (posix_memalign((void **)&matrix, 32, size) != 0) {
    perror("posix_memalign failed");
    return NULL;
  }

  return matrix;
}
void check_alignment(const void *ptr, const char *name) {
  if ((uintptr_t)ptr % 32 != 0) {
    printf("Pointer %s is not aligned to 32 bytes\n", name);
  }
}

void matrix_mult(int M, int N, int K, const float *A, const float *B,
                 float *C) {
  for (int i = 0; i < M; ++i) {
    float *c = C + i * N;

    // Обнуляем строку C
    for (int j = 0; j < N; ++j)
      c[j] = 0;

    for (int k = 0; k < K; ++k) {
      const float *b = B + k * N;
      float a = A[i * K + k];

      // Загружаем элемент a в SIMD-регистр
      __m256 a_vec = _mm256_set1_ps(a);

      for (int j = 8 - (k * N % 8); j <= N - 8; j += 8) {
        // Загружаем 8 элементов из b и c
        __m256 b_vec = _mm256_load_ps(&b[j]);
        __m256 c_vec = _mm256_load_ps(&c[j]);

        // Выполняем FMA: c[j:j+8] += a * b[j:j+8]
        c_vec = _mm256_fmadd_ps(a_vec, b_vec, c_vec);

        // Сохраняем результат обратно в c
        _mm256_store_ps(&c[j], c_vec);
      }

      for (int j = 0; j < 8 - (k * N % 8); ++j) {
        c[j] += a * b[j];
      }

      for (int j = (N / 8) * 8; j < N; ++j) {
        c[j] += a * b[j];
      }
    }
  }
}

void matrix_summ(int N, int K, float *B, float *A) {
  int size = N * N; // Размер матрицы (число элементов)
  int i = 0;

  for (; i <= size - 8; i += 8) {
    __m256 vecA = _mm256_load_ps(&A[i]);       // Загружаем 8 элементов из A
    __m256 vecB = _mm256_load_ps(&B[i]);       // Загружаем 8 элементов из B
    __m256 vecSum = _mm256_add_ps(vecA, vecB); // Складываем A и B
    _mm256_store_ps(&A[i], vecSum); // Сохраняем результат обратно в A
  }

  for (; i < size; i++) {
    A[i] += B[i];
  }
}

void printMatrix(int N, float *matrix) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f ", matrix[N * i + j]);
    }
    printf("\n");
  }
}

float norm_1(int N, float *A) {
  float max = 0;
  float max_temp = 0;
  for (int j = 0; j < N; j++) {
    max_temp = 0;
    for (int i = 0; i < N; i++) {
      max_temp += A[i * N + j];
    }
    if (max < max_temp) {
      max = max_temp;
    }
  }
  return max;
}

float norm_inf(int N, float *A) {
  float max = 0;
  float max_temp = 0;
  for (int i = 0; i < N; i++) {
    max_temp = 0;
    for (int j = 0; j < N; j++) {
      max_temp += A[i * N + j];
    }
    if (max < max_temp) {
      max = max_temp;
    }
  }
  return max;
}

void scalar_multiply(int N, float a, float *A) {
  for (int i = 0; i < N * N; ++i) {
    A[i] *= a;
  }
}

void transpose_matrix(int N, const float *input, float *output) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      output[j * N + i] = input[i * N + j];
    }
  }
}

void calculate_b(int N, float *A, float *B) {
  float n_1 = norm_1(N, A);
  float n_inf = norm_inf(N, A);
  transpose_matrix(N, A, B);
  scalar_multiply(N, 1 / (n_1 * n_inf), B);
}

void create_identity_matrix(int N, float *matrix) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      matrix[i * N + j] = (i == j) ? 1.0f : 0.0f;
    }
  }
}

void calculate_r(int N, float *A, float *B, float *R) {
  // float *I = (float *)calloc(N * N, sizeof(float));
  float *I = allocate_matrix(N);
  create_identity_matrix(N, I);
  matrix_mult(N, N, N, B, A, R);
  scalar_multiply(N, -1, R);
  matrix_summ(N, N, I, R);
  free(I);
}

float *reverse_matrix(float *A, int N, int M) {
  // float *B = (float *)calloc(N * N, sizeof(float));
  // float *R = (float *)calloc(N * N, sizeof(float));
  // float *R_degree = (float *)calloc(N * N, sizeof(float));
  // float *temp = (float *)calloc(N * N, sizeof(float));
  // float *R_temp = (float *)calloc(N * N, sizeof(float));
  float *B = allocate_matrix(N);
  float *R = allocate_matrix(N);
  float *R_degree = allocate_matrix(N);
  float *temp = allocate_matrix(N);
  create_identity_matrix(N, temp);
  calculate_b(N, A, B);
  calculate_r(N, A, B, R);
  matrix_mult(N, N, N, temp, R, R_degree);
  matrix_summ(N, N, R, temp);
  float *R_temp = allocate_matrix(N);
  for (int i = 0; i < M; i++) {
    // R_temp = (float *)calloc(N * N, sizeof(float));
    matrix_mult(N, N, N, R, R_degree, R_temp);
    auto R_degree_ptr = R_degree;
    R_degree = R_temp;
    R_temp = R_degree_ptr;
    matrix_summ(N, N, R_degree, temp);
  }
  // float *revers_A = (float *)calloc(N * N, sizeof(float));
  float *revers_A = allocate_matrix(N);
  free(R_temp);
  matrix_mult(N, N, N, temp, B, revers_A);
  free(R_degree);
  free(temp);
  free(B);
  free(R);
  return revers_A;
}

float *createRandomSquareMatrix(int size) {
  // Выделяем память для матрицы
  // float *matrix = (float *)malloc(size * size * sizeof(float));
  float *matrix = allocate_matrix(size);
  if (matrix == NULL) {
    printf("Error: Memory allocation failed!\n");
    return NULL;
  }

  // Инициализируем генератор случайных чисел
  srand(time(NULL));

  // Заполняем матрицу случайными значениями
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      matrix[i * size + j] =
          (float)rand() / RAND_MAX * 100.0f; // Случайные числа от 0 до 100
    }
  }

  return matrix;
}

int main() {
  int N = 2024;
  float *A = createRandomSquareMatrix(N);
  auto start = std::chrono::high_resolution_clock::now();
  float *T = reverse_matrix(A, N, 10);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
  // printMatrix(N, T);
}
