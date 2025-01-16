#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>



#define N 16       // 簽名波形長度（16 個晶片）
#define MAX_ITER 1000
#define ALPHA 0.1   // 根升余弦濾波器的滾降係數

// 根升余弦濾波器生成
void generate_rrc_filter(double *filter, int length, double alpha, double Ts) {
    for (int i = 0; i < length; i++) {
        double t = (i - length / 2) / Ts; // 中心化時間
        if (fabs(t) < 1e-6) {
            filter[i] = 1.0; // 避免 t = 0 的分母為零
        } else if (fabs(fabs(t) - 1.0 / (2.0 * alpha)) < 1e-6) {
            filter[i] = (alpha / sqrt(2)) * ((1 + 2 / M_PI) * sin(M_PI / (4 * alpha)) +
                                             (1 - 2 / M_PI) * cos(M_PI / (4 * alpha)));
        } else {
            filter[i] = (sin(M_PI * t * (1 - alpha)) +
                         4 * alpha * t * cos(M_PI * t * (1 + alpha))) /
                        (M_PI * t * (1 - pow(4 * alpha * t, 2)));
        }
    }
}

// 使用馬可夫過程生成隨機序列
void initialize_waveform(double *k_z, int length, double p) {
    // 初始狀態：隨機選擇 +1 或 -1
    k_z[0] = (rand() % 2 ? 1 : -1);

    // 生成馬可夫序列
    for (int i = 1; i < length; i++) {
        double rand_prob = (double)rand() / RAND_MAX; // 隨機機率 [0, 1)
        if (k_z[i - 1] == 1) {
            k_z[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 的機率 p
        } else {
            k_z[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 的機率 p
        }
    }
}

// 幅值正規化
void normalize(double *k_z, int length) {
    for (int i = 0; i < length; i++) {
        k_z[i] = copysign(1.0, k_z[i]);
    }
}

// 半正弦波形 psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
        return cos(M_PI * t / Tc);
    }
    return 0.0;
}

// 初始化簽名波形 c_k^(0)(t)
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length, int oversampling) {
    int extended_length = 2 * length; // 2N
    double step = Tc / oversampling; // 時間步長
    int total_samples = extended_length * oversampling;

    for (int i = 0; i < total_samples; i++) {
        double t = i * step; // 當前時間
        double complex value = 0.0;

        // 根據公式 (6) 計算波形
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * M_PI / 2.0; // (-1)^k n 的相位偏移
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
            value += phase * s_k[n] * psi(t - n * Tc / 2.0);
        }

        k_z[i] = value; // 保存波形
    }
}

// 主流程，並保存特定迭代次數的數據
void generate_ce_waveform(double *k_z, int length, int save_iters[], int num_saves, const char *output_file) {
    double filter[length];
    generate_rrc_filter(filter, length, ALPHA, N); // 生成根升余弦濾波器

    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "Iteration,Time,Phase\n"); // 寫入表頭

    int save_index = 0;
    for (int iter = 0; iter <= MAX_ITER; iter++) {
        // 濾波（簡化的循環卷積）
        for (int i = 0; i < length; i++) {
            double temp = 0.0;
            for (int j = 0; j < length; j++) {
                temp += k_z[j] * filter[(i - j + length) % length];
            }
            k_z[i] = temp;
        }

        normalize(k_z, length); // 正規化幅值

        // 保存數據（僅在指定迭代次數）
        if (save_index < num_saves && iter == save_iters[save_index]) {
            for (int i = 0; i < length; i++) {
                fprintf(file, "%d,%f,%f\n", iter, i * 1.0, atan2(0, k_z[i])); // 計算相位
            }
            save_index++;
        }
    }

    fclose(file);
}

int main() {
    double k_z[N];
    double p = 0.2; // 馬可夫過程的轉移概率
    int save_iters[] = {0, 1, 10, 100, 1000}; // 要保存的迭代次數
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]);

    srand((unsigned int)time(NULL)); // 初始化隨機數種子
    initialize_waveform(k_z, N, p); // 初始化波形（馬可夫過程）

    generate_ce_waveform(k_z, N, save_iters, num_saves, "ce_waveform_data.csv"); // 保存數據

    printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
