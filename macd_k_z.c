#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>



#define N 16       // ñ�W�i�Ϊ��ס]16 �Ӵ����^
#define MAX_ITER 1000
#define ALPHA 0.1   // �ڤɧE���o�i�����u���Y��

// �ڤɧE���o�i���ͦ�
void generate_rrc_filter(double *filter, int length, double alpha, double Ts) {
    for (int i = 0; i < length; i++) {
        double t = (i - length / 2) / Ts; // ���ߤƮɶ�
        if (fabs(t) < 1e-6) {
            filter[i] = 1.0; // �קK t = 0 ���������s
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

// �ϥΰ��i�ҹL�{�ͦ��H���ǦC
void initialize_waveform(double *k_z, int length, double p) {
    // ��l���A�G�H����� +1 �� -1
    k_z[0] = (rand() % 2 ? 1 : -1);

    // �ͦ����i�ҧǦC
    for (int i = 1; i < length; i++) {
        double rand_prob = (double)rand() / RAND_MAX; // �H�����v [0, 1)
        if (k_z[i - 1] == 1) {
            k_z[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 �����v p
        } else {
            k_z[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 �����v p
        }
    }
}

// �T�ȥ��W��
void normalize(double *k_z, int length) {
    for (int i = 0; i < length; i++) {
        k_z[i] = copysign(1.0, k_z[i]);
    }
}

// �b�����i�� psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
        return cos(M_PI * t / Tc);
    }
    return 0.0;
}

// ��l��ñ�W�i�� c_k^(0)(t)
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length, int oversampling) {
    int extended_length = 2 * length; // 2N
    double step = Tc / oversampling; // �ɶ��B��
    int total_samples = extended_length * oversampling;

    for (int i = 0; i < total_samples; i++) {
        double t = i * step; // ��e�ɶ�
        double complex value = 0.0;

        // �ھڤ��� (6) �p��i��
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * M_PI / 2.0; // (-1)^k n ���ۦ찾��
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
            value += phase * s_k[n] * psi(t - n * Tc / 2.0);
        }

        k_z[i] = value; // �O�s�i��
    }
}

// �D�y�{�A�ëO�s�S�w���N���ƪ��ƾ�
void generate_ce_waveform(double *k_z, int length, int save_iters[], int num_saves, const char *output_file) {
    double filter[length];
    generate_rrc_filter(filter, length, ALPHA, N); // �ͦ��ڤɧE���o�i��

    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "Iteration,Time,Phase\n"); // �g�J���Y

    int save_index = 0;
    for (int iter = 0; iter <= MAX_ITER; iter++) {
        // �o�i�]²�ƪ��`�����n�^
        for (int i = 0; i < length; i++) {
            double temp = 0.0;
            for (int j = 0; j < length; j++) {
                temp += k_z[j] * filter[(i - j + length) % length];
            }
            k_z[i] = temp;
        }

        normalize(k_z, length); // ���W�ƴT��

        // �O�s�ƾڡ]�Ȧb���w���N���ơ^
        if (save_index < num_saves && iter == save_iters[save_index]) {
            for (int i = 0; i < length; i++) {
                fprintf(file, "%d,%f,%f\n", iter, i * 1.0, atan2(0, k_z[i])); // �p��ۦ�
            }
            save_index++;
        }
    }

    fclose(file);
}

int main() {
    double k_z[N];
    double p = 0.2; // ���i�ҹL�{���ಾ���v
    int save_iters[] = {0, 1, 10, 100, 1000}; // �n�O�s�����N����
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]);

    srand((unsigned int)time(NULL)); // ��l���H���ƺؤl
    initialize_waveform(k_z, N, p); // ��l�ƪi�Ρ]���i�ҹL�{�^

    generate_ce_waveform(k_z, N, save_iters, num_saves, "ce_waveform_data.csv"); // �O�s�ƾ�

    printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
