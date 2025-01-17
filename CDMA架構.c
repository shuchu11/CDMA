### �ѼƳ]�w
```
#define N 16       // ñ�W�i�Ϊ��ס]16 �Ӵ����^
#define MAX_ITER 1000 // �̤j���N����
#define ALPHA 0.1   // �ڤɾl���o�i�����u���Y��
#define Tc 1.0     // �����g��
#define p 0.2     // ���i�Ҿ��v
#define oversampling 10     // ���˲v�]�C�Ӵ��������q�ơ^
```
### ��l��ñ�W�i�� c_k^(0)(t)
```
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length, int oversampling)
#### ��J ####
## k_z: c_k^(0)(t)        �Ű}�C(�}�C������: N*oversampling )
## s_k: ���i�ҧǦC          �Ű}�C(�}�C������: 2*N  )
## user_index:             �ϥΪ̼�(�Ȥ��Ҽ{,�]��0)
## length: N


#### ��X ####
## k_z: c_k^(0)(t)
```
>>### cos�ȭp��
```
double psi(double t)
#### ��X **
cos(PI * t / Tc)

```
### ���N(�ڤɾl���o�i��)
**1. �ͦ��ڥ;l���o�i��r(t)**
```
void generate_rrc_filter(double *filter, int length, double alpha, double Ts)
#### ��J ####
## filter              �Ű}�C(�}�C������: N*oversampling )
## length: N

#### ��X ####
## filter
```
**2. c_k�P�o�i���i�汲�n**
```
void generate_ce_waveform(double *k_z, int length, int save_iters[], int num_saves, const char *output_file)
#### ��J ####
## k_z :��l(�����N)c_k
## length: N
## save_iters[] :�S�w���N����{0,1,10,100,1000}
## num_saves : save_iters[]������
## output_file : �s��C���ƾ�,�榡 {���N����,num_saves,c_k�ۦ�}

#### ��X ####
## output_file: save_iters[]������*N*oversampling = 5*16*10=800
```
>>### ���T���W��
```
void normalize(double *k_z, int length)
#### ��J ####
## k_z : c_k                 (�}�C������: N*oversampling )
## length :N
#### ��X ####
## k_z :(�w���W��)
```
**��X������ : "output_file"** ���t��800���ƾ�
