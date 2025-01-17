### 參數設定
```
#define N 16       // 簽名波形長度（16 個晶片）
#define MAX_ITER 1000 // 最大迭代次數
#define ALPHA 0.1   // 根升餘弦濾波器的滾降係數
#define Tc 1.0     // 晶片週期
#define p 0.2     // 馬可夫機率
#define oversampling 10     // 采樣率（每個晶片的分段數）
```
### 初始化簽名波形 c_k^(0)(t)
```
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length, int oversampling)
#### 輸入 ####
## k_z: c_k^(0)(t)        空陣列(陣列元素數: N*oversampling )
## s_k: 馬可夫序列          空陣列(陣列元素數: 2*N  )
## user_index:             使用者數(暫不考慮,設為0)
## length: N


#### 輸出 ####
## k_z: c_k^(0)(t)
```
>>### cos值計算
```
double psi(double t)
#### 輸出 **
cos(PI * t / Tc)

```
### 迭代(根升餘弦濾波器)
**1. 生成根生餘弦濾波器r(t)**
```
void generate_rrc_filter(double *filter, int length, double alpha, double Ts)
#### 輸入 ####
## filter              空陣列(陣列元素數: N*oversampling )
## length: N

#### 輸出 ####
## filter
```
**2. c_k與濾波器進行捲積**
```
void generate_ce_waveform(double *k_z, int length, int save_iters[], int num_saves, const char *output_file)
#### 輸入 ####
## k_z :原始(未迭代)c_k
## length: N
## save_iters[] :特定迭代次數{0,1,10,100,1000}
## num_saves : save_iters[]元素數
## output_file : 存日每筆數據,格式 {迭代次數,num_saves,c_k相位}

#### 輸出 ####
## output_file: save_iters[]元素數*N*oversampling = 5*16*10=800
```
>>### 振幅正規化
```
void normalize(double *k_z, int length)
#### 輸入 ####
## k_z : c_k                 (陣列元素數: N*oversampling )
## length :N
#### 輸出 ####
## k_z :(已正規化)
```
**輸出完成檔 : "output_file"** 應含有800筆數據
