# Mini-Project---5-Digital-Filters
411086010 通訊三 林紘毅  

---

## Table of Contents

#### 1. *simple_filter.c* 概念說明
#### 2. *show_data.py* 概念說明
#### 3. 繪製的 impulse responses & log spectrum
#### 4. 綜合討論

---

### 1. *simple_filter.c* 概念說明

這次的 *simple_filter.c* 程式碼是基於教授提供的半成品程式碼去做功能擴充的。

濾波後的 *.wav* 檔案分別經過了左聲道和右聲道的band-pass & band-stop filter，老師已於半成品程式碼裡提供我們band-pass filter的樣式，我們僅需要透過更改濾波範圍(FL、FH)即可得到我們想要得到的band-pass效果，

至於band-stop filter我則是透過他近似於high-pass filter減掉low-pass filter來去實現的，high-pass與low-pass程式碼如下所示：
```c
float low_pass(int m, int n){
	float wc = 2 * PI * FL / FS;
	if(n == m) {// L'Hopital's Rule
		return wc / PI;
	}
	else {
		return sinf(wc * ((float)(n - m)))/PI/((float)(n - m)) * hamming(2 * m + 1, n);
	}
}

float high_pass(int m, int n) {
    float wc = 2 * PI * FH / FS;
    if (n == m) {
        return 1.0 - wc / PI;
    } else {
        return (sinf(PI * ((float)(n - m))) - sinf(wc * ((float)(n - m)))) / PI / ((float)(n - m)) * hamming(2 * m + 1, n);
    }
}
```

透過高通 & 低通濾波器的差去達成近似band-stop filter的效果，並且在主函式中透過for迴圈把濾波後的結果寫入音檔，如下：
```c
for(n=0;n<(2*M+1);n++) {
		h_L[n] = band_pass(M, n);
		fprintf(file_hL, "%.15e\n", h_L[n]); // Write to hL.txt

        h_R[n] = high_pass(M, n) - low_pass(M, n);
		fprintf(file_hR, "%.15e\n", h_R[n]); // Write to hR.txt
	}
```

另外，在這個迴圈我也順便把hL & hR的係數寫入矩陣裡，以便稍後我們產生出impulse response的相關係數。

---

接下來的任務是計算出濾波後的音檔之log spectrum係數，並且把他存到YL & YR的文字檔裡。

要得到spectrum(頻譜圖)，我們需要做的任務就是**傅立葉轉換**。
至於是使用何種方式來進行**DFT**的呢?我使用的是window method，並且參考老師程式碼裡已提供的hamming window函式，如下所示：

```c
/* hamming: for n=0,1,2,...N, length of N+1 */
float hamming(int N, int n)
{
	return 0.54 - 0.46 * cosf(2*PI*((float)(n))/((float)N));
}
```
透過應用hamming window，現在我們要對一段特定時間內的sample點進行轉換，題目提供給我們的是20.060 秒到 20.085 秒之間，也就代表N值應該為針對第**962880~964008**個取樣點來做傅立葉轉換(秒數乘上sample rate)。

我參考了網路上的演算法
[DFT演算法網址](https://www.geeksforgeeks.org/discrete-fourier-transform-and-its-inverse-using-c/amp/)

做出了一個我的版本之DFT函式，如下所示：

```c
// 傅立葉轉換函數
void LogSpectrum_DFT(short int *xn, int len, FILE *file_Y) {
    int k, n;
    float Xr[N_FFT];
    float Xi[N_FFT];

    for (k = 0; k < N_FFT; k++) {
        Xr[k] = 0;
        Xi[k] = 0;

        for (n = START_SAMPLE; n <= END_SAMPLE; n++) {
            Xr[k] += xn[n] * hamming(N_FFT - 1, n - START_SAMPLE) * cosf(2 * PI * k * (n - START_SAMPLE) / N_FFT);
            Xi[k] -= xn[n] * hamming(N_FFT - 1, n - START_SAMPLE) * sinf(2 * PI * k * (n - START_SAMPLE) / N_FFT);
        }

        // 取log後乘以20
        float magnitude = 20 * log10f(sqrtf(Xr[k] * Xr[k] + Xi[k] * Xi[k]));
        magnitude = fabsf(magnitude); // 取絕對值
        fprintf(file_Y, "%.15e\n", magnitude);
    }
}
```
如此一來，在主函式裡呼叫LogSpectrum_DFT函式，並且分別輸入濾波後的左聲道和右聲道到這個函示裡，便可以成功的把傅立葉轉換後再轉換成log spectrum的值存到YL & YR文字檔裡了。

最後加上一些除錯的程式碼以及印出當前工作狀態在每一項子任務程式，*simple_filter.c* 的工作就大功告成。

---

### 2. *show_data.py* 概念說明  

*show_data.py* 的主要功能為讀取 *simple_filter.c* 所產生的各個存取音訊檔案的文字檔，並且繪製相對應的impulse response & log spectrum圖形。這個python程式主要可以分成兩個部分：

* 脈衝響應圖表(impulse response)：

   *plot_impulse_response* 函式讀取hR & hL的各種不同M值的檔案，並且使用 python 的 matplotlib 之 stem 指令繪製產生這些 impulse responses，程式碼如下：

```python
def plot_impulse_response(filename, channel_name, scenario):
    # 讀取檔案並繪製脈衝響應圖表
    with open(filename, 'r') as file:
        lines = file.readlines()
        impulse_response = [float(line.strip()) for line in lines]

    plt.figure(figsize=(8, 6))
    plt.stem(range(len(impulse_response)), impulse_response, basefmt=' ', markerfmt=' ')
    plt.title(f"Impulse Response for {channel_name} channel (M = {scenario})")
    plt.xlabel("Sample Index")
    plt.ylabel("Amplitude")
    plt.grid(True)
    
    plt.savefig(f"{channel_name.lower()}_channel_impulse_response_{scenario}.png")
    plt.close()

file_scenarios = ['8', '32', '1024']
channel_names = ['Left', 'Right']

for scenario in file_scenarios:
    for channel_name in channel_names:
        filename = f'h{channel_name[0]}_{scenario}.txt'
        plot_impulse_response(filename, channel_name, scenario)
```
在以上我設定的規範之下，這段python腳本可以成功產生出每一個不同M值下的左右聲道脈衝響應圖。

* 頻譜圖(log spectrum)：

   程式緊接著處理左右聲道的頻譜資料，它讀取YR & YL的各種不同M值的檔案，並且使用 FFT（快速傅立葉變換）計算了左右聲道的頻譜，然後繪製了頻譜圖。程式碼如下：

```python
for scenario in file_scenarios:
    file_YL_name = f'YL_{scenario}.txt'
    file_YR_name = f'YR_{scenario}.txt'

    # 讀取左右聲道資料並繪製頻譜圖
    with open(file_YL_name, 'r') as file_YL, open(file_YR_name, 'r') as file_YR:
        data_YL_dB = [float(line.strip()) for line in file_YL]
        data_YR_dB = [float(line.strip()) for line in file_YR]

    num_samples_YL = len(data_YL_dB)
    num_samples_YR = len(data_YR_dB)
    num_samples = min(num_samples_YL, num_samples_YR)

    sampling_rate = 48000
    frequencies = np.fft.fftfreq(num_samples, d=1/sampling_rate)[:num_samples//2]
    mask = (frequencies >= 0) & (frequencies <= 20000)
    frequencies = frequencies[mask]
    data_YL_dB = np.array(data_YL_dB[:num_samples//2])[mask]
    data_YR_dB = np.array(data_YR_dB[:num_samples//2])[mask]

    # 繪製左聲道的頻譜圖
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies, data_YL_dB, 'b')
    plt.title(f'Log Spectrum of Left Channel (M = {scenario})')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (dB)')
    plt.yscale('linear')
    plt.xlim(0, 20000)

    plt.savefig(f'left_channel_spectrum_{scenario}.png')
    plt.close()

    # 繪製右聲道的頻譜圖
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies, data_YR_dB, 'b')
    plt.title(f'Log Spectrum of Right Channel (M = {scenario})')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (dB)')
    plt.yscale('linear')
    plt.xlim(0, 20000)

    plt.savefig(f'right_channel_spectrum_{scenario}.png')
    plt.close()
    
print("pngs of log spectrums generating complete!")
```

將以上兩段程式碼結合在一起，便能成功的將所要求的impulse response和log spectrum都繪製出來，並且生成png檔案儲存在當前工作區。

---
### 3. 繪製的 impulse responses & log spectrum 
* impulse responses
* 
  **左聲道：**
  ![left_channel_impulse_response_8](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/be12542d-f46e-43a6-b291-37a2e6c0b639)
  <p align="center">圖一：M=8之左聲道脈衝響應</p>  
  
  ![left_channel_impulse_response_32](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/29ba663f-2f8b-4749-a8ca-82585b929ffb)
  <p align="center">圖二：M=32之左聲道脈衝響應</p>

  ![left_channel_impulse_response_1024](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/a41d73fc-a38f-4c57-be1b-225cf5499c17)
  <p align="center">圖三：M=1024之左聲道脈衝響應</p>

  **右聲道：**
  ![right_channel_impulse_response_8](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/54d203ab-10b5-4da1-ad44-e2aef2d6b2d6)
  <p align="center">圖四：M=8之右聲道脈衝響應</p>

  ![right_channel_impulse_response_32](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/e8ac53ce-b14e-4b3c-984b-4fe7b9054508)
  <p align="center">圖五：M=32之右聲道脈衝響應</p>

  ![right_channel_impulse_response_1024](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/d9960401-d1c9-4372-8218-bd9d23acbfb1)
  <p align="center">圖六：M=1024之右聲道脈衝響應</p>  

  * log spectrums
  * 
    **左聲道：**
    ![left_channel_spectrum_8](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/aed5a10e-ce44-4295-8291-950c0c613ff8)
    <p align="center">圖七：M=8之左聲道頻譜圖</p>

    ![left_channel_spectrum_32](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/aff609f1-d47f-4e68-adbe-df506ce875a9)
    <p align="center">圖八：M=32之左聲道頻譜圖</p>

    ![left_channel_spectrum_1024](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/18defb22-f57e-4ee3-994c-23978fadcb73)
    <p align="center">圖九：M=1024之左聲道頻譜圖</p>  

    **右聲道：**
    ![right_channel_spectrum_8](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/875f0f3d-cafd-4479-99af-cb001940ef0a)
    <p align="center">圖十：M=8之右聲道頻譜圖</p>

    ![right_channel_spectrum_32](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/661f9cf7-5e43-4fa9-b34f-a1cf4a260b79)
    <p align="center">圖十一：M=32之右聲道頻譜圖</p>

    ![right_channel_spectrum_1024](https://github.com/linnnz19/Mini-Project---5-Digital-Filters/assets/128024684/1354ab0e-1572-4427-b2b6-044583c82d92)
    <p align="center">圖十二：M=1024之右聲道頻譜圖</p>


    

    


    



    


    
  


  






