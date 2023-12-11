# Mini-Project---5-Digital-Filters
411086010 通訊三 林紘毅  

---

## Table of Contents

1. simple_filter.c 概念說明
2. show_data.py 概念說明
3. 繪製的 impulse responses & log spectrum
4. 綜合討論

### simple_filter.c 概念說明

這次的 *simple_filter.c* 程式碼是基於教授提供的半成品程式碼去做功能擴充的。

濾波後的 *.wav* 檔案分別經過了左聲道和右聲道的band-pass & band-stop filter，老師已於半成品程式碼裡提供我們band-pass filter的樣式，我們僅需要透過更改濾波範圍(FL、FH)即可得到我們想要得到的band-pass效果，

至於band-stop filter我則是透過他近似於high-pass filter減掉low-pass filter來去實現的，high-pass與low-pass程式碼如下所示：
```C=
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

透過高通 & 低通濾波器的差去達成近似band-stop的效果，並且在主函式中透過for迴圈把濾波後的結果寫入音檔，如下：
```C=
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
至於是使用何種方式來進行**DFT**的呢?我使用的是window method，並且參考老師程式碼裡已提供的hamming window函示，如下所示：

```C=
/* hamming: for n=0,1,2,...N, length of N+1 */
float hamming(int N, int n)
{
	return 0.54 - 0.46 * cosf(2*PI*((float)(n))/((float)N));
}
```
透過應用hamming window，現在我們要對一段特定時間內的sample點進行轉換，題目提供給我們的是20.060 秒到 20.085 秒之間，也就代表N值應該為針對**962880~964008**的取樣點來做傅立葉轉換(秒數乘上sample rate)。

我參考了網路上的演算法
[DFT演算法網址](https://www.geeksforgeeks.org/discrete-fourier-transform-and-its-inverse-using-c/amp/)

做出了一個我的版本之DFT函示，如下所示：

```C=
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
如此一來，在主函式裡呼叫LogSpectrum_DFT函示，並且分別輸入濾波後的左聲道和右聲道到這個函示裡，便可以成功的把傅立葉轉換後再轉換成log spectrum的值存到YL & YR文字檔裡了。

最後加上一些除錯的程式碼以及印出當前工作狀態在每一項子任務程式，*simple_filter.c* 的工作就大功告成。






