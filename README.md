# Mini-Project---5-Digital-Filters
411086010 通訊三 林紘毅  

---

## Table of Contents

1. simple_filter.c 概念說明
2. show_data.py 概念說明
3. 繪製的 impulse responses & log spectrum
4. 綜合討論

### simple_filter.c 概念說明

這次的 simple_filter.c 程式碼是基於教授提供的半成品程式碼去做功能擴充的。

濾波後的.wav檔案分別經過了左聲道和右聲道的band-pass & band-stop filter，老師已於半成品程式碼裡提供我們band-pass filter的樣式，我們僅需要透過更改濾波範圍(FL、FH)即可得到我們想要得到的band-pass效果，

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






