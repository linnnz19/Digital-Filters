#!/bin/bash

gcc -o simple_filter simple_filter.c -lm
./simple_filter 8 hL_8.txt hR_8.txt YL_8.txt YR_8.txt blue_giant_fragment.wav output_8.wav
./simple_filter 32 hL_32.txt hR_32.txt YL_32.txt YR_32.txt blue_giant_fragment.wav output_32.wav
./simple_filter 1024 hL_1024.txt hR_1024.txt YL_1024.txt YR_1024.txt blue_giant_fragment.wav output_1024.wav

#不知為何我無法使用python3編譯，這邊使用python
python show_data.py

mkdir Ascii_txt
mkdir Impulse_Response
mkdir Log_Spectrum

# Move 12 .txt files to Ascii_txt
mv "hL_8.txt" Ascii_txt
mv "hL_32.txt" Ascii_txt
mv "hL_1024.txt" Ascii_txt
mv "hR_8.txt" Ascii_txt
mv "hR_32.txt" Ascii_txt
mv "hR_1024.txt" Ascii_txt
mv "YL_8.txt" Ascii_txt
mv "YL_32.txt" Ascii_txt
mv "YL_1024.txt" Ascii_txt
mv "YR_8.txt" Ascii_txt
mv "YR_32.txt" Ascii_txt
mv "YR_1024.txt" Ascii_txt

# Move impulse response pngs to Impulse_Response
mv "left_channel_impulse_response_8.png" Impulse_Response
mv "left_channel_impulse_response_32.png" Impulse_Response
mv "left_channel_impulse_response_1024.png" Impulse_Response
mv "right_channel_impulse_response_8.png" Impulse_Response
mv "right_channel_impulse_response_32.png" Impulse_Response
mv "right_channel_impulse_response_1024.png" Impulse_Response

# Move spectrum pngs to Log_Spectrum
mv "left_channel_spectrum_8.png" Log_Spectrum
mv "left_channel_spectrum_32.png" Log_Spectrum
mv "left_channel_spectrum_1024.png" Log_Spectrum
mv "right_channel_spectrum_8.png" Log_Spectrum
mv "right_channel_spectrum_32.png" Log_Spectrum
mv "right_channel_spectrum_1024.png" Log_Spectrum


echo "All mission done!"