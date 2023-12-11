import numpy as np
import matplotlib.pyplot as plt

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
        
print("pngs of impulse responses generating complete!")

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
