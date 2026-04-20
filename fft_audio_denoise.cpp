#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Complex {
    double real;
    double imag;

    Complex(double real_value = 0.0, double imag_value = 0.0)
        : real(real_value), imag(imag_value) {}
};

Complex operator+(const Complex& lhs, const Complex& rhs) {
    return Complex(lhs.real + rhs.real, lhs.imag + rhs.imag);
}

Complex operator-(const Complex& lhs, const Complex& rhs) {
    return Complex(lhs.real - rhs.real, lhs.imag - rhs.imag);
}

Complex operator*(const Complex& lhs, const Complex& rhs) {
    return Complex(lhs.real * rhs.real - lhs.imag * rhs.imag,
                   lhs.real * rhs.imag + lhs.imag * rhs.real);
}

Complex operator*(const Complex& value, double scale) {
    return Complex(value.real * scale, value.imag * scale);
}

void fft(vector<Complex>& values, bool invert) {
    // 先做 bit-reversal 重排，后面的蝶形迭代才能原地进行。
    const int n = static_cast<int>(values.size());
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            swap(values[i], values[j]);
        }
    }

    // 逐层合并长度为 2、4、8... 的子问题，直到覆盖整个序列。
    for (int len = 2; len <= n; len <<= 1) {
        const double angle = 2.0 * kPi / static_cast<double>(len) * (invert ? -1.0 : 1.0);
        const Complex wn(cos(angle), sin(angle));
        for (int start = 0; start < n; start += len) {
            Complex w(1.0, 0.0);
            const int half = len >> 1;
            for (int offset = 0; offset < half; ++offset) {
                const Complex even = values[start + offset];
                const Complex odd = w * values[start + offset + half];
                values[start + offset] = even + odd;
                values[start + offset + half] = even - odd;
                w = w * wn;
            }
        }
    }

    if (invert) {
        for (Complex& value : values) {
            value.real /= n;
            value.imag /= n;
        }
    }
}

size_t next_power_of_two(size_t value) {
    // Radix-2 FFT 需要长度为 2 的幂，不足部分后续会补零。
    size_t power = 1;
    while (power < value) {
        power <<= 1;
    }
    return power;
}

uint16_t read_u16(ifstream& input) {
    unsigned char bytes[2];
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    if (!input) {
        throw runtime_error("Failed to read uint16.");
    }
    return static_cast<uint16_t>(bytes[0] | (bytes[1] << 8));
}

uint32_t read_u32(ifstream& input) {
    unsigned char bytes[4];
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    if (!input) {
        throw runtime_error("Failed to read uint32.");
    }
    return static_cast<uint32_t>(bytes[0]) |
           (static_cast<uint32_t>(bytes[1]) << 8) |
           (static_cast<uint32_t>(bytes[2]) << 16) |
           (static_cast<uint32_t>(bytes[3]) << 24);
}

void write_u16(ofstream& output, uint16_t value) {
    const unsigned char bytes[2] = {
        static_cast<unsigned char>(value & 0xFF),
        static_cast<unsigned char>((value >> 8) & 0xFF),
    };
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void write_u32(ofstream& output, uint32_t value) {
    const unsigned char bytes[4] = {
        static_cast<unsigned char>(value & 0xFF),
        static_cast<unsigned char>((value >> 8) & 0xFF),
        static_cast<unsigned char>((value >> 16) & 0xFF),
        static_cast<unsigned char>((value >> 24) & 0xFF),
    };
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

struct WavData {
    uint16_t channels = 0;
    uint32_t sample_rate = 0;
    uint16_t bits_per_sample = 0;
    vector<vector<double>> samples;
};

WavData read_wav(const string& path) {
    ifstream input(path, ios::binary);
    if (!input) {
        throw runtime_error("Cannot open input file: " + path);
    }

    char riff[4];
    input.read(riff, sizeof(riff));
    if (string(riff, 4) != "RIFF") {
        throw runtime_error("Input is not a RIFF file.");
    }
    (void)read_u32(input);
    char wave[4];
    input.read(wave, sizeof(wave));
    if (string(wave, 4) != "WAVE") {
        throw runtime_error("Input is not a WAVE file.");
    }

    uint16_t audio_format = 0;
    uint16_t channels = 0;
    uint32_t sample_rate = 0;
    uint16_t bits_per_sample = 0;
    vector<unsigned char> raw_data;

    // 按 chunk 遍历 WAV，只提取 fmt/data 两类我们真正需要的信息。
    while (input.peek() != EOF) {
        char chunk_id[4];
        input.read(chunk_id, sizeof(chunk_id));
        if (!input) {
            break;
        }
        const uint32_t chunk_size = read_u32(input);
        const string chunk_name(chunk_id, 4);

        if (chunk_name == "fmt ") {
            audio_format = read_u16(input);
            channels = read_u16(input);
            sample_rate = read_u32(input);
            (void)read_u32(input);
            (void)read_u16(input);
            bits_per_sample = read_u16(input);
            const uint32_t remaining = chunk_size > 16 ? chunk_size - 16 : 0;
            if (remaining > 0) {
                input.ignore(remaining);
            }
        } else if (chunk_name == "data") {
            raw_data.resize(chunk_size);
            input.read(reinterpret_cast<char*>(raw_data.data()), static_cast<streamsize>(chunk_size));
            if (!input) {
                throw runtime_error("Failed to read WAV data chunk.");
            }
        } else {
            input.ignore(chunk_size);
        }

        if (chunk_size & 1U) {
            input.ignore(1);
        }
    }

    if (audio_format != 1) {
        throw runtime_error("Only PCM WAV is supported.");
    }
    if (channels == 0 || sample_rate == 0 || bits_per_sample != 16) {
        throw runtime_error("Only 16-bit PCM WAV with valid channel/sample-rate metadata is supported.");
    }
    if (raw_data.empty()) {
        throw runtime_error("WAV file has no sample data.");
    }

    const size_t bytes_per_sample = bits_per_sample / 8;
    const size_t frame_size = channels * bytes_per_sample;
    if (raw_data.size() % frame_size != 0) {
        throw runtime_error("Corrupted WAV data chunk.");
    }

    const size_t frame_count = raw_data.size() / frame_size;
    WavData wav;
    wav.channels = channels;
    wav.sample_rate = sample_rate;
    wav.bits_per_sample = bits_per_sample;
    wav.samples.assign(channels, vector<double>(frame_count, 0.0));

    // WAV 数据是交错存储的，这里拆成按声道组织的浮点样本，范围归一化到 [-1, 1]。
    for (size_t frame = 0; frame < frame_count; ++frame) {
        for (uint16_t channel = 0; channel < channels; ++channel) {
            const size_t offset = frame * frame_size + channel * bytes_per_sample;
            const int16_t sample = static_cast<int16_t>(
                static_cast<uint16_t>(raw_data[offset]) |
                (static_cast<uint16_t>(raw_data[offset + 1]) << 8));
            wav.samples[channel][frame] = max(-1.0, min(1.0, static_cast<double>(sample) / 32768.0));
        }
    }

    return wav;
}

void write_wav(const string& path, const WavData& wav) {
    if (wav.samples.empty()) {
        throw runtime_error("No audio data to write.");
    }

    const size_t frame_count = wav.samples.front().size();
    const uint16_t block_align = wav.channels * (wav.bits_per_sample / 8);
    const uint32_t byte_rate = wav.sample_rate * block_align;
    const uint32_t data_bytes = static_cast<uint32_t>(frame_count * block_align);

    ofstream output(path, ios::binary);
    if (!output) {
        throw runtime_error("Cannot open output file: " + path);
    }

    // 写出一个最小可用的 PCM WAV 头，再顺序写回采样数据。
    output.write("RIFF", 4);
    write_u32(output, 36 + data_bytes);
    output.write("WAVE", 4);

    output.write("fmt ", 4);
    write_u32(output, 16);
    write_u16(output, 1);
    write_u16(output, wav.channels);
    write_u32(output, wav.sample_rate);
    write_u32(output, byte_rate);
    write_u16(output, block_align);
    write_u16(output, wav.bits_per_sample);

    output.write("data", 4);
    write_u32(output, data_bytes);

    for (size_t frame = 0; frame < frame_count; ++frame) {
        for (uint16_t channel = 0; channel < wav.channels; ++channel) {
            const double scaled = max(-1.0, min(1.0, wav.samples[channel][frame])) * 32767.0;
            const int16_t sample = static_cast<int16_t>(llround(scaled));
            write_u16(output, static_cast<uint16_t>(sample));
        }
    }
}

double low_pass_weight(double frequency_hz, double cutoff_hz, double transition_hz) {
    // 用余弦过渡带代替硬截断，减少时域振铃。
    if (frequency_hz <= cutoff_hz) {
        return 1.0;
    }
    if (transition_hz <= 0.0 || frequency_hz >= cutoff_hz + transition_hz) {
        return 0.0;
    }
    const double ratio = (frequency_hz - cutoff_hz) / transition_hz;
    return 0.5 * (1.0 + cos(kPi * ratio));
}

void denoise_channel(vector<double>& channel, uint32_t sample_rate, double cutoff_hz, double transition_hz) {
    const size_t original_size = channel.size();
    const size_t fft_size = next_power_of_two(max<size_t>(1, original_size));

    // 把时域实信号搬到复数数组，未填满部分保持 0 作为零填充。
    vector<Complex> spectrum(fft_size);
    for (size_t i = 0; i < original_size; ++i) {
        spectrum[i].real = channel[i];
    }

    fft(spectrum, false);

    // 对正负频率对应的谱系数统一施加低通权重，保留共轭对称性。
    const double bin_hz = static_cast<double>(sample_rate) / static_cast<double>(fft_size);
    for (size_t i = 0; i < fft_size; ++i) {
        const double frequency_hz = bin_hz * static_cast<double>(min(i, fft_size - i));
        const double weight = low_pass_weight(frequency_hz, cutoff_hz, transition_hz);
        spectrum[i] = spectrum[i] * weight;
    }

    fft(spectrum, true);

    for (size_t i = 0; i < original_size; ++i) {
        channel[i] = max(-1.0, min(1.0, spectrum[i].real));
    }
}

void print_usage(const char* program) {
    cerr << "Usage: " << program
         << " <input.wav> <output.wav> <cutoff_hz> [transition_hz]\n";
    cerr << "Example: " << program
         << " noisy.wav clean.wav 4500 800\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        print_usage(argv[0]);
        return 1;
    }

    try {
        const string input_path = argv[1];
        const string output_path = argv[2];
        const double cutoff_hz = stod(argv[3]);
        const double transition_hz = argc == 5 ? stod(argv[4]) : max(200.0, cutoff_hz * 0.1);

        // 截止频率必须落在有效范围内，否则滤波器没有物理意义。
        if (cutoff_hz <= 0.0) {
            throw runtime_error("cutoff_hz must be positive.");
        }

        WavData wav = read_wav(input_path);
        const double nyquist = static_cast<double>(wav.sample_rate) / 2.0;
        if (cutoff_hz >= nyquist) {
            throw runtime_error("cutoff_hz must be smaller than the Nyquist frequency.");
        }

        // 多声道分别处理，最后再重新封装回同一个 WAV 文件。
        for (vector<double>& channel : wav.samples) {
            denoise_channel(channel, wav.sample_rate, cutoff_hz, transition_hz);
        }

        write_wav(output_path, wav);

        cout << "Processed " << input_path << " -> " << output_path
             << " with cutoff " << cutoff_hz << " Hz";
        if (argc == 5) {
            cout << " and transition " << transition_hz << " Hz";
        }
        cout << ".\n";
    } catch (const exception& error) {
        cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
