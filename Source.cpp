// This code is for a blog post on FIRs: https://blog.demofox.org/2020/01/05/fir-audio-data-filters/

#include <stdio.h>
#include <vector>
#include <complex>

static const size_t c_numFrequencies = 100;    // for frequency / phase response
static const size_t c_oscillatorSamples = 100; // for showing oscillator output

static const float c_pi = 3.14159265359f;

typedef std::complex<float> complex;

float DegreesToRadians(float degrees)
{
    return degrees * c_pi / 180.0f;
}

complex Z(int delay, float angle)
{
    return std::polar(1.0f, float(delay)*angle);

    // the below is the same as the above
    /*
    complex ret;
    ret.real(cos(float(delay)*angle));
    ret.imag(sin(float(delay)*angle));
    return ret;
    */
}

// Not called, but shows how to apply a filter
void ApplyFilter(const std::vector<float>& input, std::vector<float>& output, float a0, float alpha1, float alpha2, float b1, float b2)
{
    output.resize(input.size());

    const float a1 = alpha1 * a0;
    const float a2 = alpha2 * a0;

    float sample_n_2 = 0.0f;
    float sample_n_1 = 0.0f;
    float output_n_2 = 0.0f;
    float output_n_1 = 0.0f;
    for (size_t index = 0, count = input.size(); index < count; ++index)
    {
        float sample_n = input[index];
        float output_n = sample_n * a0 + sample_n_1 * a1 + sample_n_2 * a2 - output_n_1 * b1 - output_n_2 * b2;
        output[index] = output_n;
        sample_n_2 = sample_n_1;
        sample_n_1 = sample_n;
        output_n_2 = output_n_1;
        output_n_1 = output_n;
    }
}

// spit out a CSV with details of filter:
//  * difference equation
//  * location of zero
//  * frequency and phase response
void ReportFilter(const char* fileName, float a0, float alpha1, float alpha2, float b1, float b2)
{
    FILE* file = nullptr;
    fopen_s(&file, fileName, "wt");

    const float a1 = alpha1 * a0;
    const float a2 = alpha2 * a0;

    // calculate the 2 zeroes
    complex zero1, zero2;
    {
        complex left, right;
        left.real(-alpha1 / 2.0f);

        float discriminant = alpha1 * alpha1 - 4.0f * alpha2;
        if (discriminant < 0.0f)
            right.imag(sqrt(-discriminant) / 2.0f);
        else
            right.imag(sqrt(discriminant) / 2.0f);

        zero1 = left - right;
        zero2 = left + right;
    }

    // calculate the 2 poles
    complex pole1, pole2;
    {
        complex left, right;
        left.real(-b1 / 2.0f);

        float discriminant = b1 * b1 - 4.0f * b2;
        if (discriminant < 0.0f)
            right.imag(sqrt(-discriminant) / 2.0f);
        else
            right.real(sqrt(discriminant) / 2.0f);

        zero1 = left - right;
        zero2 = left + right;
    }

    fprintf(file, "\"Frequency\",\"Amplitude\",\"Phase\",");
    fprintf(file, "\"\",\"a0 = %f, alpha1 = %f, alpha2 = %f, b1 = %f, b2 = %f\"\n", a0, alpha1, alpha2, b1, b2);

    // calculate frequency & phase response
    std::vector<complex> response(c_numFrequencies);
    for (size_t index = 0; index < c_numFrequencies; ++index)
    {
        float percent = float(index) / float(c_numFrequencies - 1);
        float angle = percent * c_pi;

        complex response = a0 * (1.0f + alpha1 * Z(-1, angle) + alpha2 * Z(-2, angle));

        fprintf(file, "\"%f\",\"%f\",\"%f\"", percent, std::abs(response), atan2(response.imag(), response.real()));

        if (index == 1)
            fprintf(file, ",\"\",\"output[index] = input[index] * %f + input[index-1] * %f + input[index-1] * %f - output[index-1] * %f - output[index-2] * %f\"\n", a0, a1, a2, b1, b2);
        else if (index == 3)
            fprintf(file, ",\"\",\"Zeroes = %f + %fi, %f + %fi\"\n", zero1.real(), zero1.imag(), zero2.real(), zero2.imag());
        else if (index == 5)
            fprintf(file, ",\"\",\"Poles = %f + %fi, %f + %fi\"\n", pole1.real(), pole1.imag(), pole2.real(), pole2.imag());
        else
            fprintf(file, "\n");
    }

    fclose(file);
}

void ReportOscillator(const char* fileName, float radiansPerSample)
{
    // calculate zeroes and filter coefficients
    complex zero1 = std::polar<float>(1.0f, radiansPerSample);

    complex zero2;
    zero2.real(zero1.real());
    zero2.imag(-zero1.imag());

    // get the polynomial coefficients by doing FOIL in... (x-zero1)(x-zero2)
    // b0 is the x^2 coefficient
    // b1 is the x coefficient
    // b2 is the constant

    //float b0 = 1.0f;                                                                   // Always 1.0
    float b1 = -(zero1 + zero2).real();  // float b1 = -2.0f * cos(radiansPerSample)     // Same result either way you calculate it
    //float b2 = 1.0f;  //float b2 = (-zero1 * -zero2).real();                           // Always 1.0 too. Numerical issues will make them be close but not exact.

    // run oscillator
    std::vector<float> output(c_oscillatorSamples);
    {
        // Initial state of previous output samples lets you set the starting phase of the wave.
        // We are starting the phase at 0 radians here.
        float output_n_2 = cos(-radiansPerSample * 2.0f);
        float output_n_1 = cos(-radiansPerSample * 1.0f);
        for (size_t index = 0; index < c_oscillatorSamples; ++index)
        {
            float output_n = -output_n_1 * b1 - output_n_2;
            output[index] = output_n;
            output_n_2 = output_n_1;
            output_n_1 = output_n;
        }
    }

    // get actual values
    std::vector<float> outputActual(c_oscillatorSamples);
    for (size_t index = 0; index < c_oscillatorSamples; ++index)
        outputActual[index] = cosf(float(index) * radiansPerSample);

    // write data to the csv
    FILE* file = nullptr;
    fopen_s(&file, fileName, "wt");

    fprintf(file, "\"Sample Index\",\"Output\",\"Actual\"\n");

    int index = 0;
    for (float f : output)
    {
        fprintf(file, "\"%i\",\"%f\",\"%f\",\n", index, f, outputActual[index]);
        index++;
    }

    fclose(file);
}

int main(int argc, char**argv)
{
    // Oscillators
    {
        ReportOscillator("osc_0.csv", DegreesToRadians(0.0f));
        ReportOscillator("osc_90.csv", DegreesToRadians(90.0f));
        ReportOscillator("osc_180.csv", DegreesToRadians(180.0f));
        ReportOscillator("osc_45.csv", DegreesToRadians(45.0f));
        ReportOscillator("osc_n45.csv", DegreesToRadians(-45.0f));
        ReportOscillator("osc_315.csv", DegreesToRadians(315.0f));
        ReportOscillator("osc_10.csv", DegreesToRadians(10.0f));
    }

    // TODO: these below are old

    /*
    // order 2 filters
    {
        // a low pass filter
        ReportFilter("2_lpf.csv", 0.5f, 2.0f, 1.22f);

        // a high pass filter
        ReportFilter("2_hpf.csv", 0.5f, -1.6f, 0.8f);

        // a notch filter at 1/2 nyquist
        ReportFilter("2_notch.csv", 0.5f, 0.0f, 1.0f);
    }
    */

    return 0;
}

/*

TODO:

* adapt to IIR (biquad)
* find some good params that are worth generating csv's for

* are the poles and zeroes correct? Verify against demo.


BLOG NOTES:
* initial state lets you set phase of the wave
* make a super simple "here's how to make a cosine wave" oscillator thing. mention numerical drift.

*/