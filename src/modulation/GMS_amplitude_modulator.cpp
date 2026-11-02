
#include <cmath>
#include "GMS_amplitude_modulator.h"


gms::radiolocation
::amplitude_modulator_f32_t
::amplitude_modulator_f32_t(const float carrierFreqHz, 
                            const float samplingRateHz)
    : m_carrier_freq(carrierFreqHz),
      m_sampling_rate(samplingRateHz),
      m_time(0.0f),
      m_dist(-1.0f, 1.0f)
{
    // Seed the Mersenne Twister engine using a hardware random device
    std::random_device rd;
    m_rng.seed(rd());
}

gms::radiolocation
::amplitude_modulator_f32_t 
::amplitude_modulator_f32_t(const amplitude_modulator_f32_t &other) 
:
m_carrier_freq(other.m_carrier_freq),
m_sampling_rate(other.m_sampling_rate),
m_time(other.m_time),
m_dist(other.m_dist),
m_rng(other.m_rng)
{

}

gms::radiolocation
::amplitude_modulator_f32_t 
::amplitude_modulator_f32_t(amplitude_modulator_f32_t &&other) noexcept(true) 
:
m_carrier_freq(other.m_carrier_freq),
m_sampling_rate(other.m_sampling_rate),
m_time(other.m_time),
m_dist(std::move(other.m_dist)),
m_rng(std::move(other.m_rng))
{

}

gms::radiolocation
::amplitude_modulator_f32_t
::~amplitude_modulator_f32_t()
{

}

gms::radiolocation
::amplitude_modulator_f32_t &
gms::radiolocation
::amplitude_modulator_f32_t
::operator=(const amplitude_modulator_f32_t &other) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_carrier_freq = other.m_carrier_freq;
    this->m_sampling_rate= other.m_sampling_rate;
    this->m_time         = other.m_time;
    this->m_dist         = other.m_dist;
    this->m_rng          = other.m_rng; // unsafe statidtically
    return (*this);
}

gms::radiolocation
::amplitude_modulator_f32_t &
gms::radiolocation
::amplitude_modulator_f32_t
::operator=(amplitude_modulator_f32_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_carrier_freq = other.m_carrier_freq;
    this->m_sampling_rate= other.m_sampling_rate;
    this->m_time         = other.m_time;
    this->m_dist         = std::move(other.m_dist);
    this->m_rng          = std::move(other.m_rng); // unsafe statidtically
    return (*this);
}



float 
gms::radiolocation
::amplitude_modulator_f32_t
::generate_next_sample(float message_value) 
{
    if (message_value > 1.0f)  message_value = 1.0f;
    if (message_value < -1.0f) message_value = -1.0f;
    float carrier = std::cos(2.0f*3.1415926535897932384626433832795f*this->m_carrier_freq*this->m_time);
    float am_sample = (1.0f + message_value) * carrier;
    this->m_time += (1.0f / this->m_sampling_rate);
    return am_sample;
}

float 
gms::radiolocation
::amplitude_modulator_f32_t
::generate_random_message_value() 
{
    return this->m_dist(m_rng);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

gms::radiolocation
::amplitude_modulator_f64_t
::amplitude_modulator_f64_t(const double carrierFreqHz, 
                            const double samplingRateHz)
    : m_carrier_freq(carrierFreqHz),
      m_sampling_rate(samplingRateHz),
      m_time(0.0),
      m_dist(-1.0, 1.0)
{
    // Seed the Mersenne Twister engine using a hardware random device
    std::random_device rd;
    m_rng.seed(rd());
}

gms::radiolocation
::amplitude_modulator_f64_t 
::amplitude_modulator_f64_t(const amplitude_modulator_f64_t &other) 
:
m_carrier_freq(other.m_carrier_freq),
m_sampling_rate(other.m_sampling_rate),
m_time(other.m_time),
m_dist(other.m_dist),
m_rng(other.m_rng)
{

}

gms::radiolocation
::amplitude_modulator_f64_t 
::amplitude_modulator_f64_t(amplitude_modulator_f64_t &&other) noexcept(true) 
:
m_carrier_freq(other.m_carrier_freq),
m_sampling_rate(other.m_sampling_rate),
m_time(other.m_time),
m_dist(std::move(other.m_dist)),
m_rng(std::move(other.m_rng))
{

}

gms::radiolocation
::amplitude_modulator_f64_t
::~amplitude_modulator_f64_t()
{

}

gms::radiolocation
::amplitude_modulator_f64_t &
gms::radiolocation
::amplitude_modulator_f64_t
::operator=(const amplitude_modulator_f64_t &other) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_carrier_freq = other.m_carrier_freq;
    this->m_sampling_rate= other.m_sampling_rate;
    this->m_time         = other.m_time;
    this->m_dist         = other.m_dist;
    this->m_rng          = other.m_rng; // unsafe statidtically
    return (*this);
}

gms::radiolocation
::amplitude_modulator_f64_t &
gms::radiolocation
::amplitude_modulator_f64_t
::operator=(amplitude_modulator_f64_t &&other) noexcept(true) 
{
    if(__builtin_expect(this==&other,0)) { return (*this);}
    this->m_carrier_freq = other.m_carrier_freq;
    this->m_sampling_rate= other.m_sampling_rate;
    this->m_time         = other.m_time;
    this->m_dist         = std::move(other.m_dist);
    this->m_rng          = std::move(other.m_rng); // unsafe statidtically
    return (*this);
}



float 
gms::radiolocation
::amplitude_modulator_f64_t
::generate_next_sample(double message_value) 
{
    if (message_value > 1.0)  message_value = 1.0;
    if (message_value < -1.0) message_value = -1.f;
    double carrier = std::cos(2.0*3.1415926535897932384626433832795*this->m_carrier_freq*this->m_time);
    double am_sample = (1.0 + message_value) * carrier;
    this->m_time += (1.0 / this->m_sampling_rate);
    return am_sample;
}

float 
gms::radiolocation
::amplitude_modulator_f64_t
::generate_random_message_value() 
{
    return this->m_dist(m_rng);
}