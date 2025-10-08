


inline    std::size_t             radiolocation::WGaussianNoise::samples() const
{
	return this->m_samples;
}



inline    double                  radiolocation::WGaussianNoise::mean() const
{
	return this->m_mean;
}

inline   double                  radiolocation::WGaussianNoise::variance() const
{
	return this->m_variance;
}

inline    std::function<double(double)>      radiolocation::WGaussianNoise::WaveformGenerator() const
{
	return this->m_oWaveformGenerator;
}
