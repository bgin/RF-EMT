

#include <random>
#include "GMS_white_gauss_noise.h"



/*
@Brief: Creates vector of White Gaussian Noise initialized by Box-Mueller method.
*/
__thiscall radiolocation::WGaussianNoise::WGaussianNoise(_In_ struct WGNoiseParams const& params) : m_oWaveformGenerator{params.WaveformGenerator}, m_samples{ params.n_samples }, m_mean{ params.mean },
m_variance{ params.variance }
{


	double PI{ 3.14159265358979323846264 };
	this->m_oWGNoise = std::vector<std::pair<double, double>>(this->m_samples);
	//Set up radom generator.
	std::clock_t seed{ std::clock() };
	auto rand_value = std::bind(std::uniform_real_distribution<double>(0.0, 1.0), std::default_random_engine(seed));
	double rv1, rv2;
	for (std::size_t i{ 0 }; i != this->m_samples; ++i)
	{
		do
		rv1 = rand_value();
		while (rv1 == 0.0);
		rv2 = rand_value();
		double vu1{ std::sqrt(-2.0 * std::log(rv1)) * this->m_oWaveformGenerator(2.0 * PI * rv2) }; // unit normal variable
		double vr2{ this->m_mean + std::sqrt(this->m_variance) * vu1 }; // normal rand variable
		this->m_oWGNoise.operator[](i).operator=({ vu1, vr2 });
	}

}

/*
@Brief: Copy Constructor.
*/
__thiscall radiolocation::WGaussianNoise::WGaussianNoise(_In_ WGaussianNoise const& other) : m_mean{ other.m_mean }, m_samples{ other.m_samples },
m_variance{ other.m_variance }, m_oWGNoise{ other.m_oWGNoise }, m_oWaveformGenerator{ other.m_oWaveformGenerator }
{

}

/*
@Brief: Move Constructor.
*/
__thiscall radiolocation::WGaussianNoise::WGaussianNoise(_In_ WGaussianNoise &&other) : m_mean{ std::move(other.m_mean) }, m_samples{ std::move(other.m_samples) },
m_variance{ std::move(other.m_variance) }, m_oWGNoise{ std::move(other.m_oWGNoise) }, m_oWaveformGenerator{ std::move(other.m_oWaveformGenerator) }
{

}

radiolocation::WGaussianNoise &         radiolocation::WGaussianNoise::operator=(_In_ WGaussianNoise const& other)
{
	if (this == &other) return *this;

	this->m_mean = other.m_mean;
	this->m_samples = other.m_samples;
	this->m_variance = other.m_variance;
	this->m_oWaveformGenerator.operator=(other.m_oWaveformGenerator);
	this->m_oWGNoise.operator=(other.m_oWGNoise);
	return *this;
}

radiolocation::WGaussianNoise &         radiolocation::WGaussianNoise::operator=(_In_ WGaussianNoise &&other)
{
	if (this == &other) return *this;

	this->m_mean = std::move(other.m_mean);
	this->m_samples = std::move(other.m_samples);
	this->m_variance = std::move(other.m_variance);
	this->m_oWaveformGenerator.operator=(std::move(other.m_oWaveformGenerator));
	this->m_oWGNoise.operator=(other.m_oWGNoise);
	return *this;

}

std::ostream&                         radiolocation::operator<<(_In_ std::ostream& os, _In_ WGaussianNoise const& noise)
{
	os.scientific;
	std::setprecision(15);
	for (std::size_t i{ 0 }; i != noise.m_samples; ++i)
	{
		os << "Unit rand vector:" << noise.m_oWGNoise.operator[](i).first <<
			"WG Noise:" << noise.m_oWGNoise.operator[](i).second << std::endl;
	}
	return os;
}



void               radiolocation::WGaussianNoise::debug_print() const
{
	std::printf("WGaussianNoise::debug_print\n");
	std::printf("&this->m_samples=%p,this->m_samples=%u\n", &this->m_samples, this->m_samples);
	std::printf("&this->m_mean=%p,this->m_mean=%.9f\n", &this->m_mean, this->m_mean);
	std::printf("&this->m_variance=%p,this->m_variance=%.9f\n", &this->m_variance, this->m_variance);
	std::printf("&this->m_oWaveformGnerator=%p\n", &this->m_oWaveformGenerator);
	std::printf("  var=   |       WGN(var)=    \n");
	for (std::size_t i{ 0 }; i != this->m_samples; ++i)
		std::printf("%.15f, %.15f\n", this->m_oWGNoise.operator[](i).first, this->m_oWGNoise.operator[](i).second);
	std::printf("End of WGausianNoise object dump\n");
}


