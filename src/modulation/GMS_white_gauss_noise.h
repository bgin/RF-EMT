#ifndef __GMS_WHITE_GAUSS_NOISE_H__
#define __GMS_WHITE_GAUSS_NOISE_H__ 150620161456

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "GMS_config.h"

namespace  gms
{


	class WGaussianNoise
	{


	public:

		/*
		@Brief: Constructors.
		*/
		WGaussianNoise(_In_ struct WGNoiseParams const&);

		WGaussianNoise(_In_  WGaussianNoise const&);

		WGaussianNoise(_In_  WGaussianNoise&&);

		~WGaussianNoise() = default;

		/*
		@Brief: operators
		*/

		WGaussianNoise &       operator=(_In_ WGaussianNoise const&);

		WGaussianNoise &       operator=(_In_ WGaussianNoise &&);

		friend  std::ostream&    operator<<(_In_ std::ostream&, _In_ WGaussianNoise const&);

		/*
		@Brief: Class accessor methods.
		*/

		inline     std::size_t       samples() const;

		//inline     std::vector<std::pair<double, double>> WGNoise() const;

		inline     double            mean() const;

		inline     double            variance() const;

		inline     std::function<double(double)>   WaveformGenerator() const;



		void               debug_print() const;


	

		/*
		@Brief: Number of noise vector samples.
		*/
		std::size_t  m_samples;

		/*
		@Brief: White Gaussian Noise vector.
		Contains std::pair.first = normal unit random variables vector, std::pair.second = normal random variables vector.
		*/
		std::vector<std::pair<double, double>> m_oWGNoise; // White Gaussian Noise vector

		/*
		@Brief: mean value.
		*/
		double  m_mean;

		/*
		@Brief: variance value.
		*/
		double  m_variance;

		/*
        @Brief: Noise waveform generator.
		*/
		std::function<double(double)> m_oWaveformGenerator;



	};


	struct WGNoiseParams
	{
		_In_ std::function<double(double)> WaveformGenerator;
		_In_ double mean;
		_In_ double variance;
		_In_ std::size_t n_samples;
	};

#include "GMS_white_gauss_noise.inl"
}

#endif /*_WHITE_GAUSSIAN_NOISE_H_*/
