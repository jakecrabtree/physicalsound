#include <SDL2/SDL.h>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <math.h>
#include <cassert>
//Code based off answer 3 from: 

class AudioPlayer {
	public:
	int sampleID = 0;
	double vol = 0;
	static void getSound(void* userdata, unsigned char* raw_buffer, int bytes) {
		short* rb = (short*) raw_buffer;
		AudioPlayer& ap = ((AudioPlayer*)userdata)[0];
		int& currSample = ((int*)userdata)[0];
		double v = ap.vol;
		if(v > 2) {
			v = 2;
		}
		int shorts = bytes >> 1;
		for(int i = 0; i < shorts; i++) {
			double time = currSample / 44100.0;
			rb[i] = (short)(28000 * sin(v * time * 2 * 3.1415 * 441));
			currSample++;
		}
	}
	
	AudioPlayer() {
		std::cout << "CONSTRUCTED\n";
		int ival = SDL_Init(SDL_INIT_AUDIO);
		if(ival != 0) {
			std::cout << "Could not initalize audio " << ival << "\n";
		}
		SDL_AudioSpec target;
		target.freq = 44100;
		target.format = AUDIO_S16SYS;
		target.channels = 1;
		target.callback = getSound;
		target.samples = 2048;
		target.userdata = this;
		SDL_AudioSpec targot;
		ival = SDL_OpenAudio(&target, &targot);
		
		SDL_PauseAudio(0);
	}


	~AudioPlayer() {
		std::cout << "CLOSING\n";
		SDL_PauseAudio(1);
		SDL_CloseAudio();
	}

	static std::vector<double> lowPassFilter(float cutoffFrequency, float sampleFrequency, int order){
		std::vector<double> filter(order); 
		cutoffFrequency /= sampleFrequency;
		double cutoffOmega = 2 * M_PI * cutoffFrequency;
		int middle = order / 2;
		for (int i = -1*order / 2; i <= order / 2; ++i){
			if (i == 0){
				filter[i] = 2*cutoffFrequency;
			} else{
				filter[i + middle] = std::sin(cutoffOmega*i) / M_PI*i;
			}
		} 
		return filter;
	}

	static void convolutionFilter(std::vector<double>& samples, const std::vector<double>& filter, std::vector<double>& filteredSamples){
		filteredSamples.resize(samples.size());
		for (int m = 0; m < samples.size(); ++m){
			filteredSamples[m] = 0;
			for (int n = 0; n < filter.size(); ++n){
				if (n - m > 0){
					filteredSamples[m] += filter[m] * samples[n-m]; 
				}
			}
		}
	}

	static void dcBlockingFilter(const std::vector<double>& samples, std::vector<double>& filteredSamples, double lossConstant){
		filteredSamples.resize(samples.size());
		filteredSamples[0] = samples[0];
		for (int i = 1; i < samples.size(); ++i){
			filteredSamples[i] = (1.0 - lossConstant) * filteredSamples[i-1] + samples[i] - samples[i-1];
		}
	} 
};
