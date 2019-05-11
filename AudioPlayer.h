#include <SDL2/SDL.h>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <queue>
#include <vector>
#include <mutex>

class Sample {
	public:
	long time;
	double sam;
	const bool operator< (const Sample& o) const {
		return o.time < time;
	}
};

//Code based off answer 3 from: https://stackoverflow.com/questions/10110905/simple-sound-wave-generator-with-sdl-in-c

class AudioPlayer {
	public:
	long sampleID = 0;
	double vol = 0;
	static const int ssize = 44100 * 20;
	std::vector<double> samples;
	std::mutex lock;
	static void getSound(void* userdata, unsigned char* raw_buffer, int bytes) {
		short* rb = (short*) raw_buffer;
		AudioPlayer& ap = ((AudioPlayer*)userdata)[0];
		long& currSample = ((long*)userdata)[0];
		double v = ap.vol;
		if(v > 2) {
			v = 2;
		}
		int shorts = bytes >> 1;
		for(int i = 0; i < shorts; i++) {
			int index = (int)(currSample + i);
			if(index >= ap.samples.size()) {
				rb[i] = 0;
				currSample++;
				continue;
			}
			//std::cout << d << "\n";
			rb[i] = (short) (ap.samples[index] * 20);
			std::cout << rb[i] << "\n";
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
		if(ival != 0) {
			std::cout << "Audio could not be opened: " << ival << "\n";
			exit(0);
		}
		if(target.format != targot.format) {
			std::cout << "Format did not match\n";
			exit(0);
		}
	}

	void clear() {
		samples.clear();
	}

	void playSound() {
		sampleID = 0;
		SDL_PauseAudio(0);
	}

	void stopSound() {
		SDL_PauseAudio(1);
	}

	void addWithDelay(double sam, double delay) {
		int index = (int)(delay * 44100);
		while(index >= samples.size()) {
			samples.push_back(0);
		}
		samples[(index) % ssize] += sam;
	}

	void dumpAudio(std::ofstream& o) {
		o << samples.size() << "\n";
		for(int i = 0; i < samples.size(); i++) {
			o << (samples[i]) << "\n";
		}
	}

	~AudioPlayer() {
		std::cout << "CLOSING\n";
		SDL_PauseAudio(1);
		SDL_CloseAudio();
	}
};

