#include <SDL2/SDL.h>
#include <iostream>
#include <cstdint>
#include <cmath>
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
};
