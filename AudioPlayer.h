#include <SDL2/SDL.h>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <queue>
#include <vector>

class Sample {
	public:
	long time;
	double sam;
	const bool operator< (const Sample& o) const {
		return o.time < time;
	}
};

//Code based off answer 3 from: 

class AudioPlayer {
	public:
	long sampleID = 0;
	double vol = 0;
	std::priority_queue<Sample, std::vector<Sample>> pq;
	static void getSound(void* userdata, unsigned char* raw_buffer, int bytes) {
		std::cout << "SOUNDS!\n";
		short* rb = (short*) raw_buffer;
		AudioPlayer& ap = ((AudioPlayer*)userdata)[0];
		long& currSample = ((long*)userdata)[0];
		double v = ap.vol;
		if(v > 2) {
			v = 2;
		}
		int shorts = bytes >> 1;
		for(int i = 0; i < shorts; i++) {
			rb[i] = 0;
			std::cout << currSample << "\n";
			while(!ap.pq.empty() && ap.pq.top().time == currSample) {
				std::cout << ":D\n";
				rb[i] += ap.pq.top().sam;
				ap.pq.pop();
			}
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

	void addWithDelay(double sam, double delay) {
		Sample s;
		s.time = sampleID + (long)(delay * 44100);
		s.sam = sam;
		std::cout << "Adding at: " << s.time << "\n";
		pq.push(s);
	}

	~AudioPlayer() {
		std::cout << "CLOSING\n";
		SDL_PauseAudio(1);
		SDL_CloseAudio();
	}
};

