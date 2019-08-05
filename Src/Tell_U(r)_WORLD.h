
#pragma once

#include "stdafx.h"

#include "../Src/ExternalReference.h"


typedef struct
{
	std::string flagName;
	double flagValue;
	int defaultValue;
	int maxValue;
	int minValue;

}_flagValue;

namespace Tell_Ur_WORLD
{
	int makeVoiceParameterByWorld(std::string& inputWaveFile, int HarvestFlag);
	void InitializeFlags(std::vector<_flagValue>& flagValues);
	void setFlag(std::string& flags, std::vector<_flagValue>& flagValues);
	void process_breathiness(int consonantEnd, std::vector<std::vector<double>>& aperiodicity, int BflagNum, int YflagNum);
	void process_gFlag(int& fs, int& fft_size, std::vector<std::vector<double>>& spectrogram, double& gFlagNum);
	void process_PitchBend(int& consonantEnd, double pitchBendRate, double& frame_period, ARGUMENTS& Arg, std::vector<double>& f0);
	void relayToSynthesizer(int consonantEnd, double& pitchBendRate, double& velocity, WorldAnalyzedParameters& WP, std::vector<double>& wdata);
	void process_Volume(int& consonantEnd, std::vector<double>& wdata, double& volume, double& PflagNum);
}