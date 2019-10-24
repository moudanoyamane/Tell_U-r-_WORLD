
#pragma once

#include "stdafx.h"

#include "../Src/ExternalReference.h"


namespace Tell_Ur_WORLD
{
	struct _flagValue
	{
		std::string flag_name;
		double flag_value;
		int default_value;
		int max_value;
		int min_value;
	};

	int MakeVoiceParameterByWorld(std::string& input_wavefile, int harvest_flag);
	void InitializeFlags(std::vector<Tell_Ur_WORLD::_flagValue>& flag_values);
	void SetFlag(std::string& flags, std::vector<Tell_Ur_WORLD::_flagValue>& flag_values);
	void AperiodicityModifiction(int consonant_end, std::vector<std::vector<double>>& aperiodicity, int B_flag_num, int Y_flag_num);
	void SpectrumModification(int& fs, int& fft_size, std::vector<std::vector<double>>& spectrogram, double& g_flag_num);
	void F0Modification(int& consonant_end, double pitch_bend_rate, double& frame_period, eight_c_and_c::Tool2Arguments& tool2_argments, std::vector<double>& f0);
	void RelayToSynthesizer(int consonant_end, const double& pitch_bend_rate, const double& velocity, jewelwasp::VoiceParameters& voice_parameters, std::vector<double>& wdata);
	void process_Volume(const int& consonant_end, std::vector<double>& wdata, const double& volume, const double& P_flag_num);
}