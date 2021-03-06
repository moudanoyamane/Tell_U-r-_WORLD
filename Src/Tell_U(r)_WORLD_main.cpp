﻿
//Tell UTAU(r) WORLD 

/*
	ビルドにはWaveFileIO(https://github.com/moudanoyamane/Wavefile_IO),8C&C(https://github.com/moudanoyamane/Eight_C_and_C),
	JewelWasp(https://github.com/moudanoyamane/Jewelwasp),WORLD(https://github.com/mmorise/World)が必要です.
*/
#pragma once

#include "../Src/stdafx.h"
#include "../Src/Tell_U(r)_WORLD.h"


int main(int argc, char* argv[])
{
	std::cout << " Synthesis by Tell U(r) WORLD" << std::endl;
	if (argc == 0)
	{
		std::cout << "argments errer." << std::endl;

		return 0;
	}
#ifdef _DEBUG
	for (size_t i = 1; i < argc; i++)
	{
		std::cout << argv[i] << std::endl;
	}
#endif

	eight_c_and_c::Tool2Arguments tool2_argments;
	double pitch_bend_rate;
	eight_c_and_c::SetArguments(argv, tool2_argments, pitch_bend_rate);

#ifdef _DEBUG
	std::cout << "Completed seting argments" << std::endl;
#endif

	std::vector<Tell_Ur_WORLD::_flagValue> flag_values;
	Tell_Ur_WORLD::SetFlag(tool2_argments.flags, flag_values);

	std::string  vpbw_filename = tool2_argments.input_wavefile;
	vpbw_filename = vpbw_filename.erase(vpbw_filename.size() - 4, vpbw_filename.size()) + ".vpbw";

	std::vector<double> wdata;
	
	jewelwasp::VoiceParameters voice_parameters;

	int errorFlag = 1;

	if (flag_values[flag_values.size() - 1].flag_value == 0)
	{
		errorFlag = jewelwasp::ReadVoiceParameter(vpbw_filename, voice_parameters,
			tool2_argments.offset_length, -tool2_argments.end_blank_length);
	}

	if (errorFlag != 0)
	{
#ifdef _DEBUG
		std::cout << "Makeing voice parameters" << std::endl;
#endif
		if (Tell_Ur_WORLD::MakeVoiceParameterByWorld(tool2_argments.input_wavefile, flag_values[flag_values.size() - 1].flag_value) != 0)
		{
			return -1;
		}
#ifdef _DEBUG
		std::cout << "Completed Makeing voice parameters" << std::endl;
#endif
		errorFlag = jewelwasp::ReadVoiceParameter(vpbw_filename, voice_parameters,
			tool2_argments.offset_length, -tool2_argments.end_blank_length);
	}

#ifdef _DEBUG
	std::cout << "Completed loading voice parameters" << std::endl;
#endif
	// ParameterModification
	int consonant_end = tool2_argments.consonant_length / voice_parameters.frame_period;
	tool2_argments.velocity = pow(2, (tool2_argments.velocity - 100) / 100);

	Tell_Ur_WORLD::AperiodicityModifiction(consonant_end, voice_parameters.aperiodicity[0], flag_values[1].flag_value, flag_values[2].flag_value);

	if (flag_values[0].flag_value != flag_values[0].default_value)
	{
		flag_values[0].flag_value = pow(10, (-flag_values[0].flag_value) / 200);
		Tell_Ur_WORLD::SpectrumModification(voice_parameters.fs, voice_parameters.fft_size, voice_parameters.spectrogram[0], flag_values[0].flag_value);
	}

	consonant_end /= tool2_argments.velocity;
	Tell_Ur_WORLD::F0Modification(consonant_end, pitch_bend_rate, voice_parameters.frame_period, tool2_argments, voice_parameters.f0[0]);
#ifdef _DEBUG
	std::cout << "Completed parameters modification" << std::endl;
#endif

	//Synth
	wdata.resize(tool2_argments.output_length * voice_parameters.fs / 1000);
	Tell_Ur_WORLD::RelayToSynthesizer(consonant_end, pitch_bend_rate, tool2_argments.velocity, voice_parameters, wdata);
#ifdef _DEBUG
	std::cout << "Completed synthesis" << std::endl;
#endif

	//Vol
	Tell_Ur_WORLD::process_Volume(consonant_end, wdata, tool2_argments.volume, flag_values[3].flag_value);

	if (tool2_argments.end_blank_length >= 0)
	{
		wdata.resize(wdata.size() - voice_parameters.zero_padding);
	}
	
	for (int i = 0; i < wdata.size(); i++)
	{
		wdata[i] = wdata[i] * (1 << (15));
	}

	wavefile_io::SaveWave(tool2_argments.output_wavefile, wdata, 44100, 16);
#ifdef _DEBUG
	std::cout << "Completed writeing .wav file" << std::endl;
#endif

	return 0;
}
