
//Tell UTAU(r) WORLD 

/*
	ビルドにはWaveFileIO(https://github.com/moudanoyamane/WaveFileIO),8C&C(https://github.com/moudanoyamane/8C-C),
	JewelWasp(https://github.com/moudanoyamane/JewelWasp),WORLD(https://github.com/mmorise/World)が必要です.
*/
#pragma once

#include "../Src/stdafx.h"
#include "../Src/Tell_U(r)_WORLD.h"


int main(int argc, char* argv[])
{
	std::cout << "Tell U(r) WORLD Ver.1.0.0" << std::endl;
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

	ARGUMENTS Arg;
	double pitchBendRate;
	Eight_C_C::SetArguments(argv, Arg, pitchBendRate);
#ifdef _DEBUG
	std::cout << "complate SetArgments" << std::endl;
#endif

	std::vector<_flagValue> flagValues;
	Tell_Ur_WORLD::setFlag(Arg.flags, flagValues);

	std::string vpbwFN = Arg.inputWaveFile;
	vpbwFN.erase(vpbwFN.size() - 4, vpbwFN.size());
	vpbwFN += ".vpbw";

	WaveHeader WH;
	std::vector<double> wdata;

	WorldAnalyzedParameters WAP;

	FILE* pFile;

	int errorFlag = 1;

	if (flagValues[flagValues.size() - 1].flagValue == 0)
	{
		errorFlag = JewelWasp::readVoiceParameter(Arg.inputWaveFile, Arg.outputWaveFile, WAP, Arg.offset, Arg.blank);
	}

	if (errorFlag != 0)
	{
#ifdef _DEBUG
		std::cout << "make voice parameter" << std::endl;
#endif
		if (Tell_Ur_WORLD::makeVoiceParameterByWorld(Arg.inputWaveFile, flagValues[flagValues.size() - 1].flagValue) != 0)
		{
			return -1;
		}

		JewelWasp::readVoiceParameter(Arg.inputWaveFile, Arg.outputWaveFile, WAP, Arg.offset, Arg.blank);
	}

#ifdef _DEBUG
	std::cout << "complate readVoiceParameter" << std::endl;
#endif
	// ParameterModification
	int consonantEnd = Arg.consonant / WAP.frame_period;
	Arg.velocity = pow(2, (Arg.velocity - 100) / 100);

	Tell_Ur_WORLD::process_breathiness(consonantEnd, WAP.aperiodicity[0], flagValues[1].flagValue, flagValues[2].flagValue);

	if (flagValues[0].flagValue != flagValues[0].defaultValue)
	{
		flagValues[0].flagValue = pow(10, (-flagValues[0].flagValue) / 200);
		Tell_Ur_WORLD::process_gFlag(WAP.fs, WAP.fft_size, WAP.spectrogram[0], flagValues[0].flagValue);
	}

	consonantEnd /= Arg.velocity;
	Tell_Ur_WORLD::process_PitchBend(consonantEnd, pitchBendRate, WAP.frame_period, Arg, WAP.f0[0]);
#ifdef _DEBUG
	std::cout << "complate ParameterModification" << std::endl;
#endif

	//Synth
	wdata.resize(Arg.length * WAP.fs / 1000);
	Tell_Ur_WORLD::relayToSynthesizer(consonantEnd, pitchBendRate, Arg.velocity, WAP, wdata);
#ifdef _DEBUG
	std::cout << "complate Synthesis" << std::endl;
#endif

	//Vol
	Tell_Ur_WORLD::process_Volume(consonantEnd, wdata, Arg.volume, flagValues[3].flagValue);

	if (Arg.blank >= 0)
	{
		wdata.resize(wdata.size() - WAP.zeroPadding);
	}

	WH.Channels = 1;
	WH.SamplingRate = WAP.fs;
	WH.BitsPerSample = 16;

	for (int i = 0; i < wdata.size(); i++)
	{
		wdata[i] = wdata[i] * (1 << (WH.BitsPerSample - 1));
	}

	WaveFileIO::saveWave(Arg.outputWaveFile, WH, wdata, 0, 0);

	return 0;
}
