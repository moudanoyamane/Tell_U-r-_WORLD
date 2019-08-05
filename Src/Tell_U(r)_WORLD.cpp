
#pragma once

#include "../Src/stdafx.h"
#include "Tell_U(r)_WORLD.h"

#define ascii_A 65
#define A1_Frq 55

namespace {
	double lerp(double a, double b, double t)
	{
		//c++20
		return a + t * (b - a);
	}
}

namespace Tell_Ur_WORLD
{

	int makeVoiceParameterByWorld(std::string& inputWaveFile, int HarvestFlag)
	{
		WaveHeader WH;
		std::vector<double> wdata;

		if (WaveFileIO::loadWave(inputWaveFile, wdata, WH, 0, 0) != 0)
		{
			return 1;
		}

		for (int i = 0; i < wdata.size(); i++)
		{
			wdata[i] = wdata[i] / (1 << (WH.BitsPerSample - 1));
		}

		WorldAnalyzedParameters WAP;
		std::vector<double*> AP_p;
		std::vector<double*> SP_p;

		CheapTrickOption CTOP{ 0 };
		D4COption D4COP{ 0 };

		WAP.fs = WH.SamplingRate;

		std::vector<double> DioF0;
		DioOption DioOP = { 0 };
		HarvestOption HrvOP{ 0 };

		switch (HarvestFlag)
		{
		case 0:
			// Dio & StoneMask
			InitializeDioOption(&DioOP);
			DioOP.speed = 1;
			DioOP.f0_floor = 71.0;
			DioOP.allowed_range = 0.1;
			DioOP.frame_period = WAP.frame_period;

			WAP.f0[0].resize(GetSamplesForDIO(WAP.fs, wdata.size(), WAP.frame_period));

			JewelWasp::zeroPadding(WAP, wdata);

			DioF0.resize(WAP.f0[0].size());

			Dio(&wdata[0], wdata.size(), WAP.fs, &DioOP, &WAP.time_axis[0], &DioF0[0]);
			StoneMask(&wdata[0], wdata.size(), WAP.fs, &WAP.time_axis[0], &DioF0[0], WAP.f0[0].size(), &WAP.f0[0][0]);
			break;

		case 1:
			//Harvest
			InitializeHarvestOption(&HrvOP);
			HrvOP.frame_period = WAP.frame_period;
			HrvOP.f0_floor = 71.0;

			WAP.f0[0].resize(GetSamplesForHarvest(WAP.fs, wdata.size(), WAP.frame_period));

			JewelWasp::zeroPadding(WAP, wdata);

			Harvest(&wdata[0], wdata.size(), WAP.fs, &HrvOP, &WAP.time_axis[0], &WAP.f0[0][0]);

			break;

		default: break;
		}


		InitializeCheapTrickOption(WAP.fs, &CTOP);

		CTOP.f0_floor = 71.0;
		WAP.fft_size = CTOP.fft_size;

		WAP.spectrogram[0].resize(WAP.f0[0].size());
		SP_p.resize(WAP.f0[0].size());

		for (int i = 0; i < WAP.f0[0].size(); ++i)
		{
			WAP.spectrogram[0][i].resize(WAP.fft_size / 2 + 1);
			SP_p[i] = &WAP.spectrogram[0][i][0];
		}

		CheapTrick(&wdata[0], wdata.size(), WAP.fs, &WAP.time_axis[0], &WAP.f0[0][0], WAP.f0[0].size(), &CTOP, &SP_p[0]);

		InitializeD4COption(&D4COP);

		WAP.aperiodicity[0].resize(WAP.f0[0].size());
		AP_p.resize(WAP.f0[0].size());

		for (int i = 0; i < WAP.f0[0].size(); ++i)
		{
			WAP.aperiodicity[0][i].resize(WAP.fft_size / 2 + 1);
			AP_p[i] = &WAP.aperiodicity[0][i][0];
		}

		D4C(&wdata[0], wdata.size(), WAP.fs, &WAP.time_axis[0], &WAP.f0[0][0], WAP.f0[0].size(), WAP.fft_size, &D4COP, &AP_p[0]);

		return JewelWasp::writeVoiceParameter(inputWaveFile, WAP);
	}


	void InitializeFlags(std::vector<_flagValue>& flagValues)
	{
		std::vector<std::string> FLAG{ "g","B","Y","P","G" };
		std::vector<__int16> DEFAULT{ 0,50,100,86,0 };
		std::vector<__int16> MAX{ 100,100,200,100,1 };
		std::vector<__int16> MIN{ -100,0,0,0,0 };

		flagValues.resize(FLAG.size());

		for (int i = 0; i < FLAG.size(); i++)
		{
			flagValues[i].flagName = FLAG[i];
			flagValues[i].flagValue = DEFAULT[i];
			flagValues[i].defaultValue = DEFAULT[i];
			flagValues[i].maxValue = MAX[i];
			flagValues[i].minValue = MIN[i];
		}
	}


	void setFlag(std::string& flags, std::vector<_flagValue>& flagValues)
	{
		InitializeFlags(flagValues);

		std::smatch match;
		std::regex reg;
		std::string str;
		int i;

		for (i = 0; i < flagValues.size() - 1; i++)
		{
			reg = (flagValues[i].flagName + "((\\d+)|(\\+\\d+)|(-\\d+))");
			if (regex_search(flags, match, reg))
			{
				str = match[0];
				reg = ("-(\\d+)");
				if (regex_search(str, match, reg)) //音符固有のフラグ(優先)が先にに来るため後方マッチは無視
				{
					flagValues[i].flagValue = Limitter(stof(match[0]), flagValues[i].minValue, flagValues[i].maxValue);
				}
				else
				{
					reg = ("(\\d+)");
					regex_search(str, match, reg);
					flagValues[i].flagValue = Limitter(stof(match[0]), flagValues[i].minValue, flagValues[i].maxValue);
				}
			}
		}

		reg = (flagValues[i].flagName);// only large"G"flag
		if (regex_search(flags, match, reg))
		{
			flagValues[i].flagValue = flagValues[i].maxValue;
		}
		else
		{
			flagValues[i].flagValue = flagValues[i].minValue;
		}
	}


	void process_breathiness(int consonantEnd, std::vector<std::vector<double>>& aperiodicity, int BflagNum, int YflagNum)
	{
		int i;
		int l;

		BflagNum /= 50;
		YflagNum /= 100;

		consonantEnd = std::min((int)aperiodicity.size(), consonantEnd);

		for (i = 0; i < consonantEnd; i++)
		{
			for (l = 0; l < aperiodicity[0].size(); l++)
			{
				aperiodicity[i][l] *= (BflagNum);
			}
		}

		for (; i < aperiodicity.size(); i++)
		{
			for (l = 0; l < aperiodicity[0].size(); l++)
			{
				aperiodicity[i][l] *= (BflagNum * YflagNum);
			}
		}
	}


	void process_gFlag(int& fs, int& fft_size, std::vector<std::vector<double>>& spectrogram, double& gFlagNum)
	{
		std::vector<double> freq_axis1(fft_size);
		std::vector<double> freq_axis2(fft_size);
		std::vector<double> spectrum1(fft_size);
		std::vector<double> spectrum2(fft_size);

		for (int i = 0; i <= fft_size / 2; ++i) {
			freq_axis1[i] = gFlagNum * i / fft_size * fs;
			freq_axis2[i] = static_cast<double>(i) / fft_size * fs;
		}

		for (int i = 0; i < spectrogram.size(); ++i)
		{
			for (int j = 0; j <= fft_size / 2; ++j)
			{
				spectrum1[j] = log(spectrogram[i][j]);
			}

			// @world::matlabfunctions.cpp
			interp1(&freq_axis1[0], &spectrum1[0], spectrogram[0].size(), &freq_axis2[0], spectrogram[0].size(), &spectrum2[0]);

			for (int j = 0; j <= fft_size / 2; ++j)
			{
				spectrogram[i][j] = exp(spectrum2[j]);
			}

			if (gFlagNum < 1.0)
			{
				for (int j = static_cast<int>(fft_size / 2.0 * gFlagNum); j <= fft_size / 2; ++j)
				{
					spectrogram[i][j] = spectrogram[i][static_cast<int>(fft_size / 2.0 * gFlagNum) - 1];
				}
			}
		}
	}

	/*
	void process_TimeStretch()
	{

	}
	*/
	void process_PitchBend(int& consonantEnd, double pitchBendRate, double& frame_period, ARGUMENTS& Arg, std::vector<double>& f0)
	{
		std::vector<double> tempF0 = f0;
		std::string code;
		int octave;
		int noteDiff;
		int sharpFlag = 0;
		double baseFrq;
		double referenceFrq = 0;
		int offset;
		double offsetLength;
		double expansionRatio;
		double lerpRate;
		int i;
		int j;

		pitchBendRate *= 1000;

		std::smatch match;

		std::regex reg("([A-Z])");
		if (regex_search(Arg.noteNum, match, reg))
		{
			code = match[0];
		}
		reg = ("(#)");
		if (regex_search(Arg.noteNum, match, reg))
		{
			sharpFlag++;
		}
		reg = ("(\\d+)");
		if (regex_search(Arg.noteNum, match, reg))
		{
			octave = stof(match[0]);
		}

		noteDiff = (code[0] - ascii_A);

		if (noteDiff >= 2)
		{
			noteDiff -= 7;
		}
		noteDiff *= 2;
		if (noteDiff <= -6)
		{
			noteDiff++;
		}
		noteDiff += sharpFlag;

		baseFrq = (static_cast<int>(A1_Frq) << (octave - 1)) * pow(2, ((double)noteDiff / 12));
		//UTAUではB0以下は指定不可のため

		f0.resize(Arg.pitchBend.size());

		if (Arg.modulation != 0)
		{
			double weight;
			double sumWeight = 0;

			if (f0[0] != 0)
			{
				referenceFrq = f0[0];
				sumWeight++;
			}

			for (i = 1; i < f0.size(); i++)
			{
				if (f0[i] != 0)
				{
					weight = 2 - std::min((double)(abs(f0[i - 1] - f0[i]) / f0[i]), 1.0);

					referenceFrq += f0[i] * weight;
					sumWeight += weight;
				}
			}
			if (sumWeight > 0)
			{
				referenceFrq /= sumWeight;
			}

		}

		expansionRatio = (double)frame_period / (pitchBendRate * Arg.velocity);

		consonantEnd = std::min((int)f0.size(), consonantEnd);

		j = 0;
		for (i = 0; i < consonantEnd; ++i)
		{
			for (; i > (j + 1) * expansionRatio; j++);

			if (tempF0[j] * tempF0[j + 1] != 0)
			{
				lerpRate = (i - j * expansionRatio) / expansionRatio; //(((j + 1) *  expansionRatio) - j *  expansionRatio) ==  expansionRatio
				f0[i] = baseFrq * pow(2, ((double)Arg.pitchBend[i] / 1200))
					+ (lerp(tempF0[j], tempF0[j + 1], lerpRate) - referenceFrq) * Arg.modulation / 100;
			}
			else if (tempF0[j] + tempF0[j + 1] != 0)
			{
				f0[i] = baseFrq * pow(2, ((double)Arg.pitchBend[i] / 1200))
					+ (tempF0[j] + tempF0[j + 1] - referenceFrq) * Arg.modulation / 100;
			}
		}

		offset = j;
		offsetLength = j * expansionRatio;
		expansionRatio = (double)(f0.size() - 1 - consonantEnd) / (tempF0.size() - 2 - offset);

		j = 0;
		for (; i < f0.size(); i++)
		{
			for (; (i - consonantEnd) > (j + 1) * expansionRatio; j++);

			if (tempF0[j + offset] * tempF0[j + offset + 1] != 0)
			{
				lerpRate = (i - offsetLength - j * expansionRatio) / expansionRatio;
				f0[i] = baseFrq * pow(2, ((double)Arg.pitchBend[i] / 1200))
					+ (lerp(tempF0[j + offset], tempF0[j + offset + 1], lerpRate) - referenceFrq) * Arg.modulation / 100;
			}
			else if (tempF0[j + offset] + tempF0[j + offset + 1] != 0)
			{
				f0[i] = baseFrq * pow(2, ((double)Arg.pitchBend[i] / 1200))
					+ (tempF0[j] + tempF0[j + 1] - referenceFrq) * Arg.modulation / 100;
			}
		}
	}

	void relayToSynthesizer(int consonantEnd, double& pitchBendRate, double& velocity, WorldAnalyzedParameters& WAP, std::vector<double>& wdata)
	{

		WorldSynthesizer synthesizer = { 0 };
		std::vector<double*> SP_p;
		std::vector<double*> AP_p;
		int i;
		int j;
		int offset;
		double expansionRatio;

		int f0_length = WAP.f0[0].size();
		WAP.f0[0].resize(WAP.f0[0].size() + 50);//適当

		SP_p.resize(WAP.f0[0].size());
		AP_p.resize(WAP.f0[0].size());

		expansionRatio = (double)WAP.frame_period / (pitchBendRate * 1000 * velocity);

		consonantEnd = std::min(f0_length, consonantEnd);

		j = 0;
		for (i = 0; i < consonantEnd; ++i)
		{
			for (; i > j * expansionRatio; j++);

			SP_p[i] = &WAP.spectrogram[0][j][0];
			AP_p[i] = &WAP.aperiodicity[0][j][0];
		}

		offset = j;
		expansionRatio = (double)(f0_length - consonantEnd - 1) / (WAP.spectrogram[0].size() - 2 - offset);

		j = 0;
		for (; i < f0_length; ++i)
		{
			for (; i - consonantEnd > j * expansionRatio; j++);

			SP_p[i] = &WAP.spectrogram[0][j + offset][0];
			AP_p[i] = &WAP.aperiodicity[0][j + offset][0];
		}

		for (; i < WAP.f0[0].size(); ++i)
		{
			SP_p[i] = SP_p[f0_length - 1];
			AP_p[i] = AP_p[f0_length - 1];
			WAP.f0[0][i] = WAP.f0[0][f0_length - 1];
		}

		Synthesis(&WAP.f0[0][0], WAP.f0[0].size(), &SP_p[0], &AP_p[0],
			WAP.fft_size, pitchBendRate * 1000, WAP.fs, wdata.size(), &wdata[0]);
	}

	void process_Volume(int& consonantEnd, std::vector<double>& wdata, double& volume, double& PflagNum)
	{
		std::pair<double*, double*> wdata_minmax = std::minmax_element(&wdata[consonantEnd], &wdata[wdata.size() - 1]);
		double tmpd[3];
		tmpd[0] = std::abs(*wdata_minmax.first);
		tmpd[1] = std::abs(*wdata_minmax.second);
		tmpd[2] = std::max(tmpd[0], tmpd[1]);

		double amplification_factor = volume / 100.0 * lerp(1.0, (0.5 / tmpd[2]), PflagNum / 100.0);

		for (int i = 0; i < wdata.size(); i++)
		{
			wdata[i] *= amplification_factor;
		}
	}

}