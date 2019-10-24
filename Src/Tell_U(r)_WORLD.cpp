
#pragma once

#include "../Src/stdafx.h"
#include "Tell_U(r)_WORLD.h"

#define ascii_A 65
#define A1_Frq 55


namespace
{
	double lerp(double a, double b, double t)
	{
		//c++20
		return a + t * (b - a);
	}
}


namespace Tell_Ur_WORLD
{
	int MakeVoiceParameterByWorld(std::string& input_wavefile, int harvest_flag)
	{
		wavefile_io::WaveHeader wav_haader;
		std::vector<double> wdata;

		if (wavefile_io::LoadWave(input_wavefile, wdata, wav_haader) != 0)
		{
			return 1;
		}

		for (int i = 0; i < wdata.size(); i++)
		{
			wdata[i] = wdata[i] / (1 << (wav_haader.bits_per_sample - 1));
		}

		jewelwasp::VoiceParameters voice_parameters;
		std::vector<double*> p_aperiodicity;
		std::vector<double*> p_spectrogram;

		CheapTrickOption cheaptrick_option{ 0 };
		D4COption d4c_option{ 0 };

		voice_parameters.fs = wav_haader.sampling_rate;

		std::vector<double> dio_f0;
		DioOption dio_option = { 0 };
		HarvestOption harvest_option{ 0 };

		switch (harvest_flag)
		{
		case 0:
			// Dio & StoneMask
			InitializeDioOption(&dio_option);
			dio_option.speed = 1;
			dio_option.f0_floor = 71.0;
			dio_option.allowed_range = 0.1;
			dio_option.frame_period = voice_parameters.frame_period;

			voice_parameters.f0[0].resize(GetSamplesForDIO(voice_parameters.fs, wdata.size(), voice_parameters.frame_period));

			jewelwasp::ZeroPadding(voice_parameters, wdata);

			dio_f0.resize(voice_parameters.f0[0].size());

			Dio(&wdata[0], wdata.size(), voice_parameters.fs, &dio_option, &voice_parameters.time_axis[0], &dio_f0[0]);
			StoneMask(&wdata[0], wdata.size(), voice_parameters.fs, &voice_parameters.time_axis[0], &dio_f0[0],
				voice_parameters.f0[0].size(), &voice_parameters.f0[0][0]);
			break;

		case 1:
			//Harvest
			InitializeHarvestOption(&harvest_option);
			harvest_option.frame_period = voice_parameters.frame_period;
			harvest_option.f0_floor = 71.0;

			voice_parameters.f0[0].resize(GetSamplesForHarvest(voice_parameters.fs, wdata.size(), voice_parameters.frame_period));

			jewelwasp::ZeroPadding(voice_parameters, wdata);

			Harvest(&wdata[0], wdata.size(), voice_parameters.fs, &harvest_option, &voice_parameters.time_axis[0], &voice_parameters.f0[0][0]);

			break;

		default: break;
		}
		
		InitializeCheapTrickOption(voice_parameters.fs, &cheaptrick_option);

		cheaptrick_option.f0_floor = 71.0;
		voice_parameters.fft_size = cheaptrick_option.fft_size;

		voice_parameters.spectrogram[0].resize(voice_parameters.f0[0].size());
		p_spectrogram.resize(voice_parameters.f0[0].size());

		for (int i = 0; i < voice_parameters.f0[0].size(); ++i)
		{
			voice_parameters.spectrogram[0][i].resize(voice_parameters.fft_size / 2 + 1);
			p_spectrogram[i] = &voice_parameters.spectrogram[0][i][0];
		}

		CheapTrick(&wdata[0], wdata.size(), voice_parameters.fs, &voice_parameters.time_axis[0],
			&voice_parameters.f0[0][0], voice_parameters.f0[0].size(), &cheaptrick_option, &p_spectrogram[0]);

		InitializeD4COption(&d4c_option);

		voice_parameters.aperiodicity[0].resize(voice_parameters.f0[0].size());
		p_aperiodicity.resize(voice_parameters.f0[0].size());

		for (int i = 0; i < voice_parameters.f0[0].size(); ++i)
		{
			voice_parameters.aperiodicity[0][i].resize(voice_parameters.fft_size / 2 + 1);
			p_aperiodicity[i] = &voice_parameters.aperiodicity[0][i][0];
		}

		D4C(&wdata[0], wdata.size(), voice_parameters.fs, &voice_parameters.time_axis[0],
			&voice_parameters.f0[0][0], voice_parameters.f0[0].size(), voice_parameters.fft_size, &d4c_option, &p_aperiodicity[0]);

		std::string vpbw_filename = input_wavefile;

		vpbw_filename = vpbw_filename.erase(vpbw_filename.size() - 4, vpbw_filename.size()) + ".vpbw";

		return jewelwasp::WriteVoiceParameter(vpbw_filename, voice_parameters);
	}


	void InitializeFlags(std::vector<Tell_Ur_WORLD::_flagValue>& flag_values)
	{
		std::vector<std::string> FLAG{ "g","B","Y","P","G" };
		std::vector<__int16> DEFAULT{ 0,50,100,86,0 };
		std::vector<__int16> MAX{ 100,100,200,100,1 };
		std::vector<__int16> MIN{ -100,0,0,0,0 };

		flag_values.resize(FLAG.size());

		for (int i = 0; i < FLAG.size(); i++)
		{
			flag_values[i].flag_name = FLAG[i];
			flag_values[i].flag_value = DEFAULT[i];
			flag_values[i].default_value = DEFAULT[i];
			flag_values[i].max_value = MAX[i];
			flag_values[i].min_value = MIN[i];
		}
	}


	void SetFlag(std::string& flags, std::vector<Tell_Ur_WORLD::_flagValue>& flag_values)
	{
		InitializeFlags(flag_values);

		std::smatch match;
		std::regex reg;
		std::string str;
		int i;

		for (i = 0; i < flag_values.size() - 1; i++)
		{
			reg = (flag_values[i].flag_name + "((\\d+)|(\\+\\d+)|(-\\d+))");
			if (regex_search(flags, match, reg))
			{
				str = match[0];
				reg = ("-(\\d+)");
				if (regex_search(str, match, reg)) //音符固有のフラグ(優先)が先にに来るため後方マッチは無視
				{
					flag_values[i].flag_value = std::clamp(stoi(match[0]), flag_values[i].min_value, flag_values[i].max_value);
				}
				else
				{
					reg = ("(\\d+)");
					regex_search(str, match, reg);
					flag_values[i].flag_value = std::clamp(stoi(match[0]), flag_values[i].min_value, flag_values[i].max_value);
				}
			}
		}

		reg = (flag_values[i].flag_name);// only large"G"flag
		if (regex_search(flags, match, reg))
		{
			flag_values[i].flag_value = flag_values[i].max_value;
		}
		else
		{
			flag_values[i].flag_value = flag_values[i].min_value;
		}
	}


	void AperiodicityModifiction(int consonant_end, std::vector<std::vector<double>>& aperiodicity, int B_flag_num, int Y_flag_num)
	{
		int i;
		int l;

		B_flag_num /= 50;
		Y_flag_num /= 100;

		consonant_end = std::min((int)aperiodicity.size(), consonant_end);

		for (i = 0; i < consonant_end; i++)
		{
			for (l = 0; l < aperiodicity[0].size(); l++)
			{
				aperiodicity[i][l] *= (B_flag_num);
			}
		}

		for (; i < aperiodicity.size(); i++)
		{
			for (l = 0; l < aperiodicity[0].size(); l++)
			{
				aperiodicity[i][l] *= (B_flag_num * Y_flag_num);
			}
		}
	}


	void SpectrumModification(int& fs, int& fft_size, std::vector<std::vector<double>>& spectrogram, double& g_flag_num)
	{
		std::vector<double> freq_axis1(fft_size);
		std::vector<double> freq_axis2(fft_size);
		std::vector<double> spectrum1(fft_size);
		std::vector<double> spectrum2(fft_size);

		for (int i = 0; i <= fft_size / 2; ++i) {
			freq_axis1[i] = g_flag_num * i / fft_size * fs;
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

			if (g_flag_num < 1.0)
			{
				for (int j = static_cast<int>(fft_size / 2.0 * g_flag_num); j <= fft_size / 2; ++j)
				{
					spectrogram[i][j] = spectrogram[i][static_cast<int>(fft_size / 2.0 * g_flag_num) - 1];
				}
			}
		}
	}


	void F0Modification(int& consonant_end, double pitch_bend_rate, double& frame_period, eight_c_and_c::Tool2Arguments& tool2_argments, std::vector<double>& f0)
	{
		std::vector<double> temp_f0 = f0;
		std::string code;
		int octave;
		int sharp_flag = 0;
		int i;
		int j;

		pitch_bend_rate *= 1000;

		std::smatch match;

		std::regex reg("([A-Z])");
		if (regex_search(tool2_argments.note_num, match, reg))
		{
			code = match[0];
		}
		reg = ("(#)");
		if (regex_search(tool2_argments.note_num, match, reg))
		{
			sharp_flag++;
		}
		reg = ("(\\d+)");
		if (regex_search(tool2_argments.note_num, match, reg))
		{
			octave = stof(match[0]);
		}

		int note_diff = (code[0] - ascii_A);

		if (note_diff >= 2)
		{
			note_diff -= 7;
		}
		note_diff *= 2;
		if (note_diff <= -6)
		{
			note_diff++;
		}
		note_diff += sharp_flag;


		double reference_frq = 0;
		const double base_frq = (static_cast<int>(A1_Frq) << (octave - 1)) * pow(2, ((double)note_diff / 12));
		//UTAUではB0以下は指定不可のため

		f0.resize(tool2_argments.pitch_bend.size());

		if (tool2_argments.modulation != 0)
		{
			double weight;
			double sumWeight = 0;

			if (f0[0] != 0)
			{
				reference_frq = f0[0];
				sumWeight++;
			}

			for (i = 1; i < f0.size(); i++)
			{
				if (f0[i] != 0)
				{
					weight = 2 - std::min((double)(abs(f0[i - 1] - f0[i]) / f0[i]), 1.0);

					reference_frq += f0[i] * weight;
					sumWeight += weight;
				}
			}
			if (sumWeight > 0)
			{
				reference_frq /= sumWeight;
			}

		}

		double expansion_ratio = (double)frame_period / (pitch_bend_rate * tool2_argments.velocity);
		double lerp_rate;

		consonant_end = std::min((int)f0.size(), consonant_end);

		j = 0;
		for (i = 0; i < consonant_end; ++i)
		{
			for (; i > (j + 1) * expansion_ratio; j++);

			if (temp_f0[j] * temp_f0[j + 1] != 0)
			{
				lerp_rate = (i - j * expansion_ratio) / expansion_ratio; //(((j + 1) *  expansion_ratio) - j *  expansion_ratio) ==  expansion_ratio
				f0[i] = base_frq * pow(2, ((double)tool2_argments.pitch_bend[i] / 1200))
					+ (lerp(temp_f0[j], temp_f0[j + 1], lerp_rate) - reference_frq) * tool2_argments.modulation / 100;
			}
			else if (temp_f0[j] + temp_f0[j + 1] != 0)
			{
				f0[i] = base_frq * pow(2, ((double)tool2_argments.pitch_bend[i] / 1200))
					+ (temp_f0[j] + temp_f0[j + 1] - reference_frq) * tool2_argments.modulation / 100;
			}
		}

		const int offset = j;
		const double offset_length = j * expansion_ratio;
		expansion_ratio = (double)(f0.size() - 1 - consonant_end) / (temp_f0.size() - 2 - offset);

		j = 0;
		for (; i < f0.size(); i++)
		{
			for (; (i - consonant_end) > (j + 1) * expansion_ratio; j++);

			if (temp_f0[j + offset] * temp_f0[j + offset + 1] != 0)
			{
				lerp_rate = (i - offset_length - j * expansion_ratio) / expansion_ratio;
				f0[i] = base_frq * pow(2, ((double)tool2_argments.pitch_bend[i] / 1200))
					+ (lerp(temp_f0[j + offset], temp_f0[j + offset + 1], lerp_rate) - reference_frq) * tool2_argments.modulation / 100;
			}
			else if (temp_f0[j + offset] + temp_f0[j + offset + 1] != 0)
			{
				f0[i] = base_frq * pow(2, ((double)tool2_argments.pitch_bend[i] / 1200))
					+ (temp_f0[j] + temp_f0[j + 1] - reference_frq) * tool2_argments.modulation / 100;
			}
		}
	}


	void RelayToSynthesizer(int consonant_end, const double& pitch_bend_rate, const double& velocity,
		jewelwasp::VoiceParameters& voice_parameters, std::vector<double>& wdata)
	{
		WorldSynthesizer synthesizer = { 0 };
		std::vector<double*> p_spectrogram;
		std::vector<double*> p_aperiodicity;
		int i;
		int j;
		int offset;

		int f0_length = voice_parameters.f0[0].size();
		voice_parameters.f0[0].resize(voice_parameters.f0[0].size() + 50);//適当

		p_spectrogram.resize(voice_parameters.f0[0].size());
		p_aperiodicity.resize(voice_parameters.f0[0].size());

		double expansion_ratio = (double)voice_parameters.frame_period / (pitch_bend_rate * 1000 * velocity);

		consonant_end = std::min(f0_length, consonant_end);

		j = 0;
		for (i = 0; i < consonant_end; ++i)
		{
			for (; i > j * expansion_ratio; j++);

			p_spectrogram[i] = &voice_parameters.spectrogram[0][j][0];
			p_aperiodicity[i] = &voice_parameters.aperiodicity[0][j][0];
		}

		offset = j;
		expansion_ratio = (double)(f0_length - consonant_end - 1) / (voice_parameters.spectrogram[0].size() - 2 - offset);

		j = 0;
		for (; i < f0_length; ++i)
		{
			for (; i - consonant_end > j * expansion_ratio; j++);

			p_spectrogram[i] = &voice_parameters.spectrogram[0][j + offset][0];
			p_aperiodicity[i] = &voice_parameters.aperiodicity[0][j + offset][0];
		}

		for (; i < voice_parameters.f0[0].size(); ++i)
		{
			p_spectrogram[i] = p_spectrogram[f0_length - 1];
			p_aperiodicity[i] = p_aperiodicity[f0_length - 1];
			voice_parameters.f0[0][i] = voice_parameters.f0[0][f0_length - 1];
		}

		Synthesis(&voice_parameters.f0[0][0], voice_parameters.f0[0].size(), &p_spectrogram[0], &p_aperiodicity[0],
			voice_parameters.fft_size, pitch_bend_rate * 1000, voice_parameters.fs, wdata.size(), &wdata[0]);
	}

	void process_Volume(const int& consonant_end, std::vector<double>& wdata, const double& volume, const double& P_flag_num)
	{
		std::pair<double*, double*> wdata_minmax = std::minmax_element(&wdata[consonant_end], &wdata[wdata.size() - 1]);
		double tmpd[3];
		tmpd[0] = std::abs(*wdata_minmax.first);
		tmpd[1] = std::abs(*wdata_minmax.second);
		tmpd[2] = std::max(tmpd[0], tmpd[1]);

		double amplification_factor = volume / 100.0 * lerp(1.0, (0.5 / tmpd[2]), P_flag_num / 100.0);

		for (int i = 0; i < wdata.size(); i++)
		{
			wdata[i] *= amplification_factor;
		}
	}

}