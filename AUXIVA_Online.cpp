#include <stdio.h>
#include "AUXIVA_Online.h"
#include "header.h"
#include "sigproc.h"


AUXIVA::AUXIVA()
{
	nfft = nWin;
	nshift = nWin / 4;
	nol = 3 * nWin / 4;
	nfreq = nfft / 2 + 1;
	//epsi = 0.000001;
	epsi = 2.220446049250313*1E-16;
	f_alpha = 0.96;

	int i, j, k, freq, ch;
	int re, im;
	X = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		X[i] = new double[nfreq * 2];
	}
	X_r = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		X_r[i] = new double[nfreq * 2];
	}
	Y = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		Y[i] = new double[nfreq * 2];
	}
	Pwr = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		Pwr[i] = new double[nfreq];
	}
	invWDE = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		invWDE[i] = new double[nfreq * 2];
	}
	sum_Pwr = new double[Nch];
	r = new double[Nch];
	diag_WV = new double*[Nch];
	for (i = 0; i < Nch; i++)
	{
		diag_WV[i] = new double[nfreq * 2];
		sum_Pwr[i] = 0.0;
		r[i] = 0.0;
	}
	win_STFT = new double[nWin];
	for (i = 0; i < nWin; i++)
	{
		win_STFT[i] = sqrt((double)2/3) * 0.5 * (1.0 - cos(2.0 * (double)M_PI*(double)(i) / (nWin)));
	}
	W = new double **[Nch];
	for (i = 0; i < Nch; i++)
	{
		W[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			W[i][j] = new double[nfreq * 2];		
		}
	}
	V = new double ***[Nch];
	for (i = 0; i < Nch; i++)
	{
		V[i] = new double**[Nch];
		for (j = 0; j < Nch; j++)
		{		
			V[i][j] = new double*[nfreq * 2];
			for (k = 0; k < nfreq * 2; k++)
			{
				V[i][j][k] = new double[Nch];
			}
		}
	}
	U = new double ***[Nch];
	for (i = 0; i < Nch; i++)
	{
		U[i] = new double**[Nch];
		for (j = 0; j < Nch; j++)
		{
			U[i][j] = new double*[nfreq * 2];
			for (k = 0; k < nfreq * 2; k++)
			{
				U[i][j][k] = new double[Nch];
			}
		}
	}

	//W Á¤ÀÇ
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				if (i == j)
				{
					W[i][j][re] = 1.0;
					W[i][j][im] = 0.0;
				}
				else
				{
					W[i][j][re] = 0.0;
					W[i][j][im] = 0.0;
				}
			}
		}
	}

	//frameInd over 2
	p = new double[Nch];
	Unumer = new double[nfreq * 2];
	Udenom = new double**[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		Udenom[ch] = new double*[Nch];
		for (i = 0; i < Nch; i++)
		{
			Udenom[ch][i] = new double[nfreq * 2];
		}		
	}
	p_U_X = new double*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		p_U_X[ch] = new double[nfreq * 2];
	}
	X_T_U = new double*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		X_T_U[ch] = new double[nfreq * 2];
	}
	p_U_X_X = new double***[Nch];
	for ( ch = 0; ch < Nch; ch++)
	{
		p_U_X_X[ch] = new double**[Nch];
		for ( i = 0; i < Nch; i++)
		{
			p_U_X_X[ch][i] = new double*[nfreq * 2];
			for ( j = 0; j < nfreq * 2; j++)
			{
				p_U_X_X[ch][i][j] = new double[Nch];
			}
		}
	}

	//normalizing
	normCoef = new double[nfreq * 2];
	sqnorm = new double[nfreq * 2];
	A = new double**[Nch];
	for (i = 0; i < Nch; i++)
	{
		A[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			A[i][j] = new double[nfreq * 2];
		}
	}
	WDE_V = new double*[Nch];
	for ( i = 0; i < Nch; i++)
	{
		WDE_V[i] = new double[nfreq * 2];
	}
	w = new double*[Nch];
	for ( i = 0; i < Nch; i++)
	{
		w[i] = new double[nfreq * 2];
	}
	dW = new double**[Nch];
	for ( i = 0; i < Nch; i++)
	{
		dW[i] = new double*[Nch];
		for ( j = 0; j < Nch; j++)
		{
			dW[i][j] = new double[nfreq * 2];
		}
	}

	//Calculate A
	Anumer = new double[nfreq * 2];
	AdW = new double**[Nch];
	for (i = 0; i < Nch; i++)
	{
		AdW[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			AdW[i][j] = new double[nfreq * 2];
		}
	}
	Adenom = new double**[Nch];
	for ( i = 0; i < Nch; i++)
	{
		Adenom[i] = new double*[Nch];
		for ( j = 0; j < Nch; j++)
		{
			Adenom[i][j] = new double[nfreq * 2];
		}
	}

	//Result
	Wbp = new double**[Nch];
	for (i = 0; i < Nch; i++)
	{
		Wbp[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			Wbp[i][j] = new double[nfreq * 2]; 
		}
	}
	Ytmp = new double*[Nch];
	for (i = 0; i < Nch; i++)
	{
		Ytmp[i] = new double[nfreq * 2];
	}
	Ybuff = new double*[Nch];
	for ( i = 0; i < Nch; i++)
	{
		Ybuff[i] = new double[nWin];
	}
	int sample;
	for (i = 0; i < Nch; i++)
	{
		for (sample = 0; sample < nfreq * 2; sample++)
		{
			Ytmp[i][sample] = 0.0;
		}
		for (sample = 0; sample < nWin; sample++)
		{
			Ybuff[i][sample] = 0.0;
		}
	}

}

AUXIVA::~AUXIVA()
{
	int i, j, k;
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (k = 0; k < nfreq * 2; k++)
			{
				delete[] V[i][j][k];
				delete[] U[i][j][k];
     			delete[] p_U_X_X[i][j][k];
			}
			delete[] V[i][j];
			delete[] U[i][j];
			delete[] p_U_X_X[i][j];
		}
		delete[] V[i];
		delete[] U[i];
		delete[] p_U_X_X[i];
	}
	delete[] V;
	delete[] U;
	delete[] p_U_X_X;


	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] W[i][j];
		}
		delete[] X[i];
		delete[] X_r[i];
		delete[] Y[i];
		delete[] Pwr[i];
		delete[] W[i];
		delete[] invWDE[i];
		delete[] diag_WV[i];
	}
	delete[] X;
	delete[] X_r;
	delete[] Y;
	delete[] Pwr;
	delete[] sum_Pwr;
	delete[] r;
	delete[] W;
	delete[] invWDE;
	delete[] diag_WV;
	delete[] win_STFT;

	//frameInd over 2
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] Udenom[i][j];
		}
		delete[] Udenom[i];
		delete[] p_U_X[i];
		delete[] X_T_U[i];
	}
	delete[] Udenom;
	delete[] p_U_X;
	delete[] X_T_U;
	delete[] p;
	delete[] Unumer;

	//normalizing
	delete[] normCoef;
	delete[] sqnorm;
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] A[i][j];
			delete[] dW[i][j];
		}
		delete[] A[i];
		delete[] dW[i];
		delete[] WDE_V[i];
		delete[] w[i];
	}
	delete[] A;
	delete[] dW;
	delete[] WDE_V;
	delete[] w;

	//Calculate A
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] AdW[i][j];
			delete[] Adenom[i][j];
		}
		delete[] AdW[i];
		delete[] Adenom[i];
	}
	delete[] AdW;
	delete[] Adenom;
	delete[] Anumer;

	//result
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] Wbp[i][j];
		}
		delete[] Wbp[i];
		delete[] Ytmp[i];
		delete[] Ybuff[i];
	}
	delete[] Wbp;
	delete[] Ytmp;
	delete[] Ybuff;
}

void AUXIVA::AUXIVA_lemma(double **input, int frameInd, double **output)
{
	int i, j, ch, channel, freq, freqInd;
	int ch1, ch2;
	int re, im;
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < nWin; i++)
		{
			X[ch][i] = win_STFT[i] * input[ch][i];
		}
		hfft3(X[ch], nfft, 1);
	}

	for (freq = 0; freq < nfreq; freq++)
	{
		re = freq + freq;
		im = re + 1;
		for ( ch1 = 0; ch1 < Nch; ch1++)
		{
			Y[ch1][re] = W[ch1][0][re] * X[0][re] - W[ch1][0][im] * X[0][im] + W[ch1][1][re] * X[1][re] - W[ch1][1][im] * X[1][im];
			Y[ch1][im] = W[ch1][0][re] * X[0][im] + W[ch1][0][im] * X[0][re] + W[ch1][1][re] * X[1][im] + W[ch1][1][im] * X[1][re];
		}
	}

	// Pwr
	for (i = 0; i < Nch; i++)
	{
		sum_Pwr[i] = 0.0;
		for (j = 0; j < nfreq; j++)
		{
			re = j + j;
			im = re + 1;
			Pwr[i][j] = pow(Y[i][re], 2) + pow(Y[i][im], 2);
			if (Pwr[i][j] < epsi)
			{
				Pwr[i][j] = epsi;
			}
			sum_Pwr[i] = sum_Pwr[i] + Pwr[i][j];
		}
	}

	//r & p
	for (i = 0; i < Nch; i++)
	{
		r[i] = sqrt(sum_Pwr[i]);
		p[i] = (1 - f_alpha) / r[i];
	}


	for (ch = 0; ch < Nch; ch++)
	{
		if (frameInd == 3)
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd + freqInd;
				im = re + 1;
				for (freq = 0; freq < nfreq * 2; freq++)
				{
					for (channel = 0; channel < Nch; channel++)
					{
						X_r[channel][freq] = X[channel][freq] / r[ch];
					}
				}
				// Calculate V
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						V[ch1][ch2][re][ch] = X_r[ch1][re] * X[ch2][re] + X_r[ch1][im] * X[ch2][im];
						V[ch1][ch2][im][ch] = X_r[ch1][im] * X[ch2][re] - X_r[ch1][re] * X[ch2][im];
					}
				}
				// Calculate diag_WV
				for (channel = 0; channel < Nch; channel++)
				{
					diag_WV[channel][re] = W[channel][0][re] * V[0][channel][re][ch] - W[channel][0][im] * V[0][channel][im][ch] + W[channel][1][re] * V[1][channel][re][ch] - W[channel][1][im] * V[1][channel][im][ch];
					diag_WV[channel][im] = W[channel][0][re] * V[0][channel][im][ch] + W[channel][0][im] * V[0][channel][re][ch] + W[channel][1][re] * V[1][channel][im][ch] - W[channel][1][im] * V[1][channel][re][ch];
				}

				// Calculate inverse diag_WV
				if (ch == 0)
				{
					invWDE[0][re] = diag_WV[ch][re] / (pow(diag_WV[ch][re], 2) + pow(diag_WV[ch][im], 2));
					invWDE[0][im] = -diag_WV[ch][im] / (pow(diag_WV[ch][re], 2) + pow(diag_WV[ch][im], 2));
					invWDE[1][re] = 0.0;
					invWDE[1][im] = 0.0;
				}
				else
				{
					invWDE[0][re] = 0.0;
					invWDE[0][im] = 0.0;
					invWDE[1][re] = diag_WV[ch][re] / (pow(diag_WV[ch][re], 2) + pow(diag_WV[ch][im], 2));
					invWDE[1][im] = -diag_WV[ch][im] / (pow(diag_WV[ch][re], 2) + pow(diag_WV[ch][im], 2));
				}

				// Calculate U
				for ( ch1 = 0; ch1 < Nch; ch1++)
				{
					for ( ch2 = 0; ch2 < Nch; ch2++)
					{
						if (ch1 == ch2)
						{
							U[ch1][ch2][re][ch] = V[ch1][ch2][re][ch] / (pow(V[ch1][ch2][re][ch],2) + pow(V[ch1][ch2][im][ch], 2));
							U[ch1][ch2][im][ch] = -V[ch1][ch2][im][ch] / (pow(V[ch1][ch2][re][ch], 2) + pow(V[ch1][ch2][im][ch], 2));
						}
						else
						{
							U[ch1][ch2][re][ch] = 0.0;
							U[ch1][ch2][im][ch] = 0.0;
						}
					}
				}
			}
		}
		else
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd + freqInd;
				im = re + 1;
				// Calculate p_U_X
				for (channel = 0; channel < Nch; channel++)
				{
					p_U_X[channel][re] = p[ch] * (U[channel][0][re][ch] * X[0][re] - U[channel][0][im][ch] * X[0][im] + U[channel][1][re][ch] * X[1][re] - U[channel][1][im][ch] * X[1][im]);
					p_U_X[channel][im] = p[ch] * (U[channel][0][im][ch] * X[0][re] + U[channel][0][re][ch] * X[0][im] + U[channel][1][re][ch] * X[1][im] + U[channel][1][im][ch] * X[1][re]);
				}
				for ( ch1 = 0; ch1 < Nch; ch1++)
				{
					for ( ch2 = 0; ch2 < Nch; ch2++)
					{
						p_U_X_X[ch1][ch2][re][ch] = p_U_X[ch1][re] * X[ch2][re] + p_U_X[ch1][im] * X[ch2][im];
						p_U_X_X[ch1][ch2][im][ch] = p_U_X[ch1][im] * X[ch2][re] - p_U_X[ch1][re] * X[ch2][im];
					}
				}
				for ( ch1 = 0; ch1 < Nch; ch1++)
				{
					for ( ch2 = 0; ch2 < Nch; ch2++)
					{
						Udenom[ch1][ch2][re] = p_U_X_X[ch1][0][re][ch] * U[ch2][0][re][ch] + p_U_X_X[ch1][0][im][ch] * U[ch2][0][im][ch] + p_U_X_X[ch1][1][re][ch] * U[ch2][1][re][ch] + p_U_X_X[ch1][1][im][ch] * U[ch2][1][im][ch];
						Udenom[ch1][ch2][im] = p_U_X_X[ch1][0][im][ch] * U[ch2][0][re][ch] - p_U_X_X[ch1][0][re][ch] * U[ch2][0][im][ch] + p_U_X_X[ch1][1][im][ch] * U[ch2][1][re][ch] - p_U_X_X[ch1][1][re][ch] * U[ch2][1][im][ch];
					}
				}
				for ( channel = 0; channel < Nch; channel++)
				{
					X_T_U[channel][re] = X[0][re] * U[0][channel][re][ch] + X[0][im] * U[0][channel][im][ch] + X[1][re] * U[1][channel][re][ch] + X[1][im] * U[1][channel][im][ch];
					X_T_U[channel][im] = X[0][re] * U[0][channel][im][ch] - X[0][im] * U[0][channel][re][ch] + X[1][re] * U[1][channel][im][ch] - X[1][im] * U[1][channel][re][ch];
				}

				Unumer[re] = pow(f_alpha, 2) + (f_alpha * p[ch]) * (X_T_U[0][re] * X[0][re] - X_T_U[0][im] * X[0][im] + X_T_U[1][re] * X[1][re] - X_T_U[1][im] * X[1][im]);
				Unumer[im] = (f_alpha * p[ch]) * (X_T_U[0][re] * X[0][im] + X_T_U[0][im] * X[0][re] + X_T_U[1][re] * X[1][im] + X_T_U[1][im] * X[1][re]);
				
				if (sqrt(pow(Unumer[re],2) + pow(Unumer[im], 2)) < epsi)
				{
					Unumer[re] = epsi;
					Unumer[im] = 0.0;
				}

				//Calculate U
				for ( ch1 = 0; ch1 < Nch; ch1++)
				{
					for ( ch2 = 0; ch2 < Nch; ch2++)
					{
						U[ch1][ch2][re][ch] = (U[ch1][ch2][re][ch] / f_alpha) - (Udenom[ch1][ch2][re] * Unumer[re] + Udenom[ch1][ch2][im] * Unumer[im]) / (pow(Unumer[re], 2) + pow(Unumer[im], 2));
						U[ch1][ch2][im][ch] = (U[ch1][ch2][im][ch] / f_alpha) - (Udenom[ch1][ch2][im] * Unumer[re] - Udenom[ch1][ch2][re] * Unumer[im]) / (pow(Unumer[re], 2) + pow(Unumer[im], 2));
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						U[ch1][ch2][re][ch] = (U[ch1][ch2][re][ch] + U[ch2][ch1][re][ch]) / 2;
						U[ch1][ch2][im][ch] = (U[ch1][ch2][im][ch] - U[ch2][ch1][im][ch]) / 2;
					}
				}

				//Calculate invWD
				if (ch == 0)
				{
					for ( ch1 = 0; ch1 < Nch; ch1++)
					{
						invWDE[ch1][re] = U[ch1][0][re][ch] * A[0][0][re] - U[ch1][0][im][ch] * A[0][0][im] + U[ch1][1][re][ch] * A[1][0][re] - U[ch1][1][im][ch] * A[1][0][im];
						invWDE[ch1][im] = U[ch1][0][im][ch] * A[0][0][re] + U[ch1][0][re][ch] * A[0][0][im] + U[ch1][1][im][ch] * A[1][0][re] + U[ch1][1][re][ch] * A[1][0][im];
					}
				}
				else
				{
					for (ch1 = 0; ch1 < Nch; ch1++)
					{
						invWDE[ch1][re] = U[ch1][0][re][ch] * A[0][1][re] - U[ch1][0][im][ch] * A[0][1][im] + U[ch1][1][re][ch] * A[1][1][re] - U[ch1][1][im][ch] * A[1][1][im];
						invWDE[ch1][im] = U[ch1][0][im][ch] * A[0][1][re] + U[ch1][0][re][ch] * A[0][1][im] + U[ch1][1][im][ch] * A[1][1][re] + U[ch1][1][re][ch] * A[1][1][im];
					}
				}
			}
		}
		// Normalizing
		for (freqInd = 0; freqInd < nfreq; freqInd++)
		{
			re = freqInd * 2;
			im = re + 1;
			for ( ch1 = 0; ch1 < Nch; ch1++)
			{
				WDE_V[ch1][re] = invWDE[0][re] * V[0][ch1][re][ch] + invWDE[0][im] * V[0][ch1][im][ch] + invWDE[1][re] * V[1][ch1][re][ch] + invWDE[1][im] * V[1][ch1][im][ch];
				WDE_V[ch1][im] = invWDE[0][re] * V[0][ch1][im][ch] - invWDE[0][im] * V[0][ch1][re][ch] + invWDE[1][re] * V[1][ch1][im][ch] - invWDE[1][im] * V[1][ch1][re][ch];
			}

			normCoef[re] = WDE_V[0][re] * invWDE[0][re] - WDE_V[0][im] * invWDE[0][im] + WDE_V[1][re] * invWDE[1][re] - WDE_V[1][im] * invWDE[1][im];
			normCoef[im] = WDE_V[0][re] * invWDE[0][im] + WDE_V[0][im] * invWDE[0][re] + WDE_V[1][re] * invWDE[1][im] + WDE_V[1][im] * invWDE[1][re];
			
			sqnorm[re] = sqrt(normCoef[re]);
			sqnorm[im] = 0.0;
			if (sqnorm[re] < epsi)
			{
				sqnorm[re] = epsi;
			}

			for ( ch1 = 0; ch1 < Nch; ch1++)
			{
				w[ch1][re] = invWDE[ch1][re] / sqnorm[re];
				w[ch1][im] = invWDE[ch1][im] / sqnorm[re];
			}
		}

		for ( freq = 0; freq < nfreq; freq++)
		{
			re = freq * 2;
			im = re + 1;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				dW[ch][ch1][re] = w[ch1][re] - W[ch][ch1][re];
				dW[ch][ch1][im] = -w[ch1][im] - W[ch][ch1][im];
				W[ch][ch1][re] = w[ch1][re];
				W[ch][ch1][im] = -w[ch1][im];
			}
		}

		// Calculate A
		if (frameInd == 3)
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd * 2;
				im = re + 1;
				for ( ch1 = 0; ch1 < Nch; ch1++)
				{
					for ( ch2 = 0; ch2 < Nch; ch2++)
					{
						if (ch1 == ch2)
						{
							A[ch1][ch2][re] = W[ch1][ch2][re] / (pow(W[ch1][ch2][re], 2) + pow(W[ch1][ch2][im], 2));
							A[ch1][ch2][im] = -W[ch1][ch2][im] / (pow(W[ch1][ch2][re], 2) + pow(W[ch1][ch2][im], 2));
						}
						else
						{
							A[ch1][ch2][re] = 0.0;
							A[ch1][ch2][im] = 0.0;
						}
					}
				}
			}
		}
		else
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd * 2;
				im = re + 1;
				for ( ch1 = 0; ch1 < Nch; ch1++)
				{
					for ( ch2 = 0; ch2 < Nch; ch2++)
					{
						AdW[ch1][ch2][re] = A[ch1][ch][re] * dW[ch][ch2][re] - A[ch1][ch][im] * dW[ch][ch2][im];
						AdW[ch1][ch2][im] = A[ch1][ch][re] * dW[ch][ch2][im] + A[ch1][ch][im] * dW[ch][ch2][re];
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						Adenom[ch1][ch2][re] = AdW[ch1][0][re] * A[0][ch2][re] - AdW[ch1][0][im] * A[0][ch2][im] + AdW[ch1][1][re] * A[1][ch2][re] - AdW[ch1][1][im] * A[1][ch2][im];
						Adenom[ch1][ch2][im] = AdW[ch1][0][re] * A[0][ch2][im] + AdW[ch1][0][im] * A[0][ch2][re] + AdW[ch1][1][re] * A[1][ch2][im] + AdW[ch1][1][im] * A[1][ch2][re];
					}
				}
				Anumer[re] = 1 + dW[ch][0][re] * A[0][ch][re] - dW[ch][0][im] * A[0][ch][im] + dW[ch][1][re] * A[1][ch][re] - dW[ch][1][im] * A[1][ch][im];
				Anumer[im] = dW[ch][0][re] * A[0][ch][im] + dW[ch][0][im] * A[0][ch][re] + dW[ch][1][re] * A[1][ch][im] + dW[ch][1][im] * A[1][ch][re];
				if (sqrt(pow(Anumer[re],2) + pow(Anumer[im], 2)) < epsi)
				{
					Anumer[re] = epsi;
					Anumer[im] = 0.0;
				}

				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						A[ch1][ch2][re] = A[ch1][ch2][re] - (Adenom[ch1][ch2][re] * Anumer[re] + Adenom[ch1][ch2][im] * Anumer[im]) / (pow(Anumer[re], 2) + pow(Anumer[im], 2));
						A[ch1][ch2][im] = A[ch1][ch2][im] - (Adenom[ch1][ch2][im] * Anumer[re] - Adenom[ch1][ch2][re] * Anumer[im]) / (pow(Anumer[re], 2) + pow(Anumer[im], 2));
					}
				}

			}
		}
	}

	// result - back projection using A
	for ( freqInd = 0; freqInd < nfreq; freqInd++)
	{
		re = freqInd * 2;
		im = re + 1;
		for ( ch1 = 0; ch1 < Nch; ch1++)
		{
			for ( ch2 = 0; ch2 < Nch; ch2++)
			{
				Wbp[ch1][ch2][re] = A[0][ch1][re] * W[ch1][ch2][re] - A[0][ch1][im] * W[ch1][ch2][im];
				Wbp[ch1][ch2][im] = A[0][ch1][re] * W[ch1][ch2][im] + A[0][ch1][im] * W[ch1][ch2][re];
			}
		}
		
		for ( ch1 = 0; ch1 < Nch; ch1++)
		{
			Ytmp[ch1][re] = Wbp[ch1][0][re] * X[0][re] - Wbp[ch1][0][im] * X[0][im] + Wbp[ch1][1][re] * X[1][re] - Wbp[ch1][1][im] * X[1][im];
			Ytmp[ch1][im] = Wbp[ch1][0][re] * X[0][im] + Wbp[ch1][0][im] * X[0][re] + Wbp[ch1][1][re] * X[1][im] + Wbp[ch1][1][im] * X[1][re];
		}
	}

	for (ch1 = 0; ch1 < Nch; ch1++)
	{
		hfft3(Ytmp[ch1], nfft, -1);
		for (i = 0; i < nWin - BufferSize; i++)
		{
			Ybuff[ch1][i] = Ybuff[ch1][BufferSize + i];
			Ybuff[ch1][i] += Ytmp[ch1][i];
		}
		for (; i < nWin; i++)
		{
			Ybuff[ch1][i] = Ytmp[ch1][i];
		}
	}

	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			output[ch][i] = Ybuff[ch][i];
		}
	}
}