#define Nch			2
#define nWin		1024
#define BufferSize		256
#define SamplingFreq	16000


class AUXIVA {

private:
	int nfft;
	int nshift;
	int nol;
	int nfreq;
	double epsi;
	double f_alpha;

	double *win_STFT;
	double **X;
	double **X_r;
	double **Y;
	double ***W;
	double **Pwr;
	double *sum_Pwr;
	double *r;
	double ****V;
	double ****U;
	double **diag_WV;
	double **invWDE;

	//frameInd over 2
	double *p;
	double **p_U_X;
	double ****p_U_X_X;
	double ***Udenom;
	double *Unumer;
	double **X_T_U;

	//normalizing
	double *normCoef;
	double *sqnorm;
	double ***A;
	double **WDE_V;
	double unW;
	double **w;
	double ***dW;

	//Calculate A
	double ***AdW;
	double ***Adenom;
	double *Anumer;

	//result
	double ***Wbp;
	double **Ytmp;
	double **Ybuff;

public:
	AUXIVA();
	~AUXIVA();
	void AUXIVA_lemma(double **input, int frameInd, double **output);
};
