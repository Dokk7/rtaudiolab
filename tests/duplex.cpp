/******************************************/
/*
  duplex.cpp
  by Gary P. Scavone, 2006-2007.

  This program opens a duplex stream and passes
  input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <string.h>

/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8


typedef signed short MY_TYPE;
#define FORMAT RTAUDIO_SINT16


typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32
*/
typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

struct MyData{

double* impRep;   				//réponse impulsionelle
double* interBuffer;			//buffer intermédiaire, résultat de la convolution
int M;										//M taille de impRep
int L;										//L taille du buffer
int bufferMax;						//taille maximale de la convolution
double tempsMax;					//temps maximal de calcul
int fs;										//fréquence
double* impRepReal;				//réponse impulsionelle fréquentielle réelle
double* impRepIm;					//réponse impulsionelle fréquentielle imaginaire
double* interBufferReal;	//buffer intermédiaire, résultat de la convolution fréquentielle
double* interBufferIm;		//buffer intermédiaire, résultat de la convolution fréquentielle
int sizeFFT;						//taille de la fft
int convChoice;						//0=none, 1=temp, 2=freq

};

/*****************************************************************************************************/
//Some functions

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

double *_sintbl = 0;
int maxfftsize = 0;
int fft(double *x, double *y, const int m);
  int ifft(double *x, double *y, const int m);
int fftr(double *x, double *y, const int m);
int ifftr(double *x, double *y, const int l);
  static int checkm(const int m);
int get_nextpow2(int n);
char *getmem(int leng, unsigned size);
double *dgetmem(int leng);
double get_process_time();


///////////////////////////////
// FFT functions
int get_nextpow2(int n)
{
  int k = 1;
  while (k < n){
    k *= 2;
  }

  return k;
}

int fftr(double *x, double *y, const int m)
{
   int i, j;
   double *xp, *yp, *xq;
   double *yq;
   int mv2, n, tblsize;
   double xt, yt, *sinp, *cosp;
   double arg;

   mv2 = m / 2;

   /* separate even and odd  */
   xq = xp = x;
   yp = y;
   for (i = mv2; --i >= 0;) {
      *xp++ = *xq++;
      *yp++ = *xq++;
   }

   if (fft(x, y, mv2) == -1)    /* m / 2 point fft */
      return (-1);


   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = M_PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }
   //printf("Debug: m=%i, maxfftsize=%i\n",m,maxfftsize);

   n = maxfftsize / m;
   sinp = _sintbl;
   cosp = _sintbl + maxfftsize / 4;

   xp = x;
   yp = y;
   xq = xp + m;
   yq = yp + m;
   *(xp + mv2) = *xp - *yp;
   *xp = *xp + *yp;
   *(yp + mv2) = *yp = 0;

   for (i = mv2, j = mv2 - 2; --i; j -= 2) {
      ++xp;
      ++yp;
      sinp += n;
      cosp += n;
      yt = *yp + *(yp + j);
      xt = *xp - *(xp + j);
      *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
      *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
   }

   xp = x + 1;
   yp = y + 1;
   xq = x + m;
   yq = y + m;

   for (i = mv2; --i;) {
      *xp++ = *(--xq);
      *yp++ = -(*(--yq));
   }

   return (0);
}

int ifftr(double *x, double *y, const int l)
{
   int i;
   double *xp, *yp;

   fftr(x, y, l);

   xp = x;
   yp = y;
   i = l;
   while (i--) {
      *xp++ /= l;
      *yp++ /= -l;
   }

   return (0);
}


static int checkm(const int m)
{
   int k;

   for (k = 4; k <= m; k <<= 1) {
      if (k == m)
         return (0);
   }
   fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n",m);

   return (-1);
}

int fft(double *x, double *y, const int m)
{
   int j, lmx, li;
   double *xp, *yp;
   double *sinp, *cosp;
   int lf, lix, tblsize;
   int mv2, mm1;
   double t1, t2;
   double arg;
   int checkm(const int);

   /**************
   * RADIX-2 FFT *
   **************/

   if (checkm(m))
      return (-1);

   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = M_PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }

   lf = maxfftsize / m;
   lmx = m;

   for (;;) {
      lix = lmx;
      lmx /= 2;
      if (lmx <= 1)
         break;
      sinp = _sintbl;
      cosp = _sintbl + maxfftsize / 4;
      for (j = 0; j < lmx; j++) {
         xp = &x[j];
         yp = &y[j];
         for (li = lix; li <= m; li += lix) {
            t1 = *(xp) - *(xp + lmx);
            t2 = *(yp) - *(yp + lmx);
            *(xp) += *(xp + lmx);
            *(yp) += *(yp + lmx);
            *(xp + lmx) = *cosp * t1 + *sinp * t2;
            *(yp + lmx) = *cosp * t2 - *sinp * t1;
            xp += lix;
            yp += lix;
         }
         sinp += lf;
         cosp += lf;
      }
      lf += lf;
   }

   xp = x;
   yp = y;
   for (li = m / 2; li--; xp += 2, yp += 2) {
      t1 = *(xp) - *(xp + 1);
      t2 = *(yp) - *(yp + 1);
      *(xp) += *(xp + 1);
      *(yp) += *(yp + 1);
      *(xp + 1) = t1;
      *(yp + 1) = t2;
   }

   /***************
   * bit reversal *
   ***************/
   j = 0;
   xp = x;
   yp = y;
   mv2 = m / 2;
   mm1 = m - 1;
   for (lmx = 0; lmx < mm1; lmx++) {
      if ((li = lmx - j) < 0) {
         t1 = *(xp);
         t2 = *(yp);
         *(xp) = *(xp + li);
         *(yp) = *(yp + li);
         *(xp + li) = t1;
         *(yp + li) = t2;
      }
      li = mv2;
      while (li <= j) {
         j -= li;
         li /= 2;
      }
      j += li;
      xp = x + j;
      yp = y + j;
   }

   return (0);
}

int ifft(double *x, double *y, const int m)
{
   int i;

   if (fft(y, x, m) == -1)
      return (-1);

   for (i = m; --i >= 0; ++x, ++y) {
      *x /= m;
      *y /= m;
   }

   return (0);
}

double *dgetmem(int leng)
{
    return ( (double *)getmem(leng, sizeof(double)) );
}

char *getmem(int leng, unsigned size)
{
    char *p = NULL;

    if ((p = (char *)calloc(leng, size)) == NULL){
        fprintf(stderr, "Memory allocation error !\n");
        exit(3);
    }
    return (p);
}

double get_process_time() {
    struct rusage usage;
    if( 0 == getrusage(RUSAGE_SELF, &usage) ) {
        return (double)(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) +
               (double)(usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1.0e6;
    }
    return 0;
}
/**********************************************************************************************************/

void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device to use (default = 0),\n";
  std::cout << "    oDevice = optional output device to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
  exit( 0 );
}

void nothingDone( double *out, double *in, MyData* myData){
	int i;
  for(i=0;i<myData->L-1;i++){
  	out[i] = in[i];
  }
}

void convTemps( double *out, double *in, MyData* myData){

	int n;
  for (n = 0; n< myData->L + myData->M -1; n++)
  {
  	if (n< myData->M -1){
			//décalage de interBuffer
			myData->interBuffer[n] = myData->interBuffer[n+myData->L];
		}
		else{
			//retirer la fin de interBuffer
			myData->interBuffer[n] = 0;
		}
		
  	if (n < myData->bufferMax){
    	int kmin, kmax, k;
    	//kmin = premier indice pour lequel le produit dans la convolution est non nul
    	//kmax = dernier indice pour lequel le produit dans la convolution est non nul
    	//k = indice de la somme

    	kmin = (n >= myData->M - 1) ? n - (myData->M - 1) : 0;  //kmin = 0 ou kmin < 512
    	kmax = (n < myData->L - 1) ? n : myData->L - 1;					//kmax < 512 ou kmax = 512

			double tmp;
			tmp = 0;
    	for (k = kmin; k <= kmax; k++)
    	{
      	tmp  += in[k]* myData->impRep[n - k];
    	}
    	myData->interBuffer[n] += tmp;
    	if(n<myData->L){
    		out[n] = myData->interBuffer[n];
    	}
    }
  }
}

void convFreq2(double *out, double *in, MyData* myData){
	myData->sizeFFT = get_nextpow2(myData->bufferMax);
	memcpy(myData->interBufferReal,in,myData->L);
	fftr(myData->interBufferReal,myData->interBufferIm,myData->sizeFFT);
	
	int n;
	for (n=0; n<myData->sizeFFT - myData->L;n++){
		myData->interBuffer[n] = myData->interBuffer[n+myData->L];
	}
	for(n=myData->sizeFFT - myData->L;n<myData->sizeFFT;n++){
		myData->interBuffer[n] = 0;
	}
	for (n=0; n<myData->sizeFFT; n++){
		myData->interBufferReal[n] = myData->interBufferReal[n]*myData->impRepReal[n] + myData->interBufferIm[n]*myData->impRepIm[n];
		myData->interBufferIm[n] = myData->interBufferReal[n]*myData->impRepIm[n] + myData->interBufferIm[n]*myData->impRepReal[n];
	}
	ifft(myData->interBufferReal,myData->interBufferIm,myData->sizeFFT);
	for (n=0; n<myData->sizeFFT;n++){
		myData->interBuffer[n] += myData->interBufferReal[n];
	}
	for (n=0; n<myData->L;n++){
		out[n] = myData->interBuffer[n];
	}
} 

void convFreq(double *out, double *in, MyData* myData){
	memset(myData->interBufferReal, 0, myData->sizeFFT);
	memset(myData->interBufferIm, 0, myData->sizeFFT);
	memcpy(myData->interBufferReal,in,myData->L);
	fftr(myData->interBufferReal,myData->interBufferIm,myData->sizeFFT);
	
	int n;
	for (n=0; n<myData->sizeFFT; n++){
		if (n<myData->sizeFFT - myData->L){
			myData->interBuffer[n] = myData->interBuffer[n+myData->L];
		}
		else{
			myData->interBuffer[n] = 0;
		}
		myData->interBufferReal[n] = myData->interBufferReal[n]*myData->impRepReal[n] - myData->interBufferIm[n]*myData->impRepIm[n];
		myData->interBufferIm[n] = myData->interBufferReal[n]*myData->impRepIm[n] + myData->interBufferIm[n]*myData->impRepReal[n];
	}
	ifft(myData->interBufferReal,myData->interBufferIm,myData->sizeFFT);
	for (n=0; n<myData->sizeFFT;n++){
		myData->interBuffer[n] += myData->interBufferReal[n];
	}
	for (n=0; n<myData->L;n++){
		out[n] = myData->interBuffer[n];
	}
}

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
	double time;
	double* in = (double*)inputBuffer;
	double* out = (double*)outputBuffer;
	
  time = get_process_time();
  
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  MyData* myData = (MyData *) data;
  
  switch(myData->convChoice)
	{
    case 0:
      nothingDone(out, in, myData);
      break;
    case 1:
      convTemps(out, in, myData);
      break;
   	case 2:
   		convFreq(out, in, myData);
   		break;
    default :
        printf("Unknown Error !!\n");
	}
  
  time = get_process_time()-time;
  
  if(time > myData->tempsMax){
  	myData->bufferMax -= 500;
 		printf("Time too long, new bufferSize : %d, ",myData->bufferMax);
  	printf("Time : %f\n",time);
  }
  outputBuffer = out;
  return 0;
}

int main( int argc, char *argv[] )
{
  unsigned int channels, fs, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

  // Minimal command-line checking
  if (argc < 3 || argc > 7 ) usage();

  RtAudio adac;
  if ( adac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  if ( argc > 3 )
    iDevice = (unsigned int) atoi(argv[3]);
  if ( argc > 4 )
    oDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    iOffset = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    oOffset = (unsigned int) atoi(argv[6]);

  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 512;
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = iDevice;
  iParams.nChannels = channels;
  iParams.firstChannel = iOffset;
  oParams.deviceId = oDevice;
  oParams.nChannels = channels;
  oParams.firstChannel = oOffset;

  if ( iDevice == 0 )
    iParams.deviceId = adac.getDefaultInputDevice();
  if ( oDevice == 0 )
    oParams.deviceId = adac.getDefaultOutputDevice();

  RtAudio::StreamOptions options;
  
  //ici on prend la réponse impulsionelle du filtre
  
  MyData myData;
  myData.convChoice = -1;
  
  while(myData.convChoice != 0 && myData.convChoice != 1 && myData.convChoice != 2){
		printf("Choose convolution : 0=none, 1=temp, 2=freq\n : ");
		scanf("%d",&myData.convChoice);
		printf("\n\n");
	}
  
	FILE* impRepFile;
	impRepFile = fopen("../../c/impres", "rb");
	
	fseek(impRepFile, 0, SEEK_END); // seek to end of file
	long unsigned int fileSize = ftell(impRepFile)/sizeof(double); // nombre d'elements de la RI
	
	fileSize = fileSize/2;
	
	fseek(impRepFile, 0, SEEK_SET); // seek back to beginning of file
	
	double impRep[fileSize];
	int sizeFFT = get_nextpow2(fileSize+bufferFrames-1);
	
	double impRepReal[sizeFFT];
	double impRepIm[sizeFFT];
	memcpy(impRepReal,impRep,fileSize);
	
	fread(impRep,sizeof(double), fileSize, impRepFile);
	printf("Taille du fichier %lu\n\n",fileSize);
	
	fclose(impRepFile);

	myData.sizeFFT = sizeFFT;
	myData.impRep = impRep;
	memcpy(impRepReal, impRep, fileSize);
	fft(impRepReal,impRepIm,myData.sizeFFT);
	myData.M = fileSize;
	myData.L = bufferFrames;
	myData.bufferMax = myData.M+myData.L-1;
	myData.tempsMax = (double)bufferFrames/(double)fs;
	myData.fs = fs;
	myData.impRepReal = impRepReal;
	myData.impRepIm = impRepIm;
	
	int k;
	printf("impRep : [");
	for (k = 0; k<300; k++){
		printf("%f, ",impRep[k]);
	}
	
	printf("...]\n\nimpRepReal : [");
	for (k = 0; k<300; k++){
		printf("%f, ",impRepReal[k]);
	}
	
	printf("...]\n\nimpRepIm : [");
	for (k = 0; k<300; k++){
		printf("%f, ",impRepIm[k]);
	}
	printf("...]\n\n");
	
	double interBuffer[sizeFFT];
	myData.interBuffer = interBuffer;
	
	double interBufferReal[sizeFFT];
	myData.interBufferReal = interBufferReal;
	
	double interBufferIm[sizeFFT];
	myData.interBufferIm = interBufferIm;
	
	double tempsMax = (double)bufferFrames/(double)fs;
	printf("Fréquence : %d, Temps max : %f, sizeFFT : %d\n\n",fs,tempsMax,sizeFFT);
	//fin
	
  //options.flags |= RTAUDIO_NONINTERLEAVED;

  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)&myData, &options );
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    exit( 1 );
  }

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  try {
    adac.startStream();
    char input;
    std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
    std::cin.get(input);
    std::cin.get(input);
    // Stop the stream.
    adac.stopStream();
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    goto cleanup;
  }

 cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  return 0;
}
