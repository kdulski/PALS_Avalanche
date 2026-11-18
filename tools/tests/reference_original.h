// Verbatim copy of the model functions of the original code (branch master,
// PALS_avalanche_fitFunctions.cpp) used as the reference in the tests.
#ifndef PALS_REFERENCE_ORIGINAL_H
#define PALS_REFERENCE_ORIGINAL_H
#include <cmath>
namespace original {
typedef double Double_t;
#define oPsLFLimit 0.7f
inline Double_t ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 );
inline Double_t GetIntensityParameterNew( Double_t *P, int Type, unsigned startIndex, unsigned NumberOfParameter );
inline Double_t DiscreteFitFunctionNoPs( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0.;
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 0; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
		}
		if( FixedIntensity >=1 )
		{
			return sum + P[0];
		}
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
						GetIntensityParameterNew( P, 0, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)], (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)])+1 )*(1 - FixedIntensity )/* Intensity */ );
		}	
		FixedIntensity = 0.;
	}
	return sum + P[0];
}
inline Double_t DiscreteFitFunctionPs( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreepPsIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > oPsLFLimit )
			FreepPsIntensity += P[5 + 3*(int)P[1]]*GetIntensityParameterNew( P, 2, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 2, (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]) );
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 1; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j++ )
		{
			if( P[4 + 3*(int)P[1] + 2*j] > oPsLFLimit )
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
												      ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
			}
			else
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
							GetIntensityParameterNew( P, 2, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 2, (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]) )*(1 - FixedIntensity )/(1 + FreepPsIntensity)/* Intensity */ );
			}
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 +
			FreepPsIntensity)/* Intensity */ );
		}
		FixedIntensity = 0.;
	}
	return sum + P[0];
}
inline Double_t ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 )			//A - arguments, p - parameters
{
	Double_t Y = P4 / (2*P3) * exp( P1*P1 / (2*P3*P3) - ( A - P2 ) / P3 ) * ( erf( ( A - P1*P1 / P3 - P2 ) / (sqrt(2)*P1) ) - erf( ( - P1*P1 / P3 - P2 ) / (sqrt(2)*P1) ) );
	return Y;
}
inline Double_t GetIntensityParameterNew( Double_t *P, int Type, unsigned startIndex, unsigned NumberOfParameter )
{
  	unsigned dimension = 0; 
	int type = Type;
	switch(type)
	{
		case 0:
			dimension = P[2] - P[6 + 3*(int)P[1] + 2*((int)P[2]-1)];			//Free Lifetime components parameters without p-Ps
			break;
		case 1:
			dimension = P[1];								//Resolution components (type is one because between different components there are three paramateres difference) And default difference between parameters in lifetimes is two -> two plus 1 is three (For sake not having negative parameters and to unify function)
			break;	
		case 2:
			dimension = P[2] - P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] - 1;			//Free Lifetime components parameters with p-Ps
			type = 0;
			break;
		case 3:											
			dimension = P[1];
			type = 0;
			break;
	}
	double intensity = 0.1;
	double Amplitude = 1.; //
	for( int i=dimension-2; i>0; i-- )
	{
		Amplitude = Amplitude*sin( P[startIndex + (2 + type)*i] ) + cos( P[startIndex + (2 + type)*i] );
	}
	switch(dimension)
	{
		case 0:
			intensity = 1.;
			break;
		case 1:
			intensity = 1.;
			break;
		case 2:
			if( NumberOfParameter == 1 )
				intensity = P[startIndex];
			else
				intensity = 1 - P[startIndex];
			break;
		case 3:
			if( NumberOfParameter == 1 )
				intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;		// r cos(Phi) / Amplitude
			else if( NumberOfParameter == 2 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)]) /Amplitude; // r sin(Phi) / Amplitude
			else
				intensity = 1 - P[startIndex];				// 1 - r
			break;
		case 4:
			if( NumberOfParameter == 1 )
				intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
			else if( NumberOfParameter == 2 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*sin(Phi2) / Amplitude
			else if( NumberOfParameter == 3 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 5:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 6:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4]) /Amplitude;
			else
				intensity = 1 - P[startIndex + (2 + type)];
			break;
		case 7:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 8:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*cos(P[startIndex + (2 + type)*6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 9:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*cos(P[startIndex + (2 + type)*6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*cos(P[startIndex + (2 + type)*7]) /Amplitude;
			else if( NumberOfParameter == 8 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*sin(P[startIndex + (2 + type)*7]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 10:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*cos(P[startIndex + (2 + type)*6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*cos(P[startIndex + (2 + type)*7]) /Amplitude;
			else if( NumberOfParameter == 8 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*sin(P[startIndex + (2 + type)*7])*cos(P[startIndex + (2 + type)*8]) /Amplitude;
			else if( NumberOfParameter == 9 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*sin(P[startIndex + (2 + type)*7])*sin(P[startIndex + (2 + type)*8]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
	}
	return (Double_t)intensity;
}
#undef oPsLFLimit
} // namespace original
#endif
