#include "PALS_avalanche_deconvolution.h"

DCV::DCV()
{
    
};

double GaussDistr( double X, double Mean, double Sigma )
{
	return 1/sqrt(2*M_PI)/Sigma*exp(-pow(X-Mean, 2)/2/pow(Sigma, 2));
}

double LogGaussDistr( double X, double Mean, double Sigma )
{
	return 1/sqrt(2*M_PI)/Sigma/X*exp(-pow( log(X) - Mean, 2)/2/pow(Sigma, 2));
}

double DCV::FindingOptimalAreaParameter( std::vector< double > Sig )
{
	std::cout << "Area Parameter Finding" << std::endl;
	
    AreaParameter = 1;
    SSig = Sig;
	double sum = 0;
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
			sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	double MeanDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
		MeanDiff = (Values[i] - Background)/(Model(i)*Values.size());
	}
	std::cout << "Mean Difference between histogram and model " << MeanDiff << std::endl;
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	AreaParameter = MeanDiff;
	unsigned NmbOfPointToMinimize = 100;
	std::vector<double> ChiSqVector, AreaParameters, ChiSqDiffVector;
    double AreaParameterTemp = 1;
	for( unsigned k=0; k<NmbOfPointToMinimize; k++ )
	{
        AreaParameterTemp = 1 + k*AreaParameter/NmbOfPointToMinimize;//AreaParameter - ((int)NmbOfPointToMinimize/2 - (int)k)*AreaParameter/NmbOfPointToMinimize;
		AreaParameters.push_back( AreaParameterTemp );
		for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
		{
			sum = 0;
			for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
			{
				sum += AreaParameterTemp*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
			}
			Intensities(i) = sum;		//Intensities for calculations in this procedure
		}
		Model =  GradMax.t() * Intensities;
		double ChiSq2 = 0;
		for( unsigned i=0; i<Values.size(); i++ )
		{		
			Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
			ChiSq2 += pow( Diff(i), 2 ) / Values[i];
		}
		ChiSqVector.push_back(ChiSq2);
		if( k>0 )
			ChiSqDiffVector.push_back( ChiSqVector[k] - ChiSqVector[k-1] );
	}
  /*  for( unsigned i=0; i<ChiSqDiffVector.size(); i++ )
    {
        std::cout << AreaParameters[i] << " " << ChiSqVector[i] << " " << ChiSqDiffVector[i] << std::endl;
    }*/
    double OptimalParameterTest = FindPeak( AreaParameters, ChiSqDiffVector, 0, ChiSqDiffVector.size() );
	std::cout << "Area Parameter Minimum From ChiSq " << OptimalParameterTest << std::endl;
    AreaParameter = OptimalParameterTest;
    return OptimalParameterTest;
    
	/*unsigned MinElement = std::distance( std::begin(ChiSqVector), std::min_element( std::begin(ChiSqVector), std::end(ChiSqVector) ) );
	std::cout << "Area Parameter Minimal Element " << MinElement << std::endl;
	unsigned LeftSide = (MinElement > 5 ? MinElement - 5 : 0);
	unsigned RightSide = (MinElement < ChiSqVector.size() - 5 ? MinElement + 5 : ChiSqVector.size() - 1);
	double OptimalAreaParameter = FindPeak( AreaParameters, ChiSqDiffVector, LeftSide, RightSide ) + AreaParameter/100;
	std::cout << "Area Parameter Minimum From ChiSq " << OptimalAreaParameter << std::endl;
    AreaParameter = OptimalAreaParameter;
	return OptimalAreaParameter;*/
}

double DCV::FindingOptimalSigmaParameter( std::vector< double > Sig, std::string ShapeOfComponent )
{
	std::cout << "Sigma Parameter Finding" << std::endl;
    
    SSig = Sig;
	double sum = 0;
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
            if( ShapeOfComponent == "Gaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
                //Value following Gaussian distribution
            }
            else if( ShapeOfComponent == "LogGaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j]);
            }
            else
            {
                if( LifetimesFromDiscrete[j].Parameter < oPsLFLimit )
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
                }
                else		//Option for assuming logGaussian distributions
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j]);
                }
            }
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	double MeanDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
		MeanDiff = (Values[i] - Background)/(Model(i)*Values.size());
	}
	std::cout << "Mean Difference between histogram and model " << MeanDiff << std::endl;
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	unsigned NmbOfPointToMinimize = 100;
	std::vector<double> ChiSqVector, SigmaParameters, ChiSqDiffVector;
    double SigmaParametersTemp = 1;
	for( unsigned k=0; k<NmbOfPointToMinimize; k++ )
	{
        SigmaParametersTemp = 1.25 + 0.75*((int)k - (int)NmbOfPointToMinimize)/NmbOfPointToMinimize;//AreaParameter - ((int)NmbOfPointToMinimize/2 - (int)k)*AreaParameter/NmbOfPointToMinimize;
		SigmaParameters.push_back( SigmaParametersTemp );
		for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
		{
			sum = 0;
			for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
			{
                double SigTemp = SigmaParametersTemp*Sig[j];
                if( ShapeOfComponent == "Gaussian" )
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigTemp);
                    //Value following Gaussian distribution
                }
                else if( ShapeOfComponent == "LogGaussian" )
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigTemp,2)/2, SigTemp);
                }
                else
                {
                    if( LifetimesFromDiscrete[j].Parameter < oPsLFLimit )
                    {
                        sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigTemp);
                    }
                    else		//Option for assuming logGaussian distributions
                    {
                        sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigTemp,2)/2, SigTemp);
                    }
                }
            }
			Intensities(i) = sum;		//Intensities for calculations in this procedure
		}
		Model =  GradMax.t() * Intensities;
		double ChiSq2 = 0;
		for( unsigned i=0; i<Values.size(); i++ )
		{		
			Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
			ChiSq2 += pow( Diff(i), 2 ) / Values[i];
		}
		ChiSqVector.push_back(ChiSq2);
		if( k>0 )
			ChiSqDiffVector.push_back( ChiSqVector[k] - ChiSqVector[k-1] );
	}
	/*for( unsigned i=0; i<ChiSqDiffVector.size(); i++ )
    {
        std::cout << SigmaParameters[i] << " " << ChiSqVector[i] << " " << ChiSqDiffVector[i] << std::endl;
    }*/
    
	unsigned MinElement = std::distance( std::begin(ChiSqVector), std::min_element( std::begin(ChiSqVector), std::end(ChiSqVector) ) );
	std::cout << "Sigma Parameter Minimal Element Rough Esitmation " << MinElement << std::endl;
	unsigned LeftSide = (MinElement > 5 ? MinElement - 5 : 0);
	unsigned RightSide = (MinElement < ChiSqVector.size() - 5 ? MinElement + 5 : ChiSqVector.size() - 1);
	double OptimalSigmaParameter = FindPeak( SigmaParameters, ChiSqDiffVector, LeftSide, RightSide );
    std::cout << "Sigma Parameter Minimal Element Final Estimation " << OptimalSigmaParameter << std::endl;
    
    for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
            double SigTemp = OptimalSigmaParameter*Sig[j];
            if( ShapeOfComponent == "Gaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigTemp);
                //Value following Gaussian distribution
            }
            else if( ShapeOfComponent == "LogGaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigTemp,2)/2, SigTemp);
            }
            else
            {
                if( LifetimesFromDiscrete[j].Parameter < oPsLFLimit )
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigTemp);
                }
                else		//Option for assuming logGaussian distributions
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigTemp,2)/2, SigTemp);
                }
            }
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
	}
	return OptimalSigmaParameter;
}

void DCV::CopyParameter( std::vector< double > LifetimeGridOut, std::vector< DiscreteFitResult > LifetimesFromDiscreteOut, std::vector< DiscreteFitResult > IntensitiesFromDiscreteOut, std::vector< double > ValuesOut, mat GradMaxOut, double BackgroundOut, double SDBackgroundOut, double oPsLFLimitOut, unsigned NumberOfPreIterationsOut, unsigned NumberOfIterationsOut, double ChiDiffBetweenIterationsOut, std::vector< double > Intensities_afterIterOut )
{
    LifetimeGrid = LifetimeGridOut;
    LifetimesFromDiscrete = LifetimesFromDiscreteOut;
    IntensitiesFromDiscrete = IntensitiesFromDiscreteOut;
    Values = ValuesOut;
    GradMax = GradMaxOut;
    Background = BackgroundOut;
    SDBackground = SDBackgroundOut;
    oPsLFLimit = oPsLFLimitOut;
    NumberOfPreIterations = NumberOfPreIterationsOut;
    NumberOfIterations = NumberOfIterationsOut;
    Intensities_afterIter = Intensities_afterIterOut;
}

void DCV::CopyIteratorAndChiDiff( unsigned iteratorOut, double ChiDiffBetweenIterationsOut )
{
    iterator = iteratorOut;
    ChiDiffBetweenIterations = ChiDiffBetweenIterationsOut;
}

void DCV::ClearIteratorAndChiSQ( int Mode )
{
    if( Mode == 0 )
        iterator = 0;
    else if( Mode == 1 )
    {
        iterator = 0;
        ChiDiffBetweenIterations = 0;
    }
    else
        ChiDiffBetweenIterations = 0;        
}

unsigned DCV::GetInterator()
{
    return iterator;
}

std::vector< double > DCV::GetSigmas()
{
    return SSig;
}

std::vector< double > DCV::GetIntensities()
{
    return Intensities_afterIter;
}

double DCV::PreIterate( std::vector< double > Sig, double Step, double Lambda, double PreviousChi, std::string ShapeOfComponent )
{
	std::cout << "PreIteration " << ++iterator << std::endl;
	if( iterator%NumberOfPreIterations == 0)			//Every ten iterations condition is checked
	{
		std::cout << "PreIteration ChiSq " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )		//If this condition is fulfilled last Sigmas are saved to container SSig, and PreIteration are ended
		{
			for( unsigned i=0; i<Sig.size(); i++ )
			{
				SSig[i] = Sig[i];
			}
			return 1;
		}
		ChiDiffBetweenIterations = 0.;
	}
	double sum = 0;
	std::vector< double > sums;		//For gradient calculation

	//Initializing vector elements to 0 -> Maybe in future will be introduced some funky method to do that
	for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
	{
		sums.push_back(0);
	}
	mat Gradient_d( Sig.size(), LifetimeGrid.size() );		//Gradient over sigmas
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
            if( ShapeOfComponent == "Gaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
                //Value following Gaussian distribution
                sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]) * ( 2*pow( LifetimeGrid[i] - LifetimesFromDiscrete[j].Parameter ,2)/pow(Sig[j],3) - 1/Sig[j] );
                //Derivative of Gaussian distribution
            }
            else if( ShapeOfComponent == "LogGaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j]);
                sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j] ) * ( 2* pow( log(LifetimeGrid[i]) - log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2
                ,2)/pow(Sig[j],3) - 1/Sig[j] );
            }
            else
            {
                if( LifetimesFromDiscrete[j].Parameter < oPsLFLimit )
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
                    //Value following Gaussian distribution
                    sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]) * ( 2*pow( LifetimeGrid[i] - LifetimesFromDiscrete[j].Parameter ,2)/pow(Sig[j],3) - 1/Sig[j] );
                    //Derivative of Gaussian distribution
                    //Gradient is vector of all derivatives
                }
                else		//Option for assuming logGaussian distributions
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j]);
                    sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j] ) * ( 2* pow( log(LifetimeGrid[i]) - log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2
                    ,2)/pow(Sig[j],3) - 1/Sig[j] );
                }
            }
			Gradient_d(j,i) = sums[j];
			//Gradient is vector of all derivatives
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
		Intensities_afterIter[i] = sum;
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	vec Gradient = Gradient_d*GradMax*Diff;			//Gradient definition in the Levenberg-Marquardt method 
	mat Hessian = Gradient_d *GradMax*GradMax.t() * Gradient_d.t();		//Hessian definition in the Levenberg-Marquardt method 
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);			//Estimation of the real hessian
	vec Delta = solve(HessToInv, Gradient);					//Solving equation: Delta (Sigma) * Real Hessian = Gradient to obtain Delta (Sigma)
	vec IntensitiesDelta(Intensities.n_elem);				//Vector containing new Intensities shifted by values in Delta 
	std::vector< double > SigDelta;						//Vector for new sigmas
	for( unsigned j=0; j<Sig.size(); j++ )
	{
		SigDelta.push_back( (Sig[j] - Step*Delta[j] > 0.005 ? Sig[j] - Step*Delta[j] : 0.005) );		//Assuming that minila sigma is 10 ps - somehow reasonable value
		std::cout << LifetimesFromDiscrete[j].Parameter << " " << Sig[j] << " " << SigDelta[j] << std::endl;
	}
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )		//Creating new values of intensities with new sigmas
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
		{
            if( ShapeOfComponent == "Gaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigDelta[j]);
            }
            else if( ShapeOfComponent == "LogGaussian" )
            {
                sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigDelta[j],2)/2, SigDelta[j] );
            }
            else
            {
                if( LifetimesFromDiscrete[j].Parameter < oPsLFLimit )
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigDelta[j]);
                }
                else
                {
                    sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigDelta[j],2)/2, SigDelta[j] );
                }
            }
		}
		IntensitiesDelta(i) = sum;
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;		//New model with new sigmas
	double ChiSq2 = 0;					//ChiSquared value for sigmas corrected by Delta
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - Background - Model2(i), 2 ) / Values[i];
	std::cout << ChiSq2 << std::endl;
	if( ChiSq2 <= ChiSq )					//If new sigmas are better
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return PreIterate( SigDelta, Step, 0.1*Lambda, ChiSq2, ShapeOfComponent );	//Lambda lowered -> Gauss-Newton approach stronger -> New sigmas are promising
	}
	else
	{
		return PreIterate( Sig, Step, 10*Lambda, ChiSq, ShapeOfComponent );		//Lambda increased -> Maximal damping method stronger -> New sigma is not good enough
	}
}

double DCV::IterateExp( vec Intensities, double Step, double Lambda, double PreviousChi )		//Similar procedure as in PreIteration but fitting intesities not sigmas
{
	std::cout << "Iteration: " << ++iterator << std::endl;
	vec Model =  GradMax.t() * Intensities;
	vec Diff( Values.size() );
	double ChiSq = 0;
	double ModelDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model(i);
		ModelDiff += ( Values[i] - Model(i) )/Values.size();
	}
	std::cout << ModelDiff << " " << Background << " " << SDBackground << std::endl;
	double BGCorr = 0.;
	if( ModelDiff < Background - 2*SDBackground )
		BGCorr = Background + 2*SDBackground;
	else if( ModelDiff > Background + 2*SDBackground )
		BGCorr = Background - 2*SDBackground;
	else
		BGCorr = 2*Background - ModelDiff;
    //BGCorr = Background;          //Different approach
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - BGCorr - Model(i);
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	if( iterator%NumberOfIterations == 0 )
	{
		std::cout << "ChiSq Difference: " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.01 )
		{
			if( ChiSq > 2000 )
			{
				double Lama;
				std::cout << "You can give new Lambda value if You want" << std::endl;
				std::cin >> Lama;
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return -Lama;
			}
			else
			{
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return 1.;	    		
			}
		}
		else
		{
			for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
				Intensities_afterIter[i] = Intensities(i);
			return -0.01;
		}
	}
	vec Gradient = GradMax*Diff;
	mat Hessian = GradMax * GradMax.t();
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);
	vec Delta = solve(HessToInv, Gradient);
	vec IntensitiesDelta(Intensities.n_elem);
	for( unsigned j=0; j<Intensities.n_elem; j++ )
	{
		if( Intensities(j) + Step*Delta(j) < 0)
			IntensitiesDelta(j) = 0;
		else	
			IntensitiesDelta(j) = Intensities(j) + Step*Delta(j);
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;
	double ChiSq2 = 0;
	ModelDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model2(i);
		ModelDiff += ( Values[i] - Model2(i) )/Values.size();
	}
	std::cout << ModelDiff << " " << Background << " " << SDBackground << std::endl;
	BGCorr = 0.;
	if( ModelDiff < Background - 2*SDBackground )
		BGCorr = Background + 2*SDBackground;
	else if( ModelDiff > Background + 2*SDBackground )
		BGCorr = Background - 2*SDBackground;
	else
		BGCorr = 2*Background - ModelDiff;
    //BGCorr = Background;                      //Different approach
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - BGCorr - Model2(i), 2 ) / Values[i];
	if( ChiSq2 <= ChiSq )
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return IterateExp( IntensitiesDelta, Step, 0.1*Lambda, ChiSq2 );
	}
	else
	{
		return IterateExp( Intensities, Step, 10*Lambda, ChiSq );
	}
}

double DCV::Iterate( vec Intensities, double Step, double Lambda, double PreviousChi )		//Similar procedure as in PreIteration but fitting intesities not sigmas
{
	std::cout << "Iteration: " << ++iterator << std::endl;
	vec Model =  GradMax.t() * Intensities;
	vec Diff( Values.size() );
	double ChiSq = 0;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model(i);
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	if( iterator%NumberOfIterations == 0 )
	{
		std::cout << "ChiSq Difference: " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.01 )
		{
			if( ChiSq > 2000 )
			{
				double Lama;
				std::cout << "You can give new Lambda value if You want" << std::endl;
				std::cin >> Lama;
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return -Lama;
			}
			else
			{
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return 1.;	    		
			}
		}
		else
		{
			for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
				Intensities_afterIter[i] = Intensities(i);
			return -0.01;
		}
	}
	vec Gradient = GradMax*Diff;
	mat Hessian = GradMax * GradMax.t();
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);
	vec Delta = solve(HessToInv, Gradient);
	vec IntensitiesDelta(Intensities.n_elem);
	for( unsigned j=0; j<Intensities.n_elem; j++ )
	{
		if( Intensities(j) + Step*Delta(j) < 0)
			IntensitiesDelta(j) = 0;
		else	
			IntensitiesDelta(j) = Intensities(j) + Step*Delta(j);
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;
	double ChiSq2 = 0;
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - Background - Model2(i), 2 ) / Values[i];
	if( ChiSq2 <= ChiSq )
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return Iterate( IntensitiesDelta, Step, 0.1*Lambda, ChiSq2 );
	}
	else
	{
		return Iterate( Intensities, Step, 10*Lambda, ChiSq );
	}
}


double FindPeak( std::vector< double > Argument, std::vector< double > Values, unsigned Lowerfit, unsigned Higherfit )
{
	double XX=0, XY=0, YY=0, X=0, Y=0;
	unsigned size = Higherfit - Lowerfit;
    for( unsigned k = Lowerfit; k < Higherfit; k++ )
    {
        XX += Argument[k]*Argument[k];
        XY += Argument[k]*Values[k];
        YY += Values[k]*Values[k];
        X += Argument[k];
        Y += Values[k];
    }
    double a = (size*XY-X*Y) / (size*XX - X*X);
    double b = (Y-a*X)/size;
	double peak;
	if(a)
	    peak = -b / a;
	else 
	    peak = 0;
	return peak;
}
