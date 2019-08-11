#include "PALS_avalanche_libTools.h"

LifetimeComponent::LifetimeComponent()			//Initializing default values for Lifetime Component 
{
	Type = "";
	Lifetime = 0.1;
	Intensity = 0.0001;
}

LifetimeComponent::LifetimeComponent( float lifetime, float intensity, std::string type )
{							//Initializing values for Lifetime Component given by user 
	Lifetime = lifetime;
	Intensity = intensity;
	Type = type;
}

DiscreteFitResult::DiscreteFitResult()			//Initializing default values for Fitted parameter 
{
	Parameter = 1;
	Uncertainity = 0.01;
	Type = "nf";	
}

DiscreteFitResult::DiscreteFitResult( double Par, double Unc )
{							//Initializing values for Fitted parameter 
	Parameter = Par;
	Uncertainity = Unc;
	Type = "nf";
}