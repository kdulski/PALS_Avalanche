


class LifetimeComponent
{
	public:
		float Lifetime;
		float Intensity;
		std::string Type;

		LifetimeComponent( float lifetime, float intensity, std::string type );	
		LifetimeComponent();
};

class DiscreteFitResult
{
	public:
		double Parameter;
		double Uncertainity;
		std::string Type;
		
		DiscreteFitResult();
		DiscreteFitResult( double Par, double Unc );
};