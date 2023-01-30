/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

Converted to true singleton to simplify usage -- Aug. 2021 GWM
*/
#ifndef MASS_LOOKUP_H
#define MASS_LOOKUP_H

#include <fstream>
#include <string>
#include <unordered_map>

namespace AnasenSim {

	class MassLookup
	{
	public:

		struct NucData
		{
			std::string isotopicSymbol;
			double isotopicMass;
			uint32_t Z;
			uint32_t A;
		};

		~MassLookup();
		double FindMass(uint32_t Z, uint32_t A);
		double FindMassU(uint32_t Z, uint32_t A) { return FindMass(Z, A)/s_u2MeV; }
		std::string FindSymbol(uint32_t Z, uint32_t A);
	
		static bool IsInvalidSymbol(const std::string& symbol) { return symbol == s_invalidSymbol; }
		static bool IsInvalidMass(double mass) { return mass == s_invalidMass; }
		static MassLookup& GetInstance() { return *s_instance; }
	
	private:
		MassLookup();

		static MassLookup* s_instance;

		std::unordered_map<uint32_t, NucData> m_dataMap;
		//constants
		static constexpr double s_u2MeV = 931.4940954;
		static constexpr double s_electronMass = 0.000548579909;
		static constexpr double s_invalidMass = -1.0; //MeV
		static constexpr char s_invalidSymbol[] = "Invalid";

	};

}

#endif
