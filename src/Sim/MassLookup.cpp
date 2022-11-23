/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

*/
#include "MassLookup.h"
#include <sstream>

namespace AnasenSim {

	MassLookup* MassLookup::s_instance = new MassLookup();

	MassLookup::MassLookup()
	{
		std::ifstream massfile("etc/mass.txt");
		if(massfile.is_open())
		{
			std::string junk, element;
			double atomicMassBig, atomicMassSmall, isotopicMass;
			getline(massfile,junk);
			getline(massfile,junk);
			NucData data;
			while(massfile>>junk)
			{
				massfile>>data.Z>>data.A>>element>>atomicMassBig>>atomicMassSmall;
				data.isotopicMass = (atomicMassBig + atomicMassSmall*1e-6 - data.Z*s_electronMass)*s_u2MeV;
				data.isotopicSymbol = std::to_string(data.A) + element;
				m_dataMap[GetUUID(data.Z, data.A)] = data;
			}
		}
	}
	
	MassLookup::~MassLookup() {}
	
	//Returns nuclear mass in MeV
	double MassLookup::FindMass(uint32_t Z, uint32_t A)
	{
		auto data = m_dataMap.find(GetUUID(Z, A));
		if(data == m_dataMap.end())
			return s_invalidMass;
	
		return data->second.isotopicMass;
	}
	
	//returns element symbol
	std::string MassLookup::FindSymbol(uint32_t Z, uint32_t A)
	{
		auto data = m_dataMap.find(GetUUID(Z, A));
		if(data == m_dataMap.end())
			return s_invalidSymbol;
	
		return data->second.isotopicSymbol;
	}

}