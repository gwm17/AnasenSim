#include "Application.h"
#include "TFile.h"
#include "TTree.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"

#include <fstream>
#include <iostream>

namespace AnasenSim {

    Application::Application(const std::string& config) :
        m_isInit(false), m_system(nullptr), m_array(nullptr)
    {
		if(!EnforceDictionaryLinked())
		{
			std::cerr << "Dictionary error!" << std::endl;
		}
        InitConfig(config);
    }

    Application::~Application()
    {
		delete m_system;
		delete m_array;
    }

    void Application::InitConfig(const std::string& config)
    {
        std::ifstream configFile(config);
        if(!configFile.is_open())
        {
            std::cerr << "Unable to open configuration file " << config << std::endl;
            return;
        }

		std::cout << "Parsing configuration file: " << config << std::endl;

		SystemParameters params;
        std::string junk;
		configFile>>junk>>m_outputName;

		double density;
		std::vector<uint32_t> avec, zvec;
		std::vector<int> svec;
		uint32_t z, a;
		int s;
		
		configFile>>junk>>junk>>density;
		avec.clear(); zvec.clear(); svec.clear();
		while(configFile>>junk) 
		{
			if(junk == "begin_elements")
			{
				configFile>>junk>>junk>>junk;
				continue;
			} 
			else if (junk == "end_elements")
				break;
			configFile>>z>>a>>s;
			zvec.push_back(z); avec.push_back(a); svec.push_back(s);
		}
		configFile>>junk;
		params.target = Target(zvec, avec, svec, density);

		while(configFile>>junk)
		{
			if(junk == "begin_chain")
            {
                configFile >> junk >> m_nSamples >> junk >> params.initialBeamEnergy
                           >> junk >> params.rxnBeamEnergy;
				continue;
            }
			else if (junk == "end_chain")
				break;
			else if(junk == "begin_step")
			{
				StepParameters currentParams;
				configFile >> junk >> junk;
				currentParams.rxnType = StringToRxnType(junk);
				if(currentParams.rxnType == RxnType::Reaction)
				{
					configFile >> junk;
					for(int i=0; i<3; i++)
					{
						configFile >> z >> a;
						currentParams.Z.push_back(z);
						currentParams.A.push_back(a);
					}
					configFile >> junk;
					configFile >> junk >> currentParams.meanResidualEx;
					configFile >> junk >> currentParams.sigmaResidualEx;
					params.stepParams.push_back(currentParams);
				}
				else if(currentParams.rxnType == RxnType::Decay)
				{
					configFile >> junk;
					for(int i=0; i<2; i++)
					{
						configFile >> z >> a;
						currentParams.Z.push_back(z);
						currentParams.A.push_back(a);
					}
					configFile >> junk;
					configFile >> junk >> currentParams.meanResidualEx;
					configFile >> junk >> currentParams.sigmaResidualEx;
					params.stepParams.push_back(currentParams);
				}
				else
				{
					std::cerr << "Invalid reaction information at SimApp::InitConfig!" << std::endl;
					return;
				}
			}
		}

		m_system = CreateSystem(params);
		m_array = new AnasenArray(params.target);
		if(m_system == nullptr || !m_system->IsValid())
		{
			std::cerr<<"Failure to parse reaction system... configuration not loaded."<<std::endl;
			return;
		}

		std::getline(configFile, junk);
		std::getline(configFile, junk);

		std::cout<<"Reaction equation: "<<m_system->GetSystemEquation()<<std::endl;
		std::cout<<"Number of samples: "<<m_nSamples<<std::endl;

		std::cout << "Configuration loaded successfully" << std::endl;

        m_isInit = true;
    }

    void Application::Run()
    {
        if(!m_isInit)
        {
            std::cerr << "Application not initialized at Application::Run()!" << std::endl;
            return;
        }

        TFile* outputFile = TFile::Open(m_outputName.c_str(), "RECREATE");
        if(!outputFile || !outputFile->IsOpen())
        {
            std::cerr << "Could not open output file " << m_outputName << " at Application::Run() " << std::endl;
            return;
        }

        TTree* outtree = new TTree("SimTree", "SimTree");
        outtree->Branch("event", m_system->GetNuclei());

        double flushPercent = 0.01;
        uint64_t flushVal = flushPercent * m_nSamples;
        uint64_t count = 0, flushCount = 0;

		std::cout << "Starting simulation..." << std::endl;

		std::vector<Nucleus>* eventHandle = m_system->GetNuclei();

        for(uint64_t i=0; i<m_nSamples; i++)
        {
            count++;
            if(count == flushVal)
            {
                count = 0;
                flushCount++;
                std::cout << "\rPercent of data simulated: " << flushCount * flushPercent * 100 << "%" << std::flush;
            }

            m_system->RunSystem();
			for(Nucleus& nucleus : *eventHandle)
			{
				m_array->IsDetected(nucleus);
			}
            outtree->Fill();
			m_system->ResetNucleiDetected();
        }

        outputFile->cd();
        outtree->Write(outtree->GetName(), TObject::kOverwrite);
        outputFile->Close();
        delete outputFile;

		std::cout << std::endl << "Simulation complete" << std::endl;
    }
}