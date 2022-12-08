#include "Application.h"
#include "TFile.h"
#include "TTree.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"

#include <fstream>
#include <iostream>

namespace AnasenSim {

    Application::Application(const std::filesystem::path& config) :
        m_isInit(false)
    {
		if(!EnforceDictionaryLinked())
		{
			std::cerr << "Dictionary error!" << std::endl;
		}
        InitConfig(config);
    }

    Application::~Application()
    {
    }

    void Application::InitConfig(const std::filesystem::path& config)
    {
		if(std::filesystem::exists(config))
        {
            std::cerr << "Unable to open configuration file " << config << " because it does not exist" << std::endl;
            return;
        }
        std::ifstream configFile(config);

		std::cout << "Parsing configuration file: " << config << std::endl;

		SystemParameters params;
        std::string junk;
		configFile>>junk>>m_outputName;
		configFile>>junk>>m_nThreads;
		m_nThreads = m_nThreads > 0 ? m_nThreads : 1; //Enforce that we must always have one processing thread
		configFile>>junk>>m_nSamples;
		m_processingChunks.resize(m_nThreads);
		m_threadPool = std::make_unique<ThreadPool<Chunk*>>(m_nThreads);
		m_fileWriter.Open(m_outputName, "SimTree");
		
		//ensure that all samples are processed even when m_nSamples isn't evenly divisible by m_nThreads
		uint64_t quotient = m_nSamples / m_nThreads;
    	uint64_t remainder = m_nSamples % m_nThreads;
    	m_processingChunks[0].samples = quotient + remainder;
		for(uint64_t i=1; i<m_processingChunks.size(); i++)
			m_processingChunks[i].samples = quotient;

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
                configFile >> junk >> params.initialBeamEnergy
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

		for(auto& chunk : m_processingChunks)
		{
			chunk.system = CreateSystem(params);
			chunk.array = new AnasenArray(params.target);
			if(chunk.system == nullptr || !chunk.system->IsValid())
			{
				std::cerr<<"Failure to parse reaction system... configuration not loaded"<<std::endl;
				return;
			}
		}

		std::getline(configFile, junk);
		std::getline(configFile, junk);

		std::cout << "Output file: " << m_outputName << std::endl;
		std::cout << "Reaction equation: " << m_processingChunks[0].system->GetSystemEquation() << std::endl;
		std::cout << "Number of samples: " << m_nSamples << std::endl;
		std::cout << "Number of threads: " << m_nThreads << std::endl;

		std::cout << "Configuration loaded successfully" << std::endl;

        m_isInit = true;
    }

	void Application::Run()
	{
		if(!m_fileWriter.IsOpen())
		{
			std::cerr << "Could not open output file " << m_outputName << " at Application::Run() " << std::endl;
			return;
		}

		if(m_processingChunks.size() != m_nThreads)
		{
			std::cerr << "System list not equal to number of threads" << std::endl;
			return;
		}

		//Give our thread pool some tasks
		for(std::size_t i=0; i<m_processingChunks.size(); i++)
		{
			//bind a lambda to the job, taking in a ReactionSystem, and then provide a reaction system as the tuple arguments.
			m_threadPool->PushJob({[this](Chunk* processChunk) 
				{
					if(processChunk->system == nullptr || !processChunk->system->IsValid())
						return;
					
					std::vector<Nucleus>* eventHandle = processChunk->system->GetNuclei();

					for(uint64_t i=0; i<processChunk->samples; i++)
					{
						processChunk->system->RunSystem();
						for(Nucleus& nucleus : *eventHandle)
						{
							processChunk->array->IsDetected(nucleus);
						}
						m_fileWriter.PushData(*eventHandle);
						processChunk->system->ResetNucleiDetected();
					}
				}, 
				{&m_processingChunks[i]}
			});
		}

		uint64_t count = 0;
		double percent = 0.01;
		uint64_t flushVal = m_nSamples*percent;
		uint64_t flushCount = 0;
		std::cout << "Starting simulation..." << std::endl;
		while(true)
		{
			if(count == flushVal)
			{
				count = 0;
				++flushCount;
				std::cout<<"\rPercent of data written to disk: "<<percent*flushCount*100<<"%"<<std::flush;
			}

			if(m_threadPool->IsFinished() && m_fileWriter.GetQueueSize() == 0)
				break;
			else if(m_fileWriter.Write())
				++count;
		}

		std::cout<<std::endl;
		std::cout<<"Simulation complete."<<std::endl;
	}
}