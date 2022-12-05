#include "Plotter.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <sstream>

namespace AnasenSim {
    
    Plotter::Plotter(const std::string& inputname, const std::string& outputname) :
        m_inputName(inputname), m_outputName(outputname)
    {
        TH1::AddDirectory(kFALSE); //Force ROOT to let us own the histograms

        if(!EnforceDictionaryLinked())
        {
            std::cout << "Dictionary error" << std::endl;
        }
    }

    Plotter::~Plotter() {}

    void Plotter::FillHistogram1D(const Histogram1DParams& params, double value)
    {
        auto iter = m_map.find(params.name);
        if(iter == m_map.end())
        {
            std::shared_ptr<TH1F> histo = std::make_shared<TH1F>(params.name.c_str(), params.title.c_str(), params.bins, params.min, params.max);
            m_map[params.name] = std::static_pointer_cast<TObject>(histo);
            histo->Fill(value);
        }
        else
        {
            auto histo = std::static_pointer_cast<TH1>(iter->second);
            if(histo)
                histo->Fill(value);
        }
    }

    void Plotter::FillHistogram2D(const Histogram2DParams& params, double valueX, double valueY)
    {
        auto iter = m_map.find(params.name);
        if(iter == m_map.end())
        {
            std::shared_ptr<TH2F> histo = std::make_shared<TH2F>(params.name.c_str(), params.title.c_str(), params.binsX, params.minX, params.maxX, params.binsY, params.minY, params.maxY);
            m_map[params.name] = std::static_pointer_cast<TObject>(histo);
            histo->Fill(valueX, valueY);
        }
        else
        {
            auto histo = std::static_pointer_cast<TH2>(iter->second);
            if(histo)
                histo->Fill(valueX, valueY);
        }
    }

    void Plotter::FillGraph(const GraphParams& params, double valueX, double valueY)
    {
        auto iter = m_map.find(params.name);
        if(iter == m_map.end())
        {
            std::shared_ptr<TGraph> graph = std::make_shared<TGraph>(1, &valueX, &valueY);
            graph->SetName(params.name.c_str());
            graph->SetTitle(params.title.c_str());
            m_map[params.name] = std::static_pointer_cast<TObject>(graph);
        }
        else
        {
            auto graph = std::static_pointer_cast<TGraph>(iter->second);
            if(graph)
                graph->AddPoint(valueX, valueY);
        }
    }

    void Plotter::Run()
    {
        TFile* inputFile = TFile::Open(m_inputName.c_str(), "READ");
        if(!inputFile || !inputFile->IsOpen())
        {
            std::cerr << "Unable to open input file " << m_inputName << std::endl;
            return;
        }

        TTree* simTree = (TTree*) inputFile->Get("SimTree");
        if(!simTree)
        {
            std::cerr << "Unable to retrieve SimTree from input file " << m_inputName << std::endl;
            inputFile->Close();
            delete inputFile;
            return;
        }

        std::vector<Nucleus>* eventHandle = new std::vector<Nucleus>();
        simTree->SetBranchAddress("event", &eventHandle);

        TFile* outputFile = TFile::Open(m_outputName.c_str(), "RECREATE");
        if(!outputFile || !outputFile->IsOpen())
        {
            std::cerr << "Unable to open output file " << m_outputName << std::endl;
            inputFile->Close();
            delete inputFile;
            delete eventHandle;
            return;
        }

        std::cout << "Generating plots from simulation data in " << m_inputName << " and writing to " << m_outputName << std::endl;

        uint64_t nentries = simTree->GetEntries();
        uint64_t count = 0;
        double flushFrac = 0.01;
        uint64_t flushCount = 0;
        uint64_t flushVal = flushFrac * nentries;

        std::cout << "Plotting..." << std::endl;
        for(uint64_t i=0; i<nentries; i++)
        {
            ++count;
            if(count == flushVal)
            {
                ++flushCount;
                count = 0;
                std::cout << "\rPercent of data processed: " << flushCount * flushFrac * 100 << "%" << std::flush;
            }

            simTree->GetEntry(i);

            for(const Nucleus& nucleus : *eventHandle)
            {
                PlotNucleus(nucleus);
            }
        }

        inputFile->Close();
        outputFile->cd();
        for(auto& iter : m_map)
            iter.second->Write(iter.second->GetName(), TObject::kOverwrite);
        outputFile->Close();

        delete inputFile;
        delete outputFile;
        delete eventHandle;
        std::cout << std::endl << "Complete." << std::endl;
    }

    void Plotter::PlotNucleus(const Nucleus& nucleus)
    {
        std::stringstream nucleusStream;
        nucleusStream << nucleus.isotopicSymbol << "_" << ReactionRoleToString(nucleus.role);
        FillGraph({nucleusStream.str() + "_KE_theta", nucleusStream.str() + ";#theta_{lab};KE (MeV)"}, nucleus.vec4.Theta() * s_rad2deg, nucleus.GetKE());
        FillGraph({nucleusStream.str() + "_KE_phi", nucleusStream.str() + ";#phi_{lab};KE (MeV)"}, FullPhi(nucleus.vec4.Phi()) * s_rad2deg, nucleus.GetKE());
        FillGraph({nucleusStream.str() + "_rxnX_rxnY", nucleusStream.str() + ";rxnX (m);rxnY (m)"}, nucleus.rxnPoint.X(), nucleus.rxnPoint.Y());
        FillHistogram1D({nucleusStream.str() + "_rxnZ", nucleusStream.str() + ";rxnZ (m);", 554, 0.0, 0.554}, nucleus.rxnPoint.Z());
        if(nucleus.isDetected)
        {
            FillGraph({nucleusStream.str() + "_KE_theta_det", nucleusStream.str() + ";#theta_{lab};KE (MeV)"}, nucleus.vec4.Theta() * s_rad2deg, nucleus.siliconDetKE);
            FillGraph({nucleusStream.str() + "_KE_phi_det", nucleusStream.str() + ";#phi_{lab};KE (MeV)"}, FullPhi(nucleus.vec4.Phi()) * s_rad2deg, nucleus.siliconDetKE);
            FillGraph({nucleusStream.str() + "_rxnX_rxnY_det", nucleusStream.str() + ";rxnX (m);rxnY (m)"}, nucleus.rxnPoint.X(), nucleus.rxnPoint.Y());
            FillHistogram2D({"EdE_pcE_siKE", "EdE;Si KE(MeV);PC E(MeV)", 200, 0.0, 35.0, 200, 0.0, 1.5}, nucleus.siliconDetKE, nucleus.pcDetE);
            FillHistogram1D({nucleusStream.str() + "_rxnZ_det", nucleusStream.str() + ";rxnZ (m);", 554, 0.0, 0.554}, nucleus.rxnPoint.Z());
        }
    }
}