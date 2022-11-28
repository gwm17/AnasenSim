#include "Plotter.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

namespace AnasenSim {

    Plotter::Plotter(const std::string& inputname, const std::string& outputname) :
        m_inputName(inputname), m_outputName(outputname)
    {
        TH1::AddDirectory(kFALSE); //Force ROOT to let us own the histograms
    }

    void Plotter::FillHistogram1D(const Histogram1DParams& params, double value)
    {
        auto iter = m_map.find(params.name);
        if(iter == m_map.end())
        {
            std::shared_ptr<TH1F> histo = std::make_shared<TH1F>(params.name, params.title, params.bins, params.min, params.max);
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
            std::shared_ptr<TH2F> histo = std::make_shared<TH2F>(params.name, params.title, params.binsX, params.minX, params.maxX, params.binsY, params.minY, params.maxY);
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

        uint64_t nentries = simTree->GetEntries();
        uint64_t count = 0;
        double flushFrac = 0.01;
        uint64_t flushCount = 0;
        uint64_t flushVal = flushFrac * nentries;

        for(uint64_t i=0; i<nentries; i++)
        {
            ++count;
            if(count == flushVal)
            {
                ++flushCount;
                count = 0;
                std::cout << "Percent of data processed: " << flushCount * flushFrac * 100 << "%" << std::flush;
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
    }

    void Plotter::PlotNucleus(const Nucleus& nucleus)
    {
        
    }
}