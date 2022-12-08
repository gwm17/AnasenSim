#include "FileWriter.h"

namespace AnasenSim {

    FileWriter::FileWriter() :
        m_file(nullptr), m_tree(nullptr), m_queueSize(0)
    {
    }

    FileWriter::FileWriter(const std::string& filename, const std::string& treename) :
        m_file(nullptr), m_tree(nullptr)
    {
        m_file = TFile::Open(filename.c_str(), "RECREATE");
        if(m_file != nullptr && m_file->IsOpen())
        {
            m_tree = new TTree(treename.c_str(), treename.c_str());
            m_tree->Branch("event", &m_dataHandle);
        }
    }

    FileWriter::~FileWriter()
    {
        if(m_file != nullptr || m_tree != nullptr)
            Close();
    }

    void FileWriter::Open(const std::string& filename, const std::string& treename)
    {
        if(m_file != nullptr || m_tree != nullptr)
            Close();

        m_file = TFile::Open(filename.c_str(), "RECREATE");
        if(m_file != nullptr && m_file->IsOpen())
        {
            m_tree = new TTree(treename.c_str(), treename.c_str());
            m_tree->Branch("event", &m_dataHandle);
        }
    }

    void FileWriter::Close()
    {
        if(m_file->IsOpen())
        {
            if(m_tree != nullptr)
                m_tree->Write(m_tree->GetName(), TObject::kOverwrite);

            m_file->Close();
            delete m_file;
            m_file = nullptr;
        }
    }

    void FileWriter::PushData(const std::vector<Nucleus>& data)
    {
        std::scoped_lock<std::mutex> guard(m_queueMutex);
        m_queue.push(data);
        ++m_queueSize;
    }

    bool FileWriter::Write()
    {
        if(m_queueSize == 0)
            return false;

        //Aquire lock for as short a time as possible
        {
            std::scoped_lock<std::mutex> guard(m_queueMutex);
            m_dataHandle = m_queue.front();
            m_queue.pop();
        }
        --m_queueSize;

        m_tree->Fill();
        return true;
    }
}