#ifndef FILE_WRITER_H
#define FILE_WRITER_H


#include "Nucleus.h"

#include "TFile.h"
#include "TTree.h"

#include <queue>
#include <mutex>
#include <atomic>

namespace AnasenSim {

    class FileWriter
    {
    public:
        FileWriter();
        FileWriter(const std::string& filename, const std::string& treename);
        ~FileWriter();

        bool IsOpen() const { return m_file == nullptr ? false : m_file->IsOpen(); }
        bool IsTree() const { return m_tree == nullptr ? false : true; }

        std::size_t GetQueueSize() const { return m_queueSize; } //Implicitly thread-safe

        void PushData(const std::vector<Nucleus>& data); //Thread-safe
        bool Write(); //Not completely thread-safe, should only be used by main application loop (acess to queue is safe, but all else not)

        void Open(const std::string& filename, const std::string& treename); //Not thread safe!
        void Close(); //Not thread safe!

    private:
        TFile* m_file;
        TTree* m_tree;

        std::vector<Nucleus> m_dataHandle;

        std::mutex m_queueMutex;
        std::atomic<std::size_t> m_queueSize;
        std::queue<std::vector<Nucleus>> m_queue;
    };
}

#endif