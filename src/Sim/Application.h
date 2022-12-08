#ifndef SIM_APP_H
#define SIM_APP_H

#include "ReactionSystem.h"
#include "Detectors/AnasenArray.h"
#include "FileWriter.h"
#include "ThreadPool.h"

#include <string>
#include <vector>
#include <memory>

namespace AnasenSim {

    class Application
    {
    public:

        struct Chunk
        {
            ReactionSystem* system = nullptr;
            AnasenArray* array = nullptr;
            uint64_t samples = 0;
        };

        Application(const std::string& config);
        ~Application();

        void Run();
        void RunSingleThread();

        bool IsInit()  const { return m_isInit; }
    private:
        void InitConfig(const std::string& config);

        bool m_isInit;

        std::string m_outputName = "";
        uint64_t m_nSamples = 0;
        uint32_t m_nThreads = 0;

        std::vector<Chunk> m_processingChunks;
        FileWriter m_fileWriter;
        std::unique_ptr<ThreadPool<Chunk*>> m_threadPool;

        Target m_target;
        ReactionSystem* m_system;
        AnasenArray* m_array;

    };
}

#endif