#ifndef SIM_APP_H
#define SIM_APP_H

#include "ReactionSystem.h"
#include <string>
#include <vector>

namespace AnasenSim {

    class Application
    {
    public:
        Application(const std::string& config);
        ~Application();

        void Run();

    private:
        void InitConfig(const std::string& config);

        bool m_isInit;

        std::string m_outputName = "";
        uint64_t m_nSamples = 0;
        Target m_target;
        std::unique_ptr<ReactionSystem> m_system;
    };
}

#endif