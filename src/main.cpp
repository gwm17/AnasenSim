#include "Sim/Application.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cerr << "Err! AnasenSim needs a configuration file to run!" << std::endl;
        return 1;
    }

    AnasenSim::Application* myApp = new AnasenSim::Application(argv[1]);

    if(!myApp->IsInit())
    {
        std::cerr << "Configuration file was not loaded correctly!" << std::endl;
        return 1;
    }

    myApp->Run();

    delete myApp;
}