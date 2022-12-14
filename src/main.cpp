#include "Sim/Application.h"
#include "Detectors/AnasenArray.h"
#include "Utils/Timer.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cerr << "Err! AnasenSim needs a configuration file to run!" << std::endl;
        return 1;
    }

    //AnasenSim::AnasenArray array({});
    //array.DrawDetectorSystem("etc/array.txt");

    
    AnasenSim::Application* myApp = new AnasenSim::Application(argv[1]);

    if(!myApp->IsInit())
    {
        std::cerr << "Configuration file was not loaded correctly!" << std::endl;
        delete myApp;
        return 1;
    }

    AnasenSim::Timer watch;
    watch.Start();
    myApp->Run();
    watch.Stop();

    std::cout << "Elapsed time: " << watch.GetElapsedMilliseconds() << " ms" << std::endl;

    delete myApp;
    
}