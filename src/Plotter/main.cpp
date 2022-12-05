#include "Plotter.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        std::cerr << "AnasenSim Plotter requires an input simulation file and an output file!" << std::endl;
        return 1;
    }

    AnasenSim::Plotter plotter(argv[1], argv[2]);

    plotter.Run();
}