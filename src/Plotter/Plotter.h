#ifndef PLOTTER_H
#define PLOTTER_H

#include "Nucleus.h"

#include <string>
#include <memory>
#include <unordered_map>

class TObject;

namespace AnasenSim {

    struct Histogram1DParams
    {
        std::string name = "";
        std::string title = "";
        int bins = 0;
        double min = 0.0;
        double max = 0.0;
    };

    struct Histogram2DParams
    {
        std::string name = "";
        std::string title = "";
        int binsX = 0;
        double minX = 0.0;
        double maxX = 0.0;
        int binsY = 0;
        double minY = 0.0;
        double maxY = 0.0;
    };

    struct GraphParams
    {
        std::string name = "";
        std::string title = "";
    };

    class Plotter
    {
    public:
        Plotter(const std::string& inputname, const std::string& outputname);
        ~Plotter();

        void Run();

    private:
        void PlotNucleus(const Nucleus& nucleus);

        void FillHistogram1D(const Histogram1DParams& params, double value);
        void FillHistogram2D(const Histogram2DParams& params, double valueX, double valueY);
        void FillGraph(const GraphParams& params, double valueX, double valueY);

        std::string m_inputName;
        std::string m_outputName;

        std::unordered_map<std::string, std::shared_ptr<TObject>> m_map;
    };
}

#endif