#include "Params.h"
#include "ParamsGeneration.h"

int main(int argc, const char** argv) {
    std::string filename = "SimulationParams.cs";
    ParamsGeneration::makeParamsGeneratorClass(filename);
    return EXIT_SUCCESS;
}