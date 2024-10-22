// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>

//
//@brief Class to store the parameters of a particle to be simulated
//
struct ParticleSimParams
{
  int pid;
  int ekin;
  double eta;
  std::string param_file_key;
  std::string latex_label;
  std::string print_label;

  ParticleSimParams(int pid, int ekin, double eta, std::string param_file_key)
      : pid(pid)
      , ekin(ekin)
      , eta(eta)
      , param_file_key(std::move(param_file_key))
  {
    // Maps pid to LaTeX label
    std::map<int, std::string> pid_label = {
        {22, R"(\gamma)"}, {11, R"(e)"}, {211, R"(\pi)"}};

    // Maps ekin to LaTeX label
    std::map<int, std::string> ekin_label = {{1024, R"(1\,\text{GeV})"},
                                             {2048, R"(2\,\text{GeV})"},
                                             {4096, R"(4\,\text{GeV})"},
                                             {8192, R"(8\,\text{GeV})"},
                                             {16384, R"(16\,\text{GeV})"},
                                             {32768, R"(32\,\text{GeV})"},
                                             {65536, R"(65\,\text{GeV})"},
                                             {131072, R"(131\,\text{GeV})"},
                                             {262144, R"(262\,\text{GeV})"},
                                             {524288, R"(524\,\text{GeV})"},
                                             {1048576, R"(1\,\text{TeV})"},
                                             {2097152, R"(2\,\text{TeV})"},
                                             {4194304, R"(4\,\text{TeV})"}};

    // Generate the LaTeX label used for plotting
    std::ostringstream latex_stream;
    latex_stream << "$ E_{" << pid_label.at(pid) << "}=" << ekin_label.at(ekin)
                 << R"(\, , \eta=)" << std::fixed << std::setprecision(2) << eta
                 << "$";
    latex_label = latex_stream.str();

    // Generate print label used for Google test display
    // This is not allowed to contain any special characters
    std::ostringstream print_stream;
    print_stream << "pid" << pid << "_E" << ekin << "_eta"
                 << static_cast<int>(eta * 100);
    print_label = print_stream.str();
  }
};
