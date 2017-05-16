/*

  command line interface for parameter changing.

*/
#include "commandline.hpp"

namespace {
  const size_t CMD_LINE_ERROR = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}

namespace po = boost::program_options;

int CommandLineSettings(int argc, char** argv,double* timestep,uint* Nparticles ,uint* Niterations ,float* theta,float* softener, bool* stats, uint* freq, uint* fpsmax, bool* drawNodes)
{
  // Defining and parsing the program options from command line.
  // note: a,b,c options show how to capture a command line value directly
  // into a variable.
  po::options_description desc("Options");
  desc.add_options()
    ("help", "This program simulates solar systems using the barnes hut algorithm.")
    (",N", po::value<uint>(Nparticles), "Number of particles to simulate. default 4")
    ("deltat,d", po::value<double>(timestep), "Timestep to use for iterations (in seconds), default is 1 day (86400s)")
    (",T", po::value<uint>(Niterations), "Total number of iterations to do. Measured in days.")
    ("theta,t", po::value<float>(theta), "Value of theta to use for Barnes Hut Simulation.")
    ("softener,s", po::value<float>(softener), "Gravitational softening parameter to use (meters). default is 1e4")
    ("stats,m", po::value<bool>(stats), "Flag to monitor simulation stats like Energy and momentum")
    ("dumpfreq,f", po::value<uint>(freq), "Frequency with which to record stats. Measured in 1/steps")
    ("fps", po::value<uint>(fpsmax), "Maximum fps for SFML window playback. Also limits maximum program execution speed.")
    ("draw", po::value<bool>(drawNodes), "Draws Barnes Hut external nodes on the window if true.");
    // ADD OPTIONS FOR FILE OUTPUTS.
    // (",d", po::value<double>(&d)->required(), "variable 'd' is double, required.");

  po::variables_map varmap;
  try
  {
    // this line can throw an error, so we should catch it.
    po::store(po::parse_command_line(argc,argv,desc),varmap);

    /* --help  */
    if (varmap.count("help")) {
      std::cout << "Basic command line parameters application." << std::endl \
                << desc << std::endl;
      return SUCCESS;
    } // --help if-block
    // this throws on error, so error checking is a good idea here.
    po::notify(varmap);
  }// inner try
  catch (po::error& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    return CMD_LINE_ERROR;
  } //catch
}
