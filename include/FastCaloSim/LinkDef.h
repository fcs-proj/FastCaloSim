#include "FastCaloSim/MLogging.h"

#ifdef __CLING__
// Standard preamble: turn off creation of dictionaries for "everything":
// we then turn it on only for the types we are interested in.
#  pragma link off all globals;
#  pragma link off all classes;
#  pragma link off all functions;
// Turn on creation of dictionaries for nested classes
#  pragma link C++ nestedclasses;

// Custom classes
#  pragma link C++ class ISF_FCS::MLogging + ;

#endif  // __CLING__