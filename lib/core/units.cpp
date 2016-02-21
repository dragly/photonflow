#include "units.h"

#include <ostream>

using namespace std;
using namespace boost::units::photonflow;

ostream& operator<< (ostream& out, Length l)
{
    out << l.value() << " lengths";
    return out;
}
