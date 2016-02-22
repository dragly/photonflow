#include "units.h"

#include <ostream>

using namespace std;

namespace photonflow {

ostream& operator<< (ostream& out, Length l)
{
    out << l.value() << " lengths";
    return out;
}

}
