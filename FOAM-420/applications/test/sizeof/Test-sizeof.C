/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2011 OpenFOAM Foundation

Description
    Test the sizeof various classes.

\*---------------------------------------------------------------------------*/

#include "primitives/bools/bool/bool.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/strings/string/string.H"
#include "db/dictionary/dictionary.H"
#include "primitives/nil/nil.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

namespace Foam
{
   class hasBoolClass
   {
   public:
      bool b_;

      hasBoolClass(const bool val=false)
      :
          b_(false)
      {}
   };

}


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    cout<<"sizeof\n------\n";
    {
        nil x;
        cout<<"nil:" << sizeof(x) << nl;
    }
    {
        bool x(0);
        cout<<"bool:" << sizeof(x) << nl;
    }
    {
        hasBoolClass x(true);
        cout<<"hasBoolClass:" << sizeof(x) << nl;
    }

    {
        Switch x("n");
        cout<<"Switch:" << sizeof(x) << nl;
        cout<<"Switch::switchType=" << sizeof(Switch::switchType) << nl;
    }

    {
        scalar x(0);
        cout<<"scalar:" << sizeof(x) << nl;
    }

    {
        label x(0);
        cout<<"label:" << sizeof(x) << nl;
    }

    {
        cout<<"int:" << sizeof(int) << nl;
        cout<<"long:" << sizeof(long) << nl;
        cout<<"float:" << sizeof(float) << nl;
        cout<<"double:" << sizeof(double) << nl;
    }


    Info<< "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
