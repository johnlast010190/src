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
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::regionModels::pyrolysisModels::pyrolysisModelCollection

Description
    A centralized pyrolysis collection.

    Container class for a set of pyrolysis with functions implemented
    to loop over the functions for each type.

SourceFiles
    pyrolysisModelCollection.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisModelCollection_H
#define pyrolysisModelCollection_H

#include "containers/Lists/PtrList/PtrList.H"
#include "pyrolysisModel/pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class fvMesh;

namespace regionModels
{
namespace pyrolysisModels
{

/*---------------------------------------------------------------------------*\
                  Class pyrolysisModelCollection Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisModelCollection
:
    public PtrList<pyrolysisModel>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pyrolysisModelCollection(const pyrolysisModelCollection&);

        //- Disallow default bitwise assignment
        void operator=(const pyrolysisModelCollection&);


public:

    // Runtime type information
    TypeName("pyrolysisModelCollection");


    // Constructors

        //- Construct from mesh
        pyrolysisModelCollection(const fvMesh&);


    //- Destructor
    virtual ~pyrolysisModelCollection();


    // Member Functions

        //- Pre-evolve regions
        virtual void preEvolveRegion();

        //- Evolve the pyrolysis equation regions
        virtual void evolveRegion();

        //- Evolve regions
        virtual void evolve();

        //- Provide some feedback from pyrolysis regions
        virtual void info();

        //- Return max diffusivity allowed in the solid
        virtual scalar maxDiff() const;

        //- Mean diffusion number of the solid regions
        virtual scalar solidRegionDiffNo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam


#endif

// ************************************************************************* //
