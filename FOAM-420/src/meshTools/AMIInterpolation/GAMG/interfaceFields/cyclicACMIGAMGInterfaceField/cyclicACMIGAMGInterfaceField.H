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
    (c) 2013 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicACMIGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface for Arbitrarily Coupled Mesh Interface
    (ACMI) fields.

SourceFiles
    cyclicACMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIGAMGInterfaceField_H
#define cyclicACMIGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.H"
#include "AMIInterpolation/GAMG/interfaces/cyclicACMIGAMGInterface/cyclicACMIGAMGInterface.H"
#include "AMIInterpolation/patches/cyclicACMI/cyclicACMILduInterfaceField/cyclicACMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicACMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public cyclicACMILduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const cyclicACMIGAMGInterface& cyclicACMIInterface_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicACMIGAMGInterfaceField(const cyclicACMIGAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicACMIGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("cyclicACMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicACMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        cyclicACMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );


    //- Destructor
    virtual ~cyclicACMIGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicACMIInterface_.size();
            }


        // Interface matrix update

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        //- Cyclic interface functions

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicACMIInterface_.transform();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
