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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicAMIGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    cyclicAMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIGAMGInterfaceField_H
#define cyclicAMIGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.H"
#include "AMIInterpolation/GAMG/interfaces/cyclicAMIGAMGInterface/cyclicAMIGAMGInterface.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMILduInterfaceField/cyclicAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicAMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public cyclicAMILduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const cyclicAMIGAMGInterface& cyclicAMIInterface_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicAMIGAMGInterfaceField(const cyclicAMIGAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicAMIGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        cyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );


    //- Destructor
    virtual ~cyclicAMIGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicAMIInterface_.size();
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
                return cyclicAMIInterface_.transform();
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
