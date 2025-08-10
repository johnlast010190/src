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
    (c) 2010-2022 Esi Ltd.

Class
    Foam::mappedGAMGInterfaceField

Description
    GAMG agglomerated region coupled interface field.

SourceFiles
    mappedGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedGAMGInterfaceField_H
#define mappedGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.H"
#include "mappedPatches/GAMG/mappedGAMGInterface.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mappedGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class mappedGAMGInterfaceField
:
    public GAMGInterfaceField
{
    // Private data

        //- Local reference cast into the region coupled interface
        const mappedGAMGInterface& regionCoupledGAMGInterface_;

        //- Agglomerated explicit correction (for radiation)
        scalarField explicitNbrCorr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mappedGAMGInterfaceField
        (
            const mappedGAMGInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const mappedGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("mappedPatch");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        mappedGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        mappedGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );


    //- Destructor
    virtual ~mappedGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return regionCoupledGAMGInterface_.size();
            }

            virtual tmp<scalarField> explicitNbrCorr() const
            {
                return tmp<scalarField>(new scalarField(explicitNbrCorr_));
            }


    //- Interface matrix update
            virtual void updateInterfaceMatrix
            (
                scalarField&,
                const bool add,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
