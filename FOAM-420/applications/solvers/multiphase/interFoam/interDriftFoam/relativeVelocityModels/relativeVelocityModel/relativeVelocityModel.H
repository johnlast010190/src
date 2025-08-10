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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::relativeVelocityModel

Description

SourceFiles
    relativeVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVelocityModel_H
#define relativeVelocityModel_H

#include "cfdTools/general/include/fvCFD.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relativeVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class relativeVelocityModel
{
    // Private Member Functions

        //- Return the list of patchFieldTypes for Udm derived from U
        wordList UdmPatchFieldTypes() const;

        //- Disallow default bitwise copy construct
        relativeVelocityModel(const relativeVelocityModel&);

        //- Disallow default bitwise assignment
        void operator=(const relativeVelocityModel&);


protected:

    // Protected data

        //- Dispersed phase fraction
        const volScalarField& alpha_;

        //- Dispersed diffusion velocity
        volVectorField Udm_;


public:

    //- Runtime type information
    TypeName("relativeVelocityModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relativeVelocityModel,
        dictionary,
        (const word& phaseName,
         const dictionary& dict,
         const volScalarField& alpha),
        (phaseName, dict, alpha)
    );


    // Constructors

        //- Construct from components
        relativeVelocityModel
        (
            const word& phaseName,
            const dictionary& dict,
            const volScalarField& alpha
        );


    // Selector
    static autoPtr<relativeVelocityModel> New
    (
        const word& phaseName,
        const dictionary& dict,
        const volScalarField& alpha
    );


    //- Destructor
    virtual ~relativeVelocityModel();


    // Member Functions

        //- Return the diffusion velocity of the dispersed phase
        const volVectorField& Udm() const
        {
            return Udm_;
        }

        //- Update the diffusion velocity
        virtual void correct() = 0;

        //- Return the alpha value corresponding to the max
        //  drift velocity (needed for consistent discretisation)
        virtual scalar alpha23() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
