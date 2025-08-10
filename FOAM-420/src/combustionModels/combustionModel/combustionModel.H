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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::combustionModel

Group
    grpCombustionModels

Description
    Base class for combustion models

SourceFiles
    combustionModel.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModel_H
#define combustionModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class combustionModel Declaration
\*---------------------------------------------------------------------------*/

class combustionModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Disallow copy construct
        combustionModel(const combustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const combustionModel&);


protected:

    // Protected data

        //- Reference to the turbulence model
        mutable compressibleTurbulenceModel* turbulencePtr_;

        //- Reference to solution database
        const objectRegistry& obr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Active
        Switch active_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model type
        const word modelType_;

        //- Phase name
        const word phaseName_;

        const fvMesh& getMesh(const objectRegistry&) const;
public:

    //- Runtime type information
    TypeName("combustionModel");

    //- Default combustionProperties dictionary name
    static const word combustionPropertiesName;


    // Constructors

        //- Construct from components
        combustionModel
        (
            const word& modelType,
            const objectRegistry& mesh,
            const word& combustionProperties=combustionPropertiesName,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~combustionModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Set turbulence
        inline void setTurbulence(compressibleTurbulenceModel& turbModel);

        //- Return access to turbulence
        inline const compressibleTurbulenceModel& turbulence() const;

        //- Return const access to rho
        inline const volScalarField& rho() const;

        //- Return const access to phi
        inline tmp<surfaceScalarField> phi() const;

        //- Is combustion active?
        inline const Switch& active() const;

        //- Return const dictionary of the model
        inline const dictionary& coeffs() const;

        //- Correct combustion rate
        virtual void correct() = 0;

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const = 0;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "combustionModel/combustionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
