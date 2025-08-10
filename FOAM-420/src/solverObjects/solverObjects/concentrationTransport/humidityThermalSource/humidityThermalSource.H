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
    (c) 2015 OpenFOAM Foundation
    (c) 2021 Esi Ltd

Class
    Foam::fv::humidityThermalSource

Group
    grpFvOptionsSources

Description
    Calculates and applies the vapour flux contribution to the energy
    equation.

    \heading Source usage
    Example usage:
    \verbatim
    humidityThermalSourceCoeffs
    {
        fieldNames      (T);    //h              // Name of energy field
    }
    \endverbatim

SourceFiles
    humidityThermalSource.C

\*---------------------------------------------------------------------------*/
#ifndef humidityThermalSource_H
#define humidityThermalSource_H

#include "solverObject/solverObject.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "basicThermo/basicThermo.H"
#include "turbulenceModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "concentrationTransport/concentrationTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class humidityThermalSource Declaration
\*---------------------------------------------------------------------------*/


class humidityThermalSource
:
    public solverObject
{

protected:

    //- Name of temperature field
    word heTName_;

    //- Boundary source field name
    word TName_;

    //- Name of the humidity field
    word humidityName_;

    //- Name of concentration transport object
    word transportName_;

    //- Pointer to species concentration solver object
    concentrationTransport* cTransport_;

    //- Pointer to thermo
    const fluidThermo* thermo_;

    //- Pointer to transport model
    const incompressible::turbulenceModel* transport_;

    //- Pointer to humidity field
    const volScalarField* w_;

    //- Cpw
    autoPtr<dimensionedScalar> Cpw_;

    //- Returns humidity mass flux
    virtual tmp<surfaceScalarField> phiVDiff();


public:

    //- Runtime type information
    TypeName("humidityThermalSource");


    // Constructors

    //- Construct from explicit source name and mesh
    humidityThermalSource
    (
        const word& sourceName,
        const objectRegistry& obr,
        const dictionary& dict
    );

    //- Destructor
    virtual ~humidityThermalSource()
    {}

    // Member Functions

    //- Function is called after all solver objects are constructed
    virtual bool initialise() override;

    virtual void getSolveGraph
    (
        wordList& solveNames,
        HashTable<wordList>& requiredDependencies,
        HashTable<wordList>& optionalDependencies,
        HashTable<wordList>& correctorMembers
    ) override
    {}

    virtual void getSourceGraph
    (
        wordList& fieldNames,
        HashTable<wordList>& sourceDependencies
    ) override;

    //- Add contribution to equation
    virtual void addSup
    (
        fvMatrix<scalar>& eqn,
        const label fieldI
    ) override;

    //- Add contribution to compressible equation
    virtual void addSup
    (
        const volScalarField& rho,
        fvMatrix<scalar>& eqn,
        const label fieldI
    ) override;

    virtual void getBoundarySourceGraph
    (
        HashTable<labelList>& fieldPatchIDs,
        HashTable<wordList>& boundarySourceDependencies
    ) override;

    virtual void addBoundarySource
    (
        const word& fieldName,
        const label patchID,
        const scalarField& pf,
        scalarField& f,
        scalarField& df
    ) override;
};

} // End namespace fv
} // End namespace Foam

#endif
