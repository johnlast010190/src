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
    (c) 2019 Esi Ltd.
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Class
    Foam::functionObjects::scalarTransport

Group
    grpSolversFunctionObjects

Description
    Evolves a passive scalar transport equation.

    - To specify the field name set the \c field entry
    - To employ the same numerical schemes as another field set
      the \c schemesField entry,
    - The diffusivity can be set manually using the 'D' entry, retrieved
      from the turbulence model or specified nut
    - Alternatively if a turbulence model is available a turbulent diffusivity
      may be constructed from the laminar and turbulent viscosities using the
      optional diffusivity coefficients \c alphaD and \c alphaDt (which default
      to 1):
      \verbatim
          D = alphaD*nu + alphaDt*nut
      \endverbatim
    - To specify a transport quantity within a phase enter phase.
    - bounded01 bounds the transported scalar within 0 and 1.

Usage
    Example of solver object specification to solve a scalar transport
    equation:
    functions
    {
        scalar1
        {
            type            scalarTransport;
            libs            ("libsolverFunctionObjects.so");

            resetOnStartUp  no;
            region          cabin;
            field           H2O;
        }
    }
    \endverbatim

    Example of function object specification to solve a residence time
    in a two phase flow:
    equation:
    \verbatim
    functions
    {
        sTransport
        {
            type            scalarTransport;
            libs            ("libsolverFunctionObjects.so");

            enabled         true;
            writeControl    outputTime;
            writeInterval   1;

            field           s;
            bounded01       false;
            phase           alpha.water;

            write           true;

            resetOnStartUp  false;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: scalarTransport | yes      |
        field        | Name of the scalar field | no         | s
        phi          | Name of flux field      | no          | phi
        rho          | Name of density field   | no          | rho
        phase        | Name of the phase       | no          | none
        nut          | Name of the turbulence viscosity | no | none
        D            | Diffusion coefficient   | no          | auto generated
        nCorr        | Number of correctors    | no          | 0
        resetOnStartUp | Reset scalar to zero on start-up | no | no
        schemesField | Name of field to specify schemes | no | field name
        bounded01    | Bounds scalar between 0-1 for multiphase | no | true
        phasePhiCompressed | Compressed flux for VOF | no | alphaPhiUn
    \endtable

See also
    Foam::fv::solverObject

SourceFiles
    scalarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scalarTransport_H
#define functionObjects_scalarTransport_H

#include "solverObject/solverObject.H"
#include "fields/volFields/volFields.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class scalarTransport
:
    public solverObject
{
    // Private data

        //- Name of the transport field.
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Name of turbulent viscosity field (optional)
        word nutName_;

        //- Name of phase field (optional)
        word phaseName_;

        //- Name of phase field compressed flux (optional)
        word phasePhiCompressedName_;

        //- Diffusion coefficient (optional)
        scalar D_;

        //- Flag to indicate whether a constant, uniform D_ is specified
        bool constantD_;

        //- Laminar diffusion coefficient (optional)
        scalar alphaD_;

        //- Turbulent diffusion coefficient (optional)
        scalar alphaDt_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Flag to reset the scalar to zero on start-up
        bool resetOnStartUp_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Bound scalar between 0-1 using MULES for multiphase case
        bool bounded01_;

        //- Modify deltaT to speed-up passive transport (optional)
        label solveSpeedup_;

        //- Allow sub-cycling with specified number of time steps
        label nSubCycles_;

        //- Stabilize phasic transport (optional)
        dimensionedScalar residualPhaseField_;

        //- Used to specify execution order in USF framework
        word fieldDependency_;

        //- Include ddt term y/n
        bool transient_;

        //- Correct for formal continuity error in transient and steady-state y/n
        bool bounded_;

        //- Force wall diffusion to be zero for conservation y/n
        bool zeroWallDiffusion_;

        //- use old-style phasic transport (more stable, less accurate) y/n
        bool stablePhasic_;

        //- use phasic fvOption sources y/n
        bool phasicSources_;


    // Private Member Functions

        //- Return reference to registered transported field
        volScalarField& transportedField();

        //- Return the diffusivity field
        tmp<volScalarField> D
        (
            const volScalarField& s,
            const surfaceScalarField& phi
        ) const;

        //- Disallow default bitwise copy construct
        scalarTransport(const scalarTransport&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const scalarTransport&) = delete;


public:

    //- Runtime type information
    TypeName("scalarTransport");


    // Constructors

        //- Construct from Time and dictionary
        scalarTransport
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    //- Destructor
    virtual ~scalarTransport();


    // Member Functions

        //- Read the scalarTransport data
        virtual void read(const dictionary&) override;

        //- Identify any solves handled by this solverObject, their levels of
        // iteration and dependent solves
        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override
        {
            solveNames.append(fieldName_);

            if (fieldDependency_ != "none")
            {
                requiredDependencies.insert(fieldName_, {fieldDependency_});
            }
            optionalDependencies.insert(fieldName_, {"fvMesh"});
            correctorMembers.insert
            (
                solverObject::outerCorrectorName, solveNames
            );
        }

        //- Calculate the scalarTransport
        virtual void correct(const word&, const word&) override;

        //- assemble and solve inisde separate function since
        //  optional sub-cycling otherwise duplicates all code
        virtual void solveEquation
        (
            volScalarField& s,
            const surfaceScalarField& phi,
            volScalarField& D,
            fv::options& fvOptions
        );

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual void write() override;

        //- Store the field in the (sub) objectRegistry under the given name
        //  Note: sets the fieldName to tfield().name() if not already set
        template<class ObjectType>
        bool store
        (
            word& fieldName,
            const tmp<ObjectType>& tfield,
            bool cacheable = false
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
