/*---------------------------------------------------------------------------* \
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
    (c) 2019-2021 Esi Ltd.
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::fv::electrostaticPotentialSolver

Group
    grpFvOptionsSolvers

Description
    Solves an electrical potential equation

    \f[
        \grad \left( \sigma \grad V \right)
    \f]

    where \f$ V \f$ is electrical potential and \f$\sigma\f$ is the
    electrical conductivity.

Usage
    Isotropic (scalar) electrical conductivity
    \verbatim
    electrostaticPotentialSolver
    {
        type electrostaticPotentialSolver;

        anisotropicConductivity no;

        // Optionally specify the conductivity as a function of temperature
        sigma           table
        (
            (273        1e5)
            (1000       1e5)
        );
    }
    \endverbatim

    Anisotropic (vectorial) electrical conductivity
    electrostaticPotentialSolver
    {
        type electrostaticPotentialSolver;

        anisotropicConductivity yes;

        sigma           (31900 63800 127600);

        //sigma           table
        //(
        //    (0      (0 0 0))
        //    (1000   (127600 127600 127600))
        //);

        coordinateSystem
        {
            type        cartesian;
            origin      (0 0 0);

            coordinateRotation
            {
                type        axesRotation;
                e1          (1 0 0);
                e3          (0 0 1);
            }
        }
    }


    Where:
    \table
        Property     | Description               | Required  | Default value
        T            | Name of temperature field | no        | T
        sigma        | Electrical conductivity as a function of temperature |no|
        anisotropicElectricalConductivity | Anisotropic flag | yes |
    \endtable

    The electrical conductivity can be specified using either:
    - If the \c sigma entry is present the electrical conductivity is specified
      as a function of temperature using a Function1 type
    - If not present the sigma field will be read from file
    - If the anisotropicElectricalConductivity flag is set to 'true', sigma
      should be specified as a vector quantity


SourceFiles
    electrostaticPotentialSolver.C

SeeAlso
    Foam::fv::jouleHeatingSource

\*---------------------------------------------------------------------------*/

#ifndef electrostaticPotentialSolver_H
#define electrostaticPotentialSolver_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "coordinate/systems/coordinateSystem.H"
#include "referenceFrames/coordinateFrame.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class electrostaticPotentialSolver Declaration
\*---------------------------------------------------------------------------*/

class electrostaticPotentialSolver
:
    public solverObject
{

protected:

    const solutionControl* solnControlPtr_;

    //- Electrical potential field / [V]
    volScalarField V_;

    //- Flag to indicate that the electrical conductivity is anisotropic
    bool anisotropic_;

    //- Electrical conductivity as a scalar function of temperature
    autoPtr<Function1<scalar>> sigmaFn_;

    //- Isotropic conductivity field
    autoPtr<volScalarField> sigma_;

    //- Electrical conductivity as a vector function of temperature
    autoPtr<Function1<vector>> aniSigmaFn_;

    //- Anisotropic conductivity field in local co-ordinates
    autoPtr<volVectorField> aniLocalSigma_;

    //- Anisotropic conductivity field
    autoPtr<volSymmTensorField> aniSigma_;

    //- Co-ordinate system - used for anisotropic conductivity
    autoPtr<coordinateSystem> coordinates_;

    //- Pointer to reference frame (will allow updates as well)
    coordinateFrame* coorFramePtr_;


protected:

    //- Return the co-ordinate system for anisotropic electrical
    //  conductivity
    const coordinateSystem& coordSys() const;

    //- Transform the anisotropic electrical conductivity into global system
    void transformSigma
    (
        const volVectorField& sigmaLocal,
        volSymmTensorField& aniSigma
    );

    //- Initialise the electrical conductivity field
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>>
    createSigma
    (
        autoPtr<Function1<Type>>& sigmaVsTPtr
    );

    //- Update local conductivity field
    template<class Type>
    void updateSigma
    (
        const autoPtr<Function1<Type>>& sigmaVsTPtr,
        GeometricField<Type, fvPatchField, volMesh>& sigma
    ) const;

    //- Update conductivity field
    void updateSigma();

public:

//- Runtime type information
    TypeName("electrostaticPotentialSolver");

// Constructors

    // - Construct from components
    electrostaticPotentialSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

    electrostaticPotentialSolver(electrostaticPotentialSolver&) = delete;

    void operator=(electrostaticPotentialSolver&) = delete;

//- Destructor
    virtual ~electrostaticPotentialSolver()
    {}

// Member Functions

    virtual bool initialise() override;

    virtual void getSolveGraph
    (
        wordList& solveNames,
        HashTable<wordList>& requiredDependencies,
        HashTable<wordList>& optionalDependencies,
        HashTable<wordList>& correctorMembers
    ) override;

    virtual bool isFinalCorrector
    (
        const label corrector,
        const word& correctorName
    ) override;

    virtual tmp<fvScalarMatrix> assembleScalarMatrix
    (
        const word& fieldName,
        bool& finalSolve,
        word& dictName
    ) override;

    virtual void correct(const word&, const word&) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "electrostaticPotentialSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
