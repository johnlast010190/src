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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Class
    Foam::porosityModels::fixedCoeff

Description
    Fixed coefficient form of porosity model

        \f[
            S = - \rho_ref (\alpha + \beta |U|) U
        \f]

    In the case of compressible flow, a value for the reference density is
    required

SourceFiles
    fixedCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef fixedCoeff_H
#define fixedCoeff_H

#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "dimensionedTypes/dimensionedTensor/dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                          Class fixedCoeff Declaration
\*---------------------------------------------------------------------------*/

class fixedCoeff
:
    public porosityModel
{
protected:
    // Protected data

        //- Alpha coefficient XYZ components (user-supplied) [1/s]
        dimensionedVector alphaXYZ_;

        //- Beta coefficient XYZ components (user-supplied) [1/m]
        dimensionedVector betaXYZ_;

        //- Model alpha coefficient - converted from alphaXYZ [1/s]
        List<tensorField> alpha_;

        //- Model beta coefficient - converted from betaXYZ [1/m]
        List<tensorField> beta_;

private:
    // Private Member Functions

        //- Apply
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const vectorField& U,
            const scalar rho
        ) const;

        //- Apply
        void apply
        (
            tensorField& AU,
            vectorField& source,
            const vectorField& U,
            const scalar rho
        ) const;

        //- Disallow default bitwise copy construct
        fixedCoeff(const fixedCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const fixedCoeff&);


public:

    //- Runtime type information
    TypeName("fixedCoeff");

    //- Constructor
    fixedCoeff
    (
        const word& name,
        const word& modelType,
        const objectRegistry& obr,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~fixedCoeff();


    // Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData();

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const;

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;
        virtual void correct(fvBlockMatrix<vector>& UEqn) const;

        //- Add resistance
        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        virtual void correct
        (
            fvBlockMatrix<vector>& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;

        //- Add adjoint resistance
        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volVectorField& Uprimal
        ) const;

        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        ) const;

        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volVectorField&
        ) const;

        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&,
            const volVectorField&
        ) const;

    // I-O

        //- Write
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
