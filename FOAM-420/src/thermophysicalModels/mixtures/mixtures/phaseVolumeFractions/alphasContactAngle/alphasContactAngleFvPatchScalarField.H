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

Class
    Foam::alphasContactAngleFvPatchScalarField

Description
    Contact-angle boundary condition for multi-phase interface-capturing
    simulations.

SourceFiles
    alphasContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphasContactAngleFvPatchScalarField_H
#define alphasContactAngleFvPatchScalarField_H

#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "interfacePair/interfacePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class alphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphasContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    class interfaceThetaProps
    {
        //- Equilibrium contact angle
        scalar theta0_;

        //- Dynamic contact angle velocity scale
        scalar uTheta_;

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receeding contact angle
        scalar thetaR_;


    public:

        // Constructors
            interfaceThetaProps()
            {}

            interfaceThetaProps(Istream&);


        // Member functions

            //- Return the equilibrium contact angle theta0
            scalar theta0(bool matched=true) const
            {
                if (matched) return theta0_;
                else return 180.0 - theta0_;
            }

            //- Return the dynamic contact angle velocity scale
            scalar uTheta() const
            {
                return uTheta_;
            }

            //- Return the limiting advancing contact angle
            scalar thetaA(bool matched=true) const
            {
                if (matched) return thetaA_;
                else return 180.0 - thetaA_;
            }

            //- Return the limiting receeding contact angle
            scalar thetaR(bool matched=true) const
            {
                if (matched) return thetaR_;
                else return 180.0 - thetaR_;
            }


        // IO functions

            friend Istream& operator>>(Istream&, interfaceThetaProps&);
            friend Ostream& operator<<(Ostream&, const interfaceThetaProps&);
    };

    typedef HashTable<interfaceThetaProps, interfacePair, interfacePair::hash>
        thetaPropsTable;


private:

    // Private data

        thetaPropsTable thetaProps_;


public:

    //- Runtime type information
    TypeName("alphasContactAngle");


    // Constructors

        //- Construct from patch and internal field
        alphasContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphasContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphasContactAngleFvPatchScalarField
        //  onto a new patch
        alphasContactAngleFvPatchScalarField
        (
            const alphasContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphasContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphasContactAngleFvPatchScalarField
        (
            const alphasContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphasContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return the contact angle properties
        const thetaPropsTable& thetaProps() const
        {
            return thetaProps_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
