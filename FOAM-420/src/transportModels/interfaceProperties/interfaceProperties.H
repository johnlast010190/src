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
    (c) 2017 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "surfaceTensionModels/surfaceTensionModel/surfaceTensionModel.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        surfaceScalarField cAlpha_;

        //- Compression coefficient model
        word cAlphaModel_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Threshold for interface search used in
        //  interface Courant number calculation
        dimensionedScalar alphaCoThreshold_;

        //- Smoothing alpha pre interface normal calculation
        const label nSmoothAlpha_;

        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nHatf_;
        volScalarField K_;
        autoPtr<volVectorField> gradAlphaPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfaceProperties(const interfaceProperties&);
        void operator=(const interfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();

        //- Re-calculate compression coeff
        void calculateCAlpha();


public:

    // To enable the use of debug switch
    TypeName("interfaceProperties");

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );

    // Destructor
        virtual ~interfaceProperties(){ }


    // Member Functions

        const surfaceScalarField& cAlpha()
        {
            calculateCAlpha();
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const label& nSmoothAlpha() const
        {
            return nSmoothAlpha_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        tmp<volScalarField> sigmaK() const;

        const volVectorField& gradAlpha() const
        {
            return gradAlphaPtr_();
        }

        tmp<volVectorField> nHat() const
        {
            return (gradAlphaPtr_()/(mag(gradAlphaPtr_()) + deltaN_));
        }

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Interface indicator used for calculating
        //  interface Courant number (user-defined threshold)
        tmp<volScalarField> nearInterfaceCo() const;

        void correct();

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
