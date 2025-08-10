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
    (c) held by original author

Class
    Foam::interface

SourceFiles
    interface.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    All rights reserved.

Description
    Holds the all the interface properties.
        - Corrects mixed material properties
        - Marks cells for special interface handling
        - Does special contact line handling for cLPosition, cLVelocity and
          contact angle
        - Calculates interface representation
        - Calculates the interface curvature

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de> (main developer).

    Method Development and Intellectual Property :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef interface_H
#define interface_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

#include "reconstruct/reconstruct.H"
#include "curvature/curvature.H"

#include "fields/DimensionedFields/DimensionedScalarField/DimensionedScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interface Declaration
\*---------------------------------------------------------------------------*/

class interface
{
    // Private data
    //- Keep a reference to the transportProperties dictionary
    const dictionary& transpProp_;

    //- Compression coefficient
    surfaceScalarField cAlpha_;

    //- Compression coefficient model
    word cAlphaModel_;

    //- Surface tension
    dimensionedScalar sigma_;

    //- Stabilization for normalization of the interface normal
    const dimensionedScalar deltaN_;

    //- Threshold for interface search used in
    //  interface Courant number calculation
    dimensionedScalar alphaCoThreshold_;

    //- reference to the volume fraction field
    const volScalarField& alpha1_;

    //- reference to velocity field
    const volVectorField& U_;

    //- indicates if an alpha1_ patch is a wall patch
    List<bool> isWallPatch_;

    //- define interface via snGradAlpha or alpha itself
    word isInterfaceMethod_;

    //- user defined value
    scalar isInterfaceThreshold_;

    //- thicken the interface by N cells
    int isInterfaceAddN_;

    //- marker-field for interface cells (1 true, 0 false)
    volScalarField isInterface_;

    //- model dependent reconstruction of the interface
    //  holds nHatv and nHatfv
    autoPtr<reconstruct> reconI_;

    //- model dependent curvature calculation
    autoPtr<curvature> K_;

    // Private Member Functions
    List<bool> setIsWallPatch(const volScalarField& alpha) const;

    //- initializes isInterface in constructor call
    void initIsInterface();

    //- sets marker field for interface-near cells
    void calculateIsInterface(volScalarField& isInterface);

    //- disallow default bitwise copy construct and assignment
    interface(const interface&);
    void operator=(const interface&);

    //- re-calculate the interface curvature
    void calculateK(); // maybe move to separate class

    //- Correction for the boundary condition on the unit normal nHat on
    //  walls to produce the correct contact dynamic angle
    //  calculated from the component of U parallel to the wall
    void correctContactAngle
    (
        surfaceVectorField::Boundary& nHat,
        const surfaceVectorField::Boundary& gradAlphaf
    );

    //- Re-calculate compression coeff
    void calculateCAlpha();

    //- Initialise class data members
    void initialise();


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;

    // Constructors
        interface
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        const surfaceScalarField& cAlpha()
        {
            calculateCAlpha();
            return cAlpha_;

        }

        //returns a const reference of nHatfv & alpha1_.mesh().Sf()
        const surfaceScalarField& nHatf()
        {
            //sidenote: in reconstructIsoSurface the herefore used normal is
            //          calculated via gradAlpha, not the reconstructed one.
            return reconI_->nHatf();
        }

        //- returns a const reference of the curvature class
        const curvature& K() const
        {
            return K_;
        }

        //- returns a const reference of the surface tension
        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        //- Surfacetension * curvature
        tmp<volScalarField> sigmaK() const
        {
            return sigma_ * K_->kappa();
        }

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Interface indicator used for calculating
        //  interface Courant number (user-defined threshold)
        tmp<volScalarField> nearInterfaceCo() const;

        //- does the necessary stuff to update the interface in one time step
        void correct()
        {
            //- create marker field isInterface
            calculateIsInterface(isInterface_);

            //- calculate interface -> nHatv and nHatfv
            //  sets interfaceDensity
            reconI_->reconstructInterface();

            //- correct nHatf_ at boundaries
            surfaceVectorField gradAlphaf(fvc::interpolate(fvc::grad(alpha1_, "nHat")));
            correctContactAngle(reconI_->nHatfv().boundaryFieldRef(), gradAlphaf.boundaryField());

            //- calculate curvature
            K_->calculateK();
        }

        //- Read transportProperties dictionary
        bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
