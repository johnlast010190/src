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
    Foam::curvatureLevelSet

SourceFiles
    curvatureLevelSet.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Vignesh Thammanna Gurumurthy < vignesh@sla.tu-darmstadt.de >

    All rights reserved.


S-SCLSVOF - Density Balanced - Method Developement:

    Albadawi A, Donoghue DB, Robinson AJ, Murray DB, Delauré YMC.
    Influence of surface tension implementation in volume of fluid and coupled volume of fluid with level set methods for bubble growth and detachment. International
    Journal of Multiphase Flow 2013; 53:11–28.

    T. YAMAMOTO, Y. OKANO AND S. DOST
    Validation of the S-CLSVOF method with the density-scaled balanced continuum surface force model in multiphase systems coupled with thermocapillary flows
    Int. J. Numer. Meth. Fluids 2017; 83:223-244.

   Yokoi K.
    A density-scaled continuum surface force model within a balanced force formulation. Journal of
    Computational Physics 2014; 278:221–228.

Description

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


#include "curvature/curvatureLevelSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug( curvatureLevelSet, 0 );
addToRunTimeSelectionTable( curvature, curvatureLevelSet, dictionary );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvatureLevelSet::curvatureLevelSet
(
    const word& name,
    const volScalarField& alpha,
    const volVectorField& nHatv,
    const surfaceVectorField& nHatfv,
    const volScalarField& interfaceDensity,
    const List<bool>& isWallPatch,
    const dictionary& transpProp,
    const volScalarField& isInterface
)
    :
    curvature( name, alpha, nHatv, nHatfv, interfaceDensity, isWallPatch, transpProp, isInterface),

    alpha1_(alpha),

    delta_(alpha.mesh().objectRegistry::lookupObject<const volScalarField> ("delta")),

    Hscaled_(alpha.mesh().objectRegistry::lookupObject<const volScalarField> ("Hscaled")),

    psi_(alpha.mesh().objectRegistry::lookupObject<const volScalarField> ("psi")),

    scaled_(transpProp.subDict("reconstruct").lookup("scaled"))
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void curvatureLevelSet::calculateK()
{

    kappa_ = - fvc::div(nHatfv_ & mesh_.Sf())*isInterface_;
}

tmp<surfaceScalarField> curvatureLevelSet::surfaceTensionForce(dimensionedScalar sigma) const
{

    //- if surface tension force is scaled to reduce density ration dependency
    //  of suprious currents
    if (scaled_)
    {
        return    fvc::interpolate(sigma*kappa_)
              //  * fvc::snGrad(psi_)*fvc::interpolate(delta_*isInterface_);
                * fvc::snGrad(Hscaled_);
    }

    return    fvc::interpolate(sigma*kappa_)
            * fvc::snGrad(psi_)
            * fvc::interpolate(delta_*isInterface_);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
