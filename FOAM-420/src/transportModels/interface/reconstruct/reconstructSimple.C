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
    reconstruct

SourceFiles
    reconstruct.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    All rights reserved.

Description

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de> (main developer).

    Method Development and Intellectual Property :
        Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
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


#include "reconstruct/reconstructSimple.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{

    defineTypeNameAndDebug( reconstructSimple, 0 );
    addToRunTimeSelectionTable( reconstruct, reconstructSimple, dictionary );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reconstructSimple::reconstructSimple
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& transpProp,
    const List<bool>& isWallPatch,
    const volScalarField& isInterface
)
:
    reconstruct( name, alpha, transpProp, isWallPatch, isInterface )
{
    reconstructInterface();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- reconstruct a smeared interface without explicite representation
//  interFoam standart procedure to optain the interface normal nHatv
//  the interfaceDensity is calculated by
void reconstructSimple::reconstructInterface()
{
    //- Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

    //- Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    //interface normal
    nHatfv_ = gradAlphaf/(mag(gradAlphaf) + deltaN_);
    nHatv_  = gradAlpha/(mag(gradAlpha) + deltaN_);

    //- calculate interface density
    interfaceDensity_ = mag(gradAlpha); //TODO is this correct?
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
