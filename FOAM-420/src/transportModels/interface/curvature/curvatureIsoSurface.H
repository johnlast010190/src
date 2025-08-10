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
    Foam::curvatureIsoSurface

SourceFiles
    curvature.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    Christian Kunkelmann (2011)
    All rights reserved.

Description
    Curvature is calculated using the normals of an isoSurface reconstruction
    of the interface.
    Cells at wall boundary are treated explicitly (method by Kunkelmann )

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


#ifndef curvatureIsoSurface_H
#define curvatureIsoSurface_H

#include "curvature/curvature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class curvatureIsoSurface Declaration
\*---------------------------------------------------------------------------*/

class curvatureIsoSurface
:
    public curvature
{

private:


public:

    // Runtime type name
        TypeName("curvatureIsoSurface");

    // Constructors

        //- Construct from components
        curvatureIsoSurface
        (
            const word& name,
            const volScalarField& alpha,
            const volVectorField& nHatv,
            const surfaceVectorField& nHatfv,
            const volScalarField& interfaceDensity,
            const List<bool>& isWallPatch,
            const dictionary& transpProp,
            const volScalarField& isInterface
        );

    //- Destructor
    virtual ~curvatureIsoSurface(){}

    // Member Functions

        //- calculates the curvature
        void calculateK();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
