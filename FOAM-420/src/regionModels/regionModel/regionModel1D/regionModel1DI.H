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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "regionModel1D/regionModel1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::labelListList&
Foam::regionModels::regionModel1D::boundaryFaceFaces() const
{
    return boundaryFaceFaces_;
}


inline const Foam::labelListList&
Foam::regionModels::regionModel1D::boundaryFaceCells() const
{
    return boundaryFaceCells_;
}


inline const Foam::labelList&
Foam::regionModels::regionModel1D::boundaryFaceOppositeFace() const
{
    return boundaryFaceOppositeFace_;
}


inline const Foam::surfaceScalarField&
Foam::regionModels::regionModel1D::nMagSf() const
{
    if (!nMagSfPtr_.valid())
    {
        FatalErrorInFunction
         << "Face normal areas not available" << abort(FatalError);
    }

    return nMagSfPtr_();
}


inline Foam::label Foam::regionModels::regionModel1D::nLayers() const
{
    return nLayers_;
}


// ************************************************************************* //
