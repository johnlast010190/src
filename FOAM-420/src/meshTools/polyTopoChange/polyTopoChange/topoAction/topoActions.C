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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "polyTopoChange/polyTopoChange/topoAction/topoAction.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyPoint.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyCell.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoAction, 0);

    defineTypeNameAndDebug(polyAddPoint, 0);
    defineTypeNameAndDebug(polyModifyPoint, 0);
    defineTypeNameAndDebug(polyRemovePoint, 0);

    defineTypeNameAndDebug(polyAddFace, 0);
    defineTypeNameAndDebug(polyModifyFace, 0);
    defineTypeNameAndDebug(polyRemoveFace, 0);

    defineTypeNameAndDebug(polyAddCell, 0);
    defineTypeNameAndDebug(polyModifyCell, 0);
    defineTypeNameAndDebug(polyRemoveCell, 0);
}


// ************************************************************************* //
