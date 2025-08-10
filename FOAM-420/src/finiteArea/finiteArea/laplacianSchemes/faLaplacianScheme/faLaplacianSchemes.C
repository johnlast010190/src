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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description
    Abstract base class for finite area calculus laplacian schemes.

\*---------------------------------------------------------------------------*/

#include "finiteArea/laplacianSchemes/faLaplacianScheme/faLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define the constructor function hash tables
#define makeLaplacianGTypeScheme(Type, GType)                                 \
    typedef laplacianScheme<Type, GType> laplacianScheme##Type##GType;        \
    defineTemplateRunTimeSelectionTable(laplacianScheme##Type##GType, Istream);

#define makeLaplacianScheme(Type)                                             \
    makeLaplacianGTypeScheme(Type, scalar);                                   \
    makeLaplacianGTypeScheme(Type, symmTensor);                               \
    makeLaplacianGTypeScheme(Type, tensor);

makeLaplacianScheme(scalar);
makeLaplacianScheme(vector);
makeLaplacianScheme(sphericalTensor);
makeLaplacianScheme(symmTensor);
makeLaplacianScheme(tensor);

//defineTemplateRunTimeSelectionTable
//(
//    laplacianScheme<scalar>,
//    Istream
//);

//defineTemplateRunTimeSelectionTable
//(
//    laplacianScheme<vector>,
//    Istream
//);

//defineTemplateRunTimeSelectionTable
//(
//    laplacianScheme<tensor>,
//    Istream
//);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
