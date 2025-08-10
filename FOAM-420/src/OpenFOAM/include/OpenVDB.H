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
    (c) 2021 Esi Ltd.

Namespace
    OpenVDB

Description
    OpenVDB headers

\*---------------------------------------------------------------------------*/

#ifndef OpenVDB_H
#define OpenVDB_H

//clash with macro expansion in <openvdb_dir>/include/openvdb/math/Vec3.h:663
#undef Log
#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>
#define Log if (log) Info

using FloatGrid = openvdb::Grid<openvdb::tree::Tree4<float, 5, 4, 3>::Type>; //default
using Vec3SGrid = openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec3f, 5, 4, 3>::Type>; //default
using Vec3dGrid = openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec3d, 5, 4, 3>::Type>; //default
using BoolGrid = openvdb::Grid<openvdb::tree::Tree4<bool, 5, 4, 3>::Type>; //default
//using FloatGrid = openvdb::Grid<openvdb::tree::Tree4<float, 1, 1, 1>::Type>;
//using BoolGrid = openvdb::Grid<openvdb::tree::Tree4<bool, 1, 1, 1>::Type>;

using FaceGrid  = openvdb::Grid<openvdb::tree::Tree4<openvdb::math::Vec4ui, 5, 4, 3>::Type>;
using Vec4IGrid = openvdb::Grid<openvdb::tree::Tree4<openvdb::math::Vec4i, 5, 4, 3>::Type>;
using IndexGrid = openvdb::Int32Grid; //openvdb::Int64Grid;

#endif

// ************************************************************************* //
