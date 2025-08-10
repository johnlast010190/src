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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    controlBox

Description
    Defines a control lattice

\*---------------------------------------------------------------------------*/

#ifndef controlBox_H
#define controlBox_H

#include "db/dictionary/dictionary.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "fields/Fields/vectorField/vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class controlBox Declaration
\*---------------------------------------------------------------------------*/

class controlBox
{
    // Private data

    //- number of control points
    label nXCPs_, nYCPs_, nZCPs_;

    //points with min and max coords
    point minpt_, maxpt_;

    //- point positions and displacements
    //- Indexing goes as follows: for point i, j, k,  p = points(index)
    //- index = i + j*nXCPs_ + k*nXCPs_*nYCPs_
    pointField points_;

    vectorField displacements_;

    // polymesh representation
    polyMesh *cBoxMesh;

    // private functions
    face createFace
    (
        const label& index0,
        const label& index1,
        const label& index2,
        const label& index3
    ) const;

    void init(const pointField& blockPoints);

    // calculates the box
    void computePoints(const pointField& blockPoints);

    // calculate index of a point
    label calculateIndex
    (
        const label  iCP,
    const label  jCP,
    const label  kCP
    ) const;

public:

    // Constructors

    //- Construct from dictionary
    controlBox
    (
        const dictionary& dict
    );

    //- Construct as copy
    controlBox(const controlBox& box);


    // Destructor

    ~controlBox()
    {
        delete cBoxMesh;

        cBoxMesh = nullptr;
    };


    // Member Functions

    const point& getPoint
    (
        const label  iCP,
        const label  jCP,
        const label  kCP
    ) const;

    const vector& getDisplacement
    (
        const label  iCP,
    const label  jCP,
       const label  kCP
    ) const;

    void updateBoxPosition();

    void setDisplacement
    (
        const label  iCP,
    const label  jCP,
        const label  kCP,
    const vector dV
    );

    void setDisplacement
    (
        const label  index,
    const vector dV
    );

    void addDisplacement
    (
        const label  index,
    const vector dV
    );

    void boxToPolyMesh(const IOobject& io);

    const polyMesh* getPolyMesh() const;

    void smoothDisplacements();

    const pointField& points() const
    {
        return points_;
    }

    const vectorField& displacements() const
    {
        return displacements_;
    }

    inline label nbXCPs() const { return nXCPs_;}
    inline label nbYCPs() const { return nYCPs_;}
    inline label nbZCPs() const { return nZCPs_;}

    inline point getMaxPt() const { return maxpt_;}
    inline point getMinPt() const { return minpt_;}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
