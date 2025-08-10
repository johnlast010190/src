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
    (c) 2020 Esi Ltd.

Class
    Foam::sweptPatchToCell

Description
    A topoSetSource to select cells based on the rotation of a patch around an
    axis

Example
    {
        name setZoneName;
        type cellZoneSet;
        action new;
        source sweptPatchToCell;
        sourceInfo
        {
            patches         ( patchName );
            origin          ( 0 0 0 );
            axis            ( 0 0 1 );
        }
    }
SourceFiles
    sweptPatchToCell.C

\*---------------------------------------------------------------------------*/

#ifndef sweptPatchToCell_H
#define sweptPatchToCell_H

#include "sets/topoSetSource/topoSetSource.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "cfdTools/general/include/fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sweptPatchToCell Declaration
\*---------------------------------------------------------------------------*/
class sweptPatchToCell
:
    public topoSetSource
{

    fvMesh mesh_;

    wordList patchNames_;
    labelList patchLabels_;

    autoPtr<volScalarField> tracerFieldPtr_;
    autoPtr<surfaceScalarField> phiMRFPtr_;

    boolList cellBelongsToSet_;

    vector axisOfRotation_;
    vector origin_;

    scalar tracerThresholdValue_;

    label numberOfSolvingIterations_;

    // Private Functions
    void calculateTracerField();

    void setFluxNameInBoundaryPatches();
    void setBoundaryValueInSweptPatches();
    void constructRotationalFlux();

    void setFvSchemesDict();
    void setFvSolutionDict();

    void solveConvectionEquation();
    void correctZeroDiagonalElements(fvScalarMatrix& eqnMatrix);

    void reverseFluxAndInitializeTracer();

    void findCellsInSet();
    void excludeCellsWithValuesSmallerThanThreshold();
    void excludeCellsFurtherThanMaximumDistance();
    void forceAddBoundaryCellsOfSweptPatches();
    scalar getMaximumDistanceFromTheAxis();

    static addToUsageTable usage_;

    // Private Member Functions

    void combine(topoSet& set, const bool add) const;


 public:

    //- Runtime type information
    TypeName("sweptPatchToCell");

    //- Construct from dictionary
    sweptPatchToCell
    (
        const polyMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~sweptPatchToCell();


    // Member Functions

    virtual sourceType setType() const
    {
        return CELLSETSOURCE;
    }

    virtual void applyToSet
    (
        const topoSetSource::setAction action,
        topoSet&
    ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
