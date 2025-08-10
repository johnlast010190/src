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
    (c) 2020 OpenFOAM Foundation
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicTransform.H"
#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<cyclicTransform::transformTypes, 4>::names[] =
    {
        "unspecified",
        "none",
        "rotational",
        "translational"
    };

    const NamedEnum<cyclicTransform::transformTypes, 4>
        cyclicTransform::transformTypeNames;

    template<>
    const char* NamedEnum<cyclicTransform::rotationalMethods, 3>::names[] =
    {
        "userSpecifiedAngle",
        "centreOfMassPatch",
        "furthestPoint"
    };

    const NamedEnum<cyclicTransform::rotationalMethods, 3>
        cyclicTransform::rotationalMethodNames;

    defineTypeNameAndDebug(cyclicTransform, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::cyclicTransform::readTransform
(
    const dictionary& dict,
    const bool defaultIsNone
)
{
    word transformEntry;

    if (dict.found("transformType"))
    {
        // Lookup the new "transformType" keyword first
        transformEntry = dict.lookup<word>("transformType");
    }
    else if (dict.found("transform"))
    {
        DeprecationWarningInFunction("transform", "keyword", 40000)
            << "Please use 'transformType' instead." << nl;

        // For backward compatibility, try to lookup "transform"
        transformEntry = dict.lookup<word>("transform");
    }

    // For backward compatibility, support the deprecated 'unknown'
    // transformation type, now superseded by 'unspecified'
    if (transformEntry == "unknown")
    {
        DeprecationWarningInFunction("unknown", "transformType", 40000)
            << "Please use 'unspecified' instead." << nl;

        transformEntry = "unspecified";
    }

    // For backward compatibility, support the deprecated 'noOrdering'
    // transformation type, now superseded by 'none'
    if (transformEntry == "noOrdering")
    {
        DeprecationWarningInFunction("noOrdering", "transformType", 40000)
            << "Please use 'none' instead." << nl;

        transformEntry = "none";
    }

    // If no transformation entry is found, assume the default transform type
    if (transformEntry == word::null)
    {
        transformEntry = defaultIsNone ? "none" : "unspecified";
    }

    return transformEntry;
}


void Foam::cyclicTransform::update()
{
    if (!transformComplete_)
    {
        return;
    }

    switch (transformType_)
    {
        case UNSPECIFIED:
            break;

        case NONE:
            transform_ = transformer();
            break;

        case ROTATIONAL:
            if (rotationAngle_ == 0)
            {
                transform_ = transformer();
            }
            else
            {
                const tensor R =
                    quaternion(rotationAxis_, degToRad(rotationAngle_)).R();

                if (mag(rotationCentre_) == 0)
                {
                    transform_ = transformer::rotation(R);
                }
                else
                {
                    transform_ =
                        transformer::translation(rotationCentre_)
                      & transformer::rotation(R)
                      & transformer::translation(- rotationCentre_);
                }
            }
            break;

        case TRANSLATIONAL:
            if (mag(separation_) == 0)
            {
                transform_ = transformer();
            }
            else
            {
                transform_ = transformer::translation(separation_);
            }
            break;
    }
}


bool Foam::cyclicTransform::set
(
    const cyclicTransform& t,
    const scalar lengthScale,
    const scalar matchTolerance
)
{
    // If the supplied transform is unspecified then there is nothing to do
    if (t.transformType_ == UNSPECIFIED)
    {
        return true;
    }

    // If this transform is specified then we need to check that it is
    // compatible with the supplied transform
    if (transformType_ != UNSPECIFIED)
    {
        // If the transforms are of different types then they are incompatible
        if (transformType_ != t.transformType_)
        {
            return false;
        }

        // Length-tolerance
        const scalar lengthTolerance = lengthScale*matchTolerance;

        // If the transforms are both rotational then the axes must be in the
        // same direction, the centre points must lie on the same line, and the
        // angles (if available) must be opposing.
        if (transformType_ == ROTATIONAL)
        {
            const scalar dot = rotationAxis_ & t.rotationAxis_;

            if (mag(dot) < 1 - matchTolerance)
            {
                return false;
            }

            if
            (
                (rotationAxis_ & (rotationCentre_ - t.rotationCentre_))
              > lengthTolerance
            )
            {
                return false;
            }

            if (transformComplete_ && t.transformComplete_)
            {
                if
                (
                    mag(degToRad(rotationAngle_ - sign(dot)*t.rotationAngle_))
                  > matchTolerance
                )
                {
                    return false;
                }
            }
        }

        // If the transforms are both translational then the separations must
        // be opposing
        if (transformType_ == TRANSLATIONAL)
        {
            if (transformComplete_ && t.transformComplete_)
            {
                if (mag(separation_ - t.separation_) > lengthTolerance)
                {
                    return false;
                }
            }
        }
    }

    // If the supplied transform is more complete then overwrite this with it
    if (!transformComplete_ && t.transformComplete_)
    {
        *this = t;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicTransform::cyclicTransform()
:
    cyclicTransform(true)
{}


Foam::cyclicTransform::cyclicTransform
(
    const bool defaultIsNone
)
:
    transformType_(defaultIsNone ? NONE : UNSPECIFIED),
    rotationalMethod_(CENTREOFMASS),
    rotationAxis_(vector::uniform(NaN)),
    rotationCentre_(vector::uniform(NaN)),
    rotationAngle_(NaN),
    separation_(vector::uniform(NaN)),
    transformComplete_(transformType_ == NONE),
    transform_()
{}


Foam::cyclicTransform::cyclicTransform
(
    const dictionary& dict,
    const bool defaultIsNone
)
:
    transformType_(transformTypeNames[readTransform(dict, defaultIsNone)]),
    rotationalMethod_
    (
        rotationalMethodNames.readOrDefault
        (
            dict.lookupOrDefault<word>
            (
                "rotationalMethod",
                "centreOfMassPatch"
            ),
            CENTREOFMASS
        )
    ),
    rotationAxis_
    (
        transformType_ == ROTATIONAL
      ? normalised(dict.lookup<vector>("rotationAxis"))
      : vector::uniform(NaN)
    ),
    rotationCentre_
    (
        transformType_ == ROTATIONAL
      ? dict.lookup<point>("rotationCentre")
      : point::uniform(NaN)
    ),
    rotationAngle_
    (
        rotationalMethod_ == USERSPECIFIED
      ? dict.lookup<scalar>("rotationAngle")
      : dict.lookupOrDefault<scalar>("rotationAngle", NaN)
    ),
    separation_
    (
        transformType_ == TRANSLATIONAL
      ? (
            // Lookup the new "separation" keyword first, then for backwards
            // compatibility try the deprecated "separationVector", then finally
            // spit an error that the "separation" keyword is not present.
            dict.found("separation")
          ? dict.lookup<vector>("separation")
          : dict.found("separationVector")
          ? dict.lookup<vector>("separationVector")
          : dict.lookup<vector>("separation")
        )
      : vector::uniform(NaN)
    ),
    transformComplete_
    (
        transformType_ == NONE
     || (transformType_ == ROTATIONAL && dict.found("rotationAngle"))
     || (
            transformType_ == TRANSLATIONAL
         && (dict.found("separation") || dict.found("separationVector"))
        )
    ),
    transform_()
{
    // If the deprecated "separationVector" keyword is present, a different
    // notation for the separation calculation is assumed.
    if (dict.found("separationVector"))
    {
        DeprecationWarningInFunction("separationVector", "keyword", 40200)
            << "Please use 'separation' instead and flip the vector sign."
            << nl;

        separation_ = -separation_;
    }

    if (transformComplete_)
    {
        update();
    }
}


Foam::cyclicTransform::cyclicTransform
(
    const word& name,
    const vectorField& areas,
    const cyclicTransform& transform,
    const word& nbrName,
    const cyclicTransform& nbrTransform,
    const scalar matchTolerance
)
:
    cyclicTransform(transform)
{
    // Calculate patch length scale
    const scalar lengthScale = sqrt(sum(mag(areas)));

    // Copy as much data from the neighbour as possible
    if (!transformComplete_ && nbrTransform.transformType_ != UNSPECIFIED)
    {
        if (!set(inv(nbrTransform), lengthScale, matchTolerance))
        {
            FatalErrorInFunction
                << "Patch " << name
                << " and neighbour patch " << nbrName
                << " have incompatible transforms:" << nl << nl << incrIndent;

            FatalErrorInFunction
                << indent << name << nl << indent << token::BEGIN_BLOCK << nl
                << incrIndent;

            cyclicTransform::write(FatalError);

            FatalErrorInFunction
                << decrIndent << indent << token::END_BLOCK << nl << nl;

            FatalErrorInFunction
                << indent << nbrName << nl << indent << token::BEGIN_BLOCK
                << nl << incrIndent;

            nbrTransform.write(FatalError);

            FatalErrorInFunction
                << decrIndent << indent << token::END_BLOCK << nl << nl;

            FatalErrorInFunction
                << decrIndent << exit(FatalError);
        }
    }
}


Foam::cyclicTransform::cyclicTransform
(
    const word& name,
    const pointField& ctrs,
    const vectorField& areas,
    const cyclicTransform& transform,
    const word& nbrName,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const cyclicTransform& nbrTransform,
    const scalar matchTolerance
)
:
    cyclicTransform
    (
        name,
        areas,
        transform,
        nbrName,
        nbrTransform,
        matchTolerance
    )
{
    // If there is no geometry from which to calculate the transform then
    // nothing can be calculated
    const label areasSize = returnReduce(areas.size(), plusOp<label>());
    const label nbrAreasSize = returnReduce(nbrAreas.size(), plusOp<label>());
    if (areasSize == 0 || nbrAreasSize == 0) return;

    // Calculate the total (vector) areas for the supplied patch data
    const vector area = gSum(areas);
    const vector nbrArea = gSum(nbrAreas);

    // Calculate the centroids for the supplied patch data
    const scalarField magAreas(mag(areas));
    const scalarField magNbrAreas(mag(nbrAreas));
    const point ctr = gSum(ctrs*magAreas)/gSum(magAreas);
    const point nbrCtr = gSum(nbrCtrs*magNbrAreas)/gSum(magNbrAreas);

    // Calculate patch length scale
    const scalar lengthScale = sqrt(gSum(magAreas));

    // Calculate the transformation from the patch geometry
    if (!transformComplete_)
    {
        // Store the old transformation type
        const transformTypes oldTransformType = transformType_;

        // Calculate the average patch normals
        const vector normal = normalised(area);
        const vector negNbrNormal = - normalised(nbrArea);

        // Calculate the angle and distance separations
        const scalar dot = normal & negNbrNormal;
        const vector delta = ctr - nbrCtr;

        // Determine the type of transformation if it has not been specified
        if (transformType_ == UNSPECIFIED)
        {
            transformType_ =
                dot < 1 - ROOTSMALL
              ? ROTATIONAL
              : mag(delta) > lengthScale*ROOTSMALL
              ? TRANSLATIONAL
              : NONE;
        }

        // If the transformation is known to be rotational, then we need to
        // calculate the angle according to the rotational method specified.
        // If the transformation was previously unspecified, then we also need
        // to calculate the axis and the centre of rotation.
        if (transformType_ == ROTATIONAL)
        {
            // Calculate the axis, if necessary
            if (transformType_ != oldTransformType)
            {
                const vector midNormal = normalised(normal + negNbrNormal);
                const vector axis =
                    (ctr - nbrCtr)
                  ^ (
                        normal*(negNbrNormal & midNormal)
                      - negNbrNormal*(normal & midNormal)
                    );
                const vector axis180 =
                    (ctr - nbrCtr)
                  ^ (normal - negNbrNormal);

                rotationAxis_ =
                    normalised
                    (
                        mag(axis) > lengthScale*ROOTSMALL
                      ? axis
                      : axis180
                    );
            }

            const tensor PerpA = tensor::I - sqr(rotationAxis_);
            const vector normalPerpA = normalised(PerpA & normal);
            const vector negNbrNormalPerpA = normalised(PerpA & negNbrNormal);
            const scalar theta =
                acos(min(max(normalPerpA & negNbrNormalPerpA, -1), 1));

            // Calculate the centre of rotation, if necessary
            if (transformType_ != oldTransformType)
            {
                const tensor R = quaternion(rotationAxis_, theta).R();
                tensor A = tensor::I - R;
                vector b = ctr - (R & nbrCtr);
                const label i = findMax(cmptMag(rotationAxis_));
                forAll(b, j)
                {
                    A(i, j) = j == i;
                }
                b[i] = 0;
                rotationCentre_ = inv(A) & b;
            }

            // Calculate the angle
            if (rotationalMethod_ == CENTREOFMASS)
            {
                const point ctrTang =
                    normalised(ctr - (ctr & rotationAxis_*rotationAxis_));

                const point nbrCtrTang =
                    normalised
                    (
                        nbrCtr - (nbrCtr & rotationAxis_*rotationAxis_)
                    );

                const scalar theta =
                    acos(min(max(ctrTang & nbrCtrTang, -1), 1));

                rotationAngle_ =
                  - sign((normalPerpA ^ negNbrNormalPerpA) & rotationAxis_)
                   *radToDeg(theta);
            }
            else
            {
                // Rotational method == "furthestPoint"
                rotationAngle_ =
                  - sign((normalPerpA ^ negNbrNormalPerpA) & rotationAxis_)
                   *radToDeg(theta);
            }
        }

        // If the transformation is known to be translational, then we just
        // need to set the separation
        if (transformType_ == TRANSLATIONAL)
        {
            separation_ = delta;
        }

        // Update the transform object
        transformComplete_ = true;
        update();

        // Print the results of calculation
        if (debug)
        {
            Info<< "Transformation calculated between patches " << name
                << " and " << nbrName << ":" << nl << token::BEGIN_BLOCK << nl
                << incrIndent;

            cyclicTransform::write(Info);

            Info<< decrIndent << token::END_BLOCK << nl << endl;
        }
    }

    // Check the transformation is correct to within the matching tolerance
    const point nbrCtrT =
        transform_.transformPosition(nbrCtr);

    const scalar ctrNbrCtrTDistance = mag(ctr - nbrCtrT);

    if (ctrNbrCtrTDistance > lengthScale*matchTolerance)
    {
        FatalErrorInFunction
            << "The distance between the centre of patch " << name
            << " and the transformed centre of patch " << nbrName << " is "
            << ctrNbrCtrTDistance << "."
            << nl
            << "This is greater than the match tolerance of "
            << lengthScale*matchTolerance << " for the patch."
            << nl
            << "Check that the patches are geometrically similar and that any "
            << "transformations defined between them are correct"
            << nl
            << "It might be possible to fix this problem by increasing the "
            << "\"matchTolerance\" setting for this patch in the boundary "
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicTransform::~cyclicTransform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicTransform::write(Ostream& os) const
{
    const label oldPrecision = os.precision();

    os.precision(16);

    if (transformType_ != UNSPECIFIED)
    {
        os.writeEntry("transformType", transformTypeNames[transformType_]);
    }

    if (transformType_ == ROTATIONAL)
    {
        os.writeEntry("rotationAxis", rotationAxis_);
        os.writeEntry("rotationCentre", rotationCentre_);

        os.writeEntry
        (
            "rotationalMethod",
            rotationalMethodNames[rotationalMethod_]
        );

        if (transformComplete_)
        {
            os.writeEntry("rotationAngle", rotationAngle_);
        }
    }

    if (transformType_ == TRANSLATIONAL)
    {
        if (transformComplete_)
        {
            os.writeEntry("separation", separation_);
        }
    }

    os.precision(oldPrecision);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::cyclicTransform Foam::operator&
(
    const transformer& t,
    const cyclicTransform& c0
)
{
    cyclicTransform c1(c0);

    if (c1.transformType_ == cyclicTransform::ROTATIONAL)
    {
        c1.rotationAxis_ = normalised(t.transform(c1.rotationAxis_));
        c1.rotationCentre_ = t.transformPosition(c1.rotationCentre_);
    }

    if (c1.transformType_ == cyclicTransform::TRANSLATIONAL)
    {
        if (c1.transformComplete_)
        {
            c1.separation_ = t.transform(c1.separation_);
        }
    }

    if (c1.transformComplete_)
    {
        c1.update();
    }

    return c1;
}


Foam::cyclicTransform Foam::inv(const cyclicTransform& c0)
{
    cyclicTransform c1(c0);

    if (c1.transformType_ == cyclicTransform::ROTATIONAL)
    {
        if (c1.transformComplete_)
        {
            c1.rotationAngle_ = - c1.rotationAngle_;
        }
    }

    if (c1.transformType_ == cyclicTransform::TRANSLATIONAL)
    {
        if (c1.transformComplete_)
        {
            c1.separation_ = - c1.separation_;
        }
    }

    if (c1.transformComplete_)
    {
        c1.update();
    }

    return c1;
}


// ************************************************************************* //
