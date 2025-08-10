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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/layerCellHandling/cellMetric.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::cellMetric::calculateTransformTensors()
{
    const vectorField& direction = stencilPtr_->anisotropicDirection();
    const scalarField& anisotropy = stencilPtr_->anisotropy();
    forAll(mesh_.cells(), cI)
    {
        vector deltaFace = direction[cI];

        if (anisotropy[cI] == 1)
        {
            transformationTensor_[cI] = tensor::identity();
        }
        else
        {
            deltaFace /= mag(deltaFace);

            vector xAxis(1, 0, 0);
            vector yAxis(0, 1, 0);
            vector zAxis(0, 0, 1);
            vector dir(vector::zero);

            scalar maxAx =
                max(mag(zAxis&deltaFace),max(mag(yAxis&deltaFace),mag(xAxis&deltaFace)));
            if (mag(zAxis&deltaFace)==maxAx)
            {
                dir = zAxis;
            }
            else if (mag(yAxis&deltaFace)==maxAx)
            {
                dir = yAxis;
            }
            else
            {
                dir = xAxis;
            }

            tensor R(tensor::identity());
            vector rotationAxis = dir^deltaFace;

            if ((mag(rotationAxis) > ROOTVSMALL))
            {
                rotationAxis /= mag(rotationAxis);

                scalar theta = Foam::acos(deltaFace&dir);
                tensor crossProductTensor
                (
                    0, -rotationAxis.z(), rotationAxis.y(),
                    rotationAxis.z(), 0, -rotationAxis.x(),
                    -rotationAxis.y(), rotationAxis.x(), 0
                );
                R = cos(theta)*tensor::identity() + sin(theta)*crossProductTensor +
                    (1-cos(theta))*rotationAxis*rotationAxis;
            }

            tensor metric(R&xAxis, R&yAxis, R&zAxis);

            scalar AR = anisotropy[cI];

            if (mag(zAxis&deltaFace)==maxAx)
            {
                metric.zz() *= AR;
                metric.zy() *= AR;
                metric.zx() *= AR;
            }
            else if (mag(yAxis&deltaFace)==maxAx)
            {
                metric.yz() *= AR;
                metric.yx() *= AR;
                metric.yy() *= AR;
            }
            else
            {
                metric.xx() *= AR;
                metric.xy() *= AR;
                metric.xz() *= AR;
            }

            transformationTensor_[cI] = R&metric;
        }
    }
}

void Foam::cellMetric::externalProductTransform()
{
    forAll(mesh_.cells(), cI)
    {
        const tensor& tT = transformationTensor_[cI];
        tensor& extT = extProductTransformation_[cI];

        extT = hinv(tT);
           extT = det(tT)*extT.T();
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::cellMetric::cellMetric(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    stencilPtr_
    (
        stencilBuilder::New(dict_, mesh_)
    ),
    transformationTensor_(mesh_.cells().size()),
    extProductTransformation_(mesh_.cells().size())
{
    calculateTransformTensors();
    externalProductTransform();
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //
Foam::cellMetric::~cellMetric(){}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::tensorField& Foam::cellMetric::tTensor() const
{
    return transformationTensor_;
}

const Foam::tensorField& Foam::cellMetric::extProductTensor() const
{
    return extProductTransformation_;
}

const Foam::scalarField& Foam::cellMetric::aspectRatio() const
{
    return stencilPtr_->anisotropy();
}

const Foam::autoPtr<Foam::stencilBuilder> Foam::cellMetric::stencilPtr() const
{
    return stencilPtr_;
}
