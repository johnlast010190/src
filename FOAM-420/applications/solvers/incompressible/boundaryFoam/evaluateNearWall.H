/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

    // Evaluate near-wall behaviour
    // parallel and 2D capable
{
    scalar vsum(0);
    scalar Rey(0);
    scalar uTau(0);
    scalar nutPlus(0);
    scalar yPlus(0);
    scalar uPlus(0);
    scalar kPlus(0);
    scalar epsilonPlus(0);

    tmp<volScalarField> nu(turbulence->nu());
    tmp<volScalarField> nut(turbulence->nut());
    tmp<volSymmTensorField> R(turbulence->devReff());
    tmp<volScalarField> epsilon(turbulence->epsilon());
    tmp<volScalarField> k(turbulence->k());

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const labelList& faceCells(currPatch.faceCells());
            const vectorField wallNormal(-1*currPatch.nf());

            forAll(faceCells, fci)
            {
                label fcellI(faceCells[fci]);
                scalar Vc(mesh.V()[fcellI]);
                vsum += Vc;

                scalar magUp
                (
                    mag(U[fcellI] - U.boundaryField()[patchi][fci])
                );

                scalar tauwc = (flowDirection
                     & R->boundaryField()[patchi][fci]
                     & wallNormal[fci]);

                scalar uTauc = ::sqrt(mag(tauwc));
                uTau += Vc*uTauc;

                yPlus += Vc*uTauc*y[fcellI]
                       /(nu()[fcellI] + ROOTVSMALL);

                uPlus += Vc*magUp/(uTauc + ROOTVSMALL);

                nutPlus += Vc*nut()[fcellI]/nu()[fcellI];

                kPlus += Vc*k()[fcellI]/(sqr(uTauc) + ROOTVSMALL);

                epsilonPlus += Vc*epsilon()[fcellI]*nu()[fcellI]
                             /(pow4(uTauc) + ROOTVSMALL);

                Rey += Vc*magUp*y[fcellI]/nu()[fcellI];
            }
        }
    }

    reduceToMaster(
        std::tie(vsum, uTau, yPlus, uPlus, nutPlus, kPlus, epsilonPlus, Rey),
        UniformParallelOp<sumOp<scalar>, 8>{}
    );

    uTau /= vsum;
    yPlus /= vsum;
    uPlus /= vsum;
    nutPlus /= vsum;
    kPlus /= vsum;
    epsilonPlus /= vsum;
    Rey /= vsum;

    Info<< "Rey = " << Rey
        << ", uTau = " << uTau
        << ", nut+ = " << nutPlus
        << ", y+ = " << yPlus
        << ", u+ = " << uPlus
        << ", k+ = " << kPlus
        << ", epsilon+ = " << epsilonPlus
        << endl;
}
