for (label patchi : snowPatchList)
{
    const tmp<volSymmTensorField> &Reff = turbulence->devRhoReff();

    const vectorField &Sfp = mesh.Sf().boundaryField()[patchi];       // 面积矢量
    const scalarField &magSfp = mesh.magSf().boundaryField()[patchi]; // 面积矢量模长
    const scalarField &rhop = rho.boundaryField()[patchi];

    const symmTensorField &Reffp = Reff().boundaryField()[patchi];
    
    vectorField &ssp = myWallShearStress.boundaryFieldRef()[patchi];
    ssp = (-Sfp / magSfp) & Reffp; // 剪切应力

    ssp /= rhop;

    forAll(rhop, i)
    {
        Info << "rho: "<< rhop[i] << endl;
    }
}
