const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

bool MULESCorr(alphaControls.getOrDefault("MULESCorr", false));

bool turbCorr(alphaControls.getOrDefault("turbCorr", true));

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
bool alphaApplyPrevCorr
(
    alphaControls.getOrDefault("alphaApplyPrevCorr", false)
);
