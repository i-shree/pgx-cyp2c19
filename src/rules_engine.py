# The "Red Flag" Letters (Alleles)
#     rs4244285 (CYP2C19*2): The "risk" letter here is A.
# If your result is AG or AA, the drug will have reduced effectiveness.
#     rs4986893 (CYP2C19*3): The "risk" letter here is A
# If your result is GA or AA, the drug will have reduced effectiveness.
#.   rs12248560 (CYP2C19*17): Rapid reaction.

def classify_metabolizer(rs4244285: int, rs12248560: int, rs4986893: int) -> str:
    # LOF variants (*2, *3)
    lof = rs4244285 + rs4986893
    # GOF variant (*17)
    gof = rs12248560

    if lof >= 2:
        return "Poor"
    if lof == 1:
        return "Intermediate"
    if gof >= 1:
        return "Rapid"
    return "Normal"


def recommendation(metabolizer: str) -> str:
    if metabolizer in {"Poor", "Intermediate"}:
        return "Consider alternative antiplatelet therapy."
    if metabolizer == "Rapid":
        return "Monitor for bleeding risk."
    return "Standard clopidogrel therapy appropriate."
