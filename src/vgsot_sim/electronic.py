from .configs import PhysicalConstantsConfig


def electronic(V1, V2, V3, R_MTJ, constants: PhysicalConstantsConfig):
    """ 
    Input: 
    V1 = terminal voltage at MTJ
    V2 = terminal voltage at AFM strip
    V3 = right terminal voltage
    R_MTJ = MTJ resistance

    Output:
    I_SOT, V_MTJ as tuple
    """
    V_MTJ = R_MTJ * (4*V1 - 2*V2 - 2*V3)/(4*R_MTJ + constants.R_W)
    I_SOT = (V2-V3)/constants.R_W

    return (I_SOT, V_MTJ)
