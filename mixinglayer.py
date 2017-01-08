""" Mixing layer Model - K_L parameters & Reynolds (Shi et al., 2011, p.1220)"""

""" Global parameters"""
V_k = 1.004e-6 * 10 ** 6 * 60  # (mm2/min) # kinematic water viscosity at 20 deg. Celsius = 1.004 x 10-6 (m2/s)
D_w = 5.09677e-06 * 60 * 10**2  # mm2/min solute diff usivity in water -> Metolachlor = 5.0967719112e-006 (cm2 / s)
Sc = V_k / (D_w / 60)  # (-) Schmidt number


def reynolds(runoffvelocity, diameter, v_k=V_k):
    """ Overall Reynolds Number """
    re = (runoffvelocity * diameter) / v_k  # [-]
    return re  # (dimensionless)
    # intensityM = rainfall intensity (P: mm/min)
    # diameter = Length (L: mm)


def kfilm(diameter, runoffvelocity):
    kl = 0.664 * (D_w/diameter)*reynolds(runoffvelocity, diameter)**(float(1)/2)*Sc**(float(1)/3)  # mm/min
    return kl  # mm/min

"""
Kinematic Viscosity (to calculate reynolds overall number (Re_L) -
Correct, based on Wallach et al., 1988.
v_k = 1.004 x 10-6 (m2/s)
D_w = Solute diffusivity / Diffusion coefficient in water (cm2 / s)

Chemical parameter source:
http://www.gsi-net.com/en/publications/gsi-chemical-database/single/377.html
"""

"""
A second interpretation of based on Shi et al. (2011) uses Dynamic Viscosity instead (I think incorrect):

Viscocity of Water at 20 deg. Celcius - Dynamic Viscosity (to calculate Reynolds number)
- mu = 1.002 x 10-3 (Pa s or N s/m2)
- mu = 1.002 x 10-2 (g per cm second) -> used here#
0.01 poise = 0.01 g per cm second = 0.001 Pascal second = 0.001 N s/m2

mu = 1.002e-2 # dynamic water viscocity at 20 deg. Celcius = 1.002 x 10-2 (g per cm second)

# Water properties at 20 deg. Celcius
p_w = 0.99823 # water density at 20 deg. Celcius = 0.99823 g/cm^3

#
Re_1 = (ksat/60)*(p_w*diameter)/mu # (-), not correct, based on Shi et al., 2011; ksat in cm/s : 0.0167* 1min/60s
"""
