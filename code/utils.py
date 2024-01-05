import numpy as  np
import astropy.units as u
from astropy.constants import k_B, h, c, m_p, m_e, eps0, e
k_B = k_B.to(u.erg / u.K)
h = h.to(u.erg * u.s)
c = c.to(u.cm / u.s)
m_p = m_p.to(u.g)
m_e = m_e.to(u.g)
eps0 = eps0.to(u.F / u.cm)
e.to(u.erg / u.V)

# SAHA

C = (2 * np.pi * m_e * k_B / (h**2))**(3/2)
C = C.to(u.K**(-3/2) * u.cm**(-3))

def saha(N_e, T, U_p, U_m, chi_lu):
    exp_fact = - (chi_lu / (k_B * T)).to(u.dimensionless_unscaled)
    exp = np.exp(exp_fact)
    fact1 = 2 * U_p / (N_e * U_m)
    result = fact1 * C * T**(3/2) * exp
    return result.to(u.dimensionless_unscaled)

def matrix_saha(a,b):
    arr = [[a, -1, 0],
           [-1, 0, b],
           [0, 1, -1]]
    return np.array(arr)

# BOLTZMANN

def boltzmann(gl, gu, chi_l, chi_u, T):
    rat_g = gu/gl
    exp_fact = (-(chi_u-chi_l)/(k_B*T)).to(u.dimensionless_unscaled)
    exp = np.exp(exp_fact)
    result = rat_g * exp
    return result.to(u.dimensionless_unscaled)

def matrix_boltz(a,b):
    arr = [[1, 1, 1],
           [a, -1, 0],
           [0, b, -1]]
    return np.array(arr)

# CROSS SECTIONS

# Rydberg constant
R = 1.0968e5 * u.cm**(-1)

# Electron scattering 

def sigma_es():
    fact1 = 8 * np.pi / 3
    fact2 = (e.esu.value)**2 / (m_e.value * c.value**2)
    result = fact1 * fact2**2
    return result 

def kappa_es(N_e):
    return sigma_es() * N_e
    

# Free-Free for H

def HI_sigma_ff(lamb, T):

    cte1 = 3.7e8
    Z = 1
    nu = c.value / lamb

    fact1 = cte1 * Z**2 
    fact2 = T**(1/2) * nu**3

    def gaunt_ff():

        f1 = 0.3456
        f2 = lamb * k_B.value * T / (h.value * c.value)
        f3 = 1/2
        f4 = (lamb*R.value)**(-1/3)

        return 1 + f1*f4*(f2+f3)
    
    result = gaunt_ff() * fact1 / fact2

    return result

def HI_kappa_ff(lamb, T, NHII, Ne):
    return NHII * Ne * HI_sigma_ff(lamb, T)

# Bound Free for H

def HI_sigma_bf(lamb, n=1):
    
    lamb_0 = (R.value/n**2)**(-1)


    lamb = lamb
    nu = c.value/lamb

    cte1 = 2.815e29 
    cte2 = 0.3456
    Z=1

    g_bf = 1 - cte2 * (lamb*R.value)**(-1/3) * (lamb*R.value / (n**2) - 1/2)
    sigma_bf = cte1 * Z**4 * g_bf / (n**5 * nu**3)
    sigma_bf = sigma_bf

    # use np.wwhere to set values below lambda_0 to zero
    sigma_bf = np.where(lamb > lamb_0, 0, sigma_bf)

    return sigma_bf

def HI_kappa_bf(lamb, NHI, n=1):
    return NHI * HI_sigma_bf(lamb, n)


# Bound-Free for H-

a_coefs = [1.99654,
           -1.18267e-5,
           2.64243e-6,
           -4.40524e-10,
           3.23992e-14,
           -1.39568e-18,
           2.78701e-23]

def Hm_sigma_bf(lamb):
    # Convert to Armstrong (lamb is in cm)
    lamb = lamb  * 10**8

    # Evaluating the polynomial using numpy.polyval
    result = np.polyval(a_coefs[::-1], lamb)

    return result * 10**-18

def Hm_kappa_bf(lamb, NHm):
    return NHm * Hm_sigma_bf(lamb)


# Free-Free for H-

# Coefficients of the f_n polynomials
f0cs = [-2.2763, -1.6850, 0.76661, -0.053346]
f1cs = [15.2827, -9.2846, 1.99381, -0.142632]
f2cs = [-197.789, 190.266, -67.9775, 10.6913, -0.625151]
fcs = [f0cs, f1cs, f2cs]

def Hm_sigma_ff(lamb, T):
    # Convert to Armstrong (lamb is in cm)
    lamb = lamb * 10**8

    # Using numpy.polyval for polynomial evaluation
    fs = [np.polyval(fc[::-1], np.log10(lamb)) for fc in fcs]

    theta = 5040 / T
    exponent = -26 + np.polyval(fs[::-1], np.log10(theta))

    result = 10**exponent
    return result


def Hm_kappa_ff(lamb, T, Pe, NHI):
    return Pe * Hm_sigma_ff(lamb, T) * NHI










# FUNCTION FOR FORMATTING

def sl(value, precision=2):
    """
    Convert a number into scientific notation with LaTeX formatting.

    :param value: The number to convert.
    :param precision: The number of decimal places to include. Default is 2.
    :return: A string representing the number in scientific notation.
    """
    # Convert the number to scientific notation
    format_string = f"{{:.{precision}e}}"
    sci_notation = format_string.format(value)

    # Split the string into the coefficient and exponent
    coefficient, exponent = sci_notation.split('e')
    
    # Remove leading '+' and unnecessary zeros from the exponent

    # Format into LaTeX notation
    if exponent != '+00':
        exponent = exponent.replace('+', '').lstrip('0')

        # Handle the case where exponent is negative and zero
        if exponent.startswith('-0'):
            exponent = '-' + exponent[2:]

        latex_string = f"${coefficient} \\times 10^{{{exponent}}}$"
    else:
        latex_string = f"${coefficient}$"

    return latex_string