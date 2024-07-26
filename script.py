
from typing import Tuple, Optional;
import py_ecc as ecc
import pdb 

from py_ecc.fields import (
    optimized_bls12_381_FQ as FQ,
    bn128_FQ2 as FQ2,
)


from py_ecc.bn128 import (
    Z1,
    Z2,
    b,
    b2,
    field_modulus as q,
    is_inf,
    is_on_curve
)

POW_2_381 = 2**381
POW_2_383 = 2**383

def get_flags(z: int) -> Tuple[bool, bool, bool]:
    c_flag = bool((z >> 383) & 1)  # The most significant bit.
    b_flag = bool((z >> 382) & 1)  # The second-most significant bit.
    a_flag = bool((z >> 381) & 1)  # The third-most significant bit.
    return c_flag, b_flag, a_flag

def compress():
    x = [20456271874651344948253341295392711905477135017044966153115686896847760857935, 3750429983862860823742319760054991864945596549379822962956770646884312979236]
    y = [703569469678413069677292343660238496680944522722992077714950198346770215535, 8744331013967169952726744659258818632296620395222738052973422057164413318811]

    x_re , x_im = x
    y_re, y_im = y

    if y_im > 0 :
        a_flag1 = (y_im * 2) // q 
    else: 
        a_flag1 = (y_re * 2) // q


    # Imaginary part of x goes to z1, real part goes to z2
    # c_flag1 = 1, b_flag1 = 0
    z1 = x_im + a_flag1 * POW_2_381 + POW_2_383
    # a_flag2 = b_flag2 = c_flag2 = 0
    z2 = x_re
    return (z1, z2)

out = compress()

def is_point_at_infinity(z1: int, z2: Optional[int] = None) -> bool:
    """
    If z2 is None, the given z1 is a G1 point.
    Else, (z1, z2) is a G2 point.
    """
    return (z1 % POW_2_381 == 0) and (z2 is None or z2 == 0)

FQ2_ORDER = q**2 - 1
EIGHTH_ROOTS_OF_UNITY = tuple(FQ2([1, 1]) ** ((FQ2_ORDER * k) // 8) for k in range(8))


#
# G2
#
def modular_squareroot_in_FQ2(value: FQ2) -> Optional[FQ2]:
    """
    Given value=``x``, returns the value ``y`` such that ``y**2 % q == x``,
    and None if this is not possible. In cases where there are two solutions,
    the value with higher imaginary component is favored;
    if both solutions have equal imaginary component the value with higher real
    component is favored.
    """
    candidate_squareroot = value ** ((FQ2_ORDER + 8) // 16)
    check = candidate_squareroot**2 / value
    if check in EIGHTH_ROOTS_OF_UNITY[::2]:
        x1 = (
            candidate_squareroot
            / EIGHTH_ROOTS_OF_UNITY[EIGHTH_ROOTS_OF_UNITY.index(check) // 2]
        )
        x2 = -x1
        x1_re, x1_im = x1.coeffs
        x2_re, x2_im = x2.coeffs
        return x1 if (x1_im > x2_im or (x1_im == x2_im and x1_re > x2_re)) else x2
    return None


def decompress(p):
    z1, z2 = p
    c_flag1, b_flag1, a_flag1 = get_flags(z1)

    # c_flag == 1 indicates the compressed form
    # MSB should be 1
    if not c_flag1:
        raise ValueError("c_flag should be 1")

    is_inf_pt = is_point_at_infinity(z1, z2)

    if b_flag1 != is_inf_pt:
        raise ValueError(f"b_flag should be {int(is_inf_pt)}")

    if is_inf_pt:
        # 3 MSBs should be 110
        if a_flag1:
            raise ValueError("a point at infinity should have a_flag == 0")
        return Z2

    # Else, not point at infinity
    # 3 MSBs should be 100 or 101
    x1 = z1 % POW_2_381
    # Ensure that x1 is less than the field modulus.
    if x1 >= q:
        raise ValueError(f"x1 value should be less than field modulus. Got {x1}")

    # Ensure that z2 is less than the field modulus.
    if z2 >= q:
        raise ValueError(f"z2 point value should be less than field modulus. Got {z2}")

    x2 = z2
    # x1 is the imaginary part, x2 is the real part
    x = FQ2([x2, x1])
    y = modular_squareroot_in_FQ2(x**3 + b2)
    if y is None:
        raise ValueError("Failed to find a modular squareroot")

    # Choose the y whose leftmost bit of the imaginary part is equal to the a_flag1
    # If y_im happens to be zero, then use the bit of y_re
    y_re, y_im = y.coeffs
    if (y_im > 0 and (y_im * 2) // q != int(a_flag1)) or (
        y_im == 0 and (y_re * 2) // q != int(a_flag1)
    ):
        y = FQ2((y * -1).coeffs)

    if not is_on_curve((x, y, FQ2([1, 0])), b2):
        raise ValueError("The given point is not on the twisted curve over FQ**2")
    return (x, y, FQ2([1, 0]))

print(out)
print(decompress(out))
# pdb.set_trace()

