from itertools import combinations

def polynomial_division(f, divisors):
    f_original = f
    s = len(divisors)
    quotients = [0]*s
    remainder = 0
    while f != 0:
        i = 0
        division_occured = False
        while i < s and not division_occured:
            divisable = False
            try:
                divisable = divisors[i].lt().divides(f.lt())
            except NotImplementedError as e:
                pass # _beautiful_ solution
            if divisable:
                q, _ = f.lt().quo_rem(divisors[i].lt())
                quotients[i] += q
                f = f - q * divisors[i]
                division_occured = True
            else:
                i += 1
        if not division_occured:
            r = f.lt()
            remainder += r
            f -= r
    # assert f_original == sum([q * d for (q, d) in zip(quotients, divisors)]) + remainder
    return quotients, remainder

def s_poly(f, g):
    l = f.lm().lcm(g.lm())
    factor_f = l // f.lt()
    factor_g = l // g.lt()
    return factor_f * f - factor_g * g

def is_regular_sequence(poly_system, give_reason=False):
    if len(poly_system) <= 0:
        if give_reason: return True, "trivial system"
        return True
    if poly_system[0].is_unit() and len(poly_system) > 1:
        if give_reason: return False, "first poly spans ring"
        return False
    ring = poly_system[0].parent()
    for i in range(1, len(poly_system)):
        quo_ring = ring.quo(Ideal(poly_system[:i]))
        f = quo_ring(poly_system[i])
        if f == 0:
            if give_reason: return False, f"f_{i} is 0 mod f_0{', …' if i>2 else ''}{f', f_{i-1}' if i>1 else ''}"
            return False
        if magma.IsZeroDivisor(f).sage():
            if give_reason: return False, f"f_{i} is in R/<f_0{', …' if i>2 else ''}{f', f_{i-1}' if i>1 else ''}>"
            return False
    if give_reason: return True, "is regular sequence"
    return True

def is_regular_sequence_m2(poly_system, give_reason=False):
    macaulay2('loadPackage "Depth"')
    return macaulay2.isRegularSequence(poly_system).sage()

def hilbert_regularity(I):
    '''
    Compute the Hilbert regularity of R/I where R = I.ring() and I.dimension() <= 0.
    This is done by iterating through all n-tuples of the Gröbner basis' leading monomials,
    computing their lcm, then determining if that lcm is actually a corner of the staircase.
    The corner that is the furthest from the origin determines the Hilbert regularity.
    '''
    if I.dimension() > 0:
        raise NotImplementedError(f"Ideal must not be of positive dimension, but is {I.dimension()}.")
    gens = I.ring().gens() # all variables
    n = len(gens)
    xyz = reduce(operator.mul, gens, 1)
    gb_lm = [f.lm() for f in I.groebner_basis()]
    I_lm = Ideal(gb_lm)
    hil_reg = 0
    for lms in combinations(gb_lm, n):
        m = lcm(lms)
        # are we considering a meaningful combination of lm's?
        # i.e., does every variable make an appearance in m?
        if len(m.degrees()) != n or not all(m.degrees()):
            continue
        m = m / xyz # 1 step towards origin along all axes
        assert I.ring()(m) == m.numerator() # no negative exponents, please
        m = I.ring()(m)
        # are we in a corner of the staircase?
        # i.e., not in the ideal, but moving 1 step along any axis, we end up in the ideal?
        if not m in I_lm and all([v*m in I_lm for v in gens]):
            hil_reg = max(hil_reg, m.degree())
    return hil_reg