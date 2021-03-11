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

def is_regular_sequence(poly_system):
    raise NotImplementedError("The following code is incorrect.")
    if len(poly_system) <= 0: return True
    if poly_system[0].is_unit() and len(poly_system) > 1: return False
    for i in range(1, len(poly_system)):
        if poly_system[i] in Ideal(poly_system[:i]): return False
    return True
