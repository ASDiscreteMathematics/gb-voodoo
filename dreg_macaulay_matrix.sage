load('common.sage')

def all_monoms_upto_deg(ring, d):
    all_monoms = set()
    last_monoms = [ring(1)]
    for i in range(d):
        all_monoms.update(last_monoms)
        last_monoms = [l*v for l in last_monoms for v in ring.gens()]
    all_monoms.update(last_monoms)
    return sorted(all_monoms)[::-1]

def macaulay_matrix(polys, d):
    ring = polys[0].parent()
    columns = all_monoms_upto_deg(ring, d)
    mm_pols = []
    for p in polys:
        factors = [f for f in columns if f.lm()*p.lm() in columns]
        factors = factors[::-1] # original polynomial in highest possible row
        mm_pols += [f*p for f in factors]
    mm_rows = [[p.monomial_coefficient(c) for c in columns] for p in mm_pols]
    return matrix(mm_rows)

def gauss_elimination(A):
    A = copy(A)
    nrows, ncols = A.nrows(), A.ncols()
    for c in range(ncols):
        for r in range(nrows):
            if A[r,c] and not any(A[r,i] for i in range(c)):
                a_inverse = ~A[r,c]
                A.rescale_row(r, a_inverse)
                for i in range(nrows):
                    if i != r and A[i,c]:
                        minus_b = -A[i,c]
                        A.add_multiple_of_row(i, r, minus_b)
                break
    empty_rows = [i for i in range(nrows) if not any(A.row(i))]
    A = A.delete_rows(empty_rows)
    A = matrix(sorted(A)[::-1])
    return A

def interreduce(polys):
    i = 0
    while i < len(polys):
        reductee = polys[i]
        polys_wo_reductee = [p for p in polys if p != reductee]
        _, rem = polynomial_division(reductee, polys_wo_reductee)
        if not rem:
            polys = polys_wo_reductee
            i = 0
        else:
            i += 1
    return polys

def buchberger_criterion(gb):
    for j in range(len(gb)):
        for i in range(j):
            s = s_poly(gb[i], gb[j])
            _, rem = polynomial_division(s, gb)
            if rem:
                return False
    return True

def d_reg_mm(polys):
    ring = polys[0].parent()
    d = 0
    is_gb = False
    while not is_gb:
        columns = vector(all_monoms_upto_deg(ring, d))
        mm = macaulay_matrix(polys, d)
        ge = gauss_elimination(mm)
        gb = [row * columns for row in ge]
        gb_red = interreduce(gb)
        quos_rems = [polynomial_division(p, gb_red) for p in polys]
        rems = [r for _, r in quos_rems]
        is_gb = buchberger_criterion(gb_red) and not any(rems)
        if not is_gb: d += 1
    assert sorted(gb_red) == sorted(Ideal(polys).groebner_basis()),\
        f"We did not compute the reduced Gröbner basis"
    return d
