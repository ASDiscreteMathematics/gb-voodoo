load('analyze.sage')
load('poly_systems.sage')
load('dreg_macaulay_matrix.sage')
from sage.rings.polynomial.hilbert import first_hilbert_series

def max_spoly_deg(gb):
    d = 0
    for j in range(len(gb)):
        for i in range(j):
            s = s_poly(gb[i], gb[j])
            d = max(d, s.degree())
    return d

def max_lcm_deg(gb):
    d = 0
    for j in range(len(gb)):
        for i in range(j):
            lcm = gb[i].lm().lcm(gb[j].lm())
            d = max(d, lcm.degree())
    return d

set_verbose(-1)
magma_installed = False
f5 = F5()
sys_gen = PolynomialSystem(101)
eq_to_mx = vector([0]*9)
eq_to_sp = vector([0]*9)
eq_to_lc = vector([0]*9)
eq_to_f5 = vector([0]*9)
eq_to_mm = vector([0]*9)
eq_to_ma = vector([0]*9)
eq_to_hp = vector([0]*9)
eq_to_fh = vector([0]*9)
eq_to_hr = vector([0]*9)
for _ in range(100):
    polys = sys_gen.random(num_vars=3, degree=3, terms=6)
    gb0, voos = f5(polys, homogenize=False)
    if magma_installed:
        gb1, degs = magma.GroebnerBasis(polys, Faugere=True, nvals=2)
        gb1 = [g.sage() for g in gb1] # convert from MagmaElement to sage object
    else:
        gb1 = magma_free(f'r<{str(gb0[0].parent().gens())[1:-1]}>:=PolynomialRing(GF({len(gb0[0].base_ring())}),{len(gb0[0].parent().gens())},"grevlex");GroebnerBasis({polys}:Faugere:=true)')
        gb1 = [gb0[0].parent()(g.strip()) for g in gb1[2:-2].split(',')]
        degs = [0] # current (2021-03-18) magma version is V2.25-8, but we need V2.25-9 for GroebneBasis to return all intermediate working degrees.
    gb2 = Ideal(polys).groebner_basis()
    assert sorted(gb0) == sorted(gb1), f"Custom F5 differs from Magma"
    assert sorted(gb0) == sorted(gb2), f"Custom F5 differs from Singular"
    I = Ideal(gb0)
    dreg_mx = max(p.degree() for p in gb0)
    dreg_sp = max_spoly_deg(gb0)
    dreg_lc = max_lcm_deg(gb0)
    dreg_f5 = f5.dreg
    dreg_mm = d_reg_mm(polys)
    dreg_ma = max([0] + list(degs))
    dreg_hp = I.homogenize().hilbert_series().numerator().degree()
    dreg_fh = [deg for deg, coeff in enumerate(list(first_hilbert_series(I))) if coeff <= 0][0]
    if magma_installed:
        _, hil_reg = magma.HilbertPolynomial(I.homogenize(), nvals=2)
    else:
        hil_reg = magma_free(f'r<{str(gb0[0].parent().gens())[1:-1]}>:=PolynomialRing(GF({len(gb0[0].base_ring())}),{len(gb0[0].parent().gens())},"grevlex");HilbertPolynomial(Homogenization(Ideal({gb0})))')
        hil_reg = int(hil_reg.split()[-1])
    eq_to_mx += vector([dreg_mx == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_sp += vector([dreg_sp == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_lc += vector([dreg_lc == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_f5 += vector([dreg_f5 == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_mm += vector([dreg_mm == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_ma += vector([dreg_ma == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_hp += vector([dreg_hp == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_fh += vector([dreg_fh == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    eq_to_hr += vector([hil_reg == other for other in [dreg_mx, dreg_sp, dreg_lc, dreg_f5, dreg_mm, dreg_ma, dreg_hp, dreg_fh, hil_reg]])
    print(f"{dreg_mx} {dreg_sp} {dreg_lc} {dreg_f5} {dreg_mm} {dreg_ma} {dreg_hp} {dreg_fh} {hil_reg}", end="\r") # liveness indicator
print()
# maximum degree in GB, maximum degree of any S-poly, maximum degree of any two GB poly's lcm, dreg reported by custom F5, …
# … dreg uiing Macaulay Matruces, highest working degree of Magma, Hibert Poincare series, degree of first term with non-positive…
# … coefficient of hilbert series, Hilbert Regularity as reported by magma
print(f" max spl lcm  f5 mac mag hps fhs  hr")
print(matrix([eq_to_mx, eq_to_sp, eq_to_lc, eq_to_f5, eq_to_mm, eq_to_ma, eq_to_hp, eq_to_fh, eq_to_hr]))
