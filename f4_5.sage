# -*- coding: utf-8 -*-

"""
Jean-Charles Faugère's F5 Algorithm.

These implementations are heavily inspired by John Perry's pseudocode
and Singular implementation of these algorithms.

See http://www.math.usm.edu/perry/Research/ for details.

The docstrings are almost verbatim copies from Just Gash's
explanations for each F5 function in his thesis: "On Efficient
Computation of Gröbner Bases". Note that Justin begins at f_m while
we begin at f_0, e.g. the first GB we calculate is <f_0> while Justin
calculates <f_m> first.

AUTHOR:
    -- 20081013 Martin Albrecht (initial version based on John Perry's pseudocode)
    -- 20081013 John Perry (loop from 0 to m-1 instead of m-1 to 0)
    -- 20090112 Martin Albrecht (F5SansRewriting)
    -- 20090124 Martin Albrecht and John Perry (F4F5)
    -- 20090126 John Perry (correction to compute_spols)
    -- 20210201 Jan Ferdinand Sauer (keep vectors of origin and syzygies)

EXAMPLE:
    sage: load('f4_5.py')
    sage: R.<x,y,z> = PolynomialRing(GF(29))
    sage: I =  R* [3*x^4*y + 18*x*y^4 + 4*x^3*y*z + 20*x*y^3*z + 3*x^2*z^3, \
                   3*x^3*y^2 + 7*x^2*y^3 + 24*y^2*z^3, \
                   12*x*y^4 + 17*x^4*z + 27*y^4*z + 11*x^3*z^2]
    sage: J = I.homogenize()

    sage: f5 = F5() # original F5
    sage: gb = f5(J)
    Increment 1
    Processing 1 pairs of degree 7
    Processing 1 pairs of degree 9
    Processing 1 pairs of degree 11
    Increment 2
    Processing 1 pairs of degree 6
    Processing 2 pairs of degree 7
    Processing 5 pairs of degree 8
    Processing 7 pairs of degree 9
    Processing 11 pairs of degree 10
    Processing 10 pairs of degree 11
    verbose 0 (...: f5.py, top_reduction) Reduction of 29 to zero.
    verbose 0 (...: f5.py, top_reduction) Reduction of 27 to zero.
    verbose 0 (...: f5.py, top_reduction) Reduction of 25 to zero.
    Processing 3 pairs of degree 12
    Processing 5 pairs of degree 13
    Processing 2 pairs of degree 14
    Processing 1 pairs of degree 16
    sage: f5.zero_reductions, len(gb)
    (3, 18)
    sage: Ideal(gb).basis_is_groebner()
    True

    sage: f5 = F5R() # F5 with interreduced B
    sage: gb = f5(J)
    Increment 1
    Processing 1 pairs of degree 7
    Processing 1 pairs of degree 9
    Processing 1 pairs of degree 11
    Increment 2
    Processing 1 pairs of degree 6
    Processing 2 pairs of degree 7
    Processing 5 pairs of degree 8
    Processing 7 pairs of degree 9
    Processing 11 pairs of degree 10
    Processing 10 pairs of degree 11
    verbose 0 (...: f5.py, top_reduction) Reduction of 29 to zero.
    verbose 0 (...: f5.py, top_reduction) Reduction of 27 to zero.
    verbose 0 (...: f5.py, top_reduction) Reduction of 25 to zero.
    Processing 3 pairs of degree 12
    Processing 5 pairs of degree 13
    Processing 2 pairs of degree 14
    Processing 1 pairs of degree 16
    sage: f5.zero_reductions, len(gb)
    (3, 18)
    sage: Ideal(gb).basis_is_groebner()
    True

    sage: f5 = F5C() # F5 with interreduced B and Gprev
    sage: gb = f5(J)
    Increment 1
    Processing 1 pairs of degree 7
    Processing 1 pairs of degree 9
    Processing 1 pairs of degree 11
    Increment 2
    Processing 1 pairs of degree 6
    Processing 2 pairs of degree 7
    Processing 5 pairs of degree 8
    Processing 7 pairs of degree 9
    Processing 11 pairs of degree 10
    Processing 10 pairs of degree 11
    verbose 0 (...: f5.py, top_reduction) Reduction of 29 to zero.
    verbose 0 (...: f5.py, top_reduction) Reduction of 27 to zero.
    verbose 0 (...: f5.py, top_reduction) Reduction of 25 to zero.
    Processing 3 pairs of degree 12
    Processing 4 pairs of degree 13
    Processing 1 pairs of degree 14

    sage: f5.zero_reductions, len(gb)
    (3, 18)
    sage: Ideal(gb).basis_is_groebner()
    True

    sage: f5 = F4F5() # F5-style F5
    sage: gb = f5(J)
    Increment 1
    Processing 1 pairs of degree 7
       5 x   13,    5,    0
    Processing 1 pairs of degree 9
      14 x   29,   14,    0
    Processing 1 pairs of degree 11
    Increment 2
    Processing 1 pairs of degree 6
       6 x   18,    6,    0
    Processing 2 pairs of degree 7
      11 x   23,   11,    0
    Processing 5 pairs of degree 8
      18 x   27,   18,    0
    Processing 7 pairs of degree 9
      19 x   23,   19,    0
    Processing 11 pairs of degree 10
      15 x   15,   15,    0
    Processing 10 pairs of degree 11
      14 x   11,   11,    3
    Processing 3 pairs of degree 12
    Processing 4 pairs of degree 13
    Processing 1 pairs of degree 14
    sage: f5.zero_reductions, len(gb)
    (3, 18)
    sage: Ideal(gb).basis_is_groebner()
    True

NOTE:
    For additional diagnostics there are a number of commented
    commands.  To count the number of reductions, one can uncomment
    commands to "interreduce" and comment out commands with
    "reduced_basis"; also uncomment commands with "normal_form" and
    comment out commands with "reduce".
"""

from functools import cmp_to_key

divides = lambda x,y: x.parent().monomial_divides(x,y)
LCM = lambda f,g: f.parent().monomial_lcm(f,g)

def compare_by_degree(f,g):
    if f.total_degree() > g.total_degree():
        return 1
    elif f.total_degree() < g.total_degree():
        return -1
    else:
        return -1 if f < g else ( 1 if f > g else 0 )

def unit_vec(ring, i, length):
    assert i < length
    vec = zero_vector(ring, length)
    vec[i] = 1
    return vec

class F5:
    """
    Jean-Charles Faugère's F5 Algorithm.
    """
    def __init__(self, F=None):
        if F:
            self.F = F
            self.R = F[0].parent()
            self.Rules = [[]]
            self.L = [0]
            self.syzygies = []
            self.zero_reductions = 0
            self.reductions = 0

    def phi(self, v):
        """
        Maps vector of origin v or voo at index 'i' to its polynomial
        under the input system. Always results in the polynomial at 'i'.
        Retuns the polynomial, allowing to check consistency of voos.
        """
        if isinstance(v, (int, sage.rings.integer.Integer)):
            v = self.voo(v)
        return v * vector(self.R, self.F)

    def voo(self, i):
        return self.L[i][2]

    def poly(self, i):
        return self.L[i][1]

    def sig(self, i):
        return self.L[i][0]

    def __call__(self, F, homogenize=True):
        if isinstance(F, sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal):
            F = F.gens()
        if homogenize and not all(f.is_homogeneous() for f in F):
            F = Ideal(F).homogenize()
            F = F.gens()
        return self.basis(F)

    def basis(self, F):
        """
        F5's main routine. Computes a Gröbner basis for F.

        INPUT:
            F -- a list of polynomials

        OUTPUT:
            G -- a list of polynomials; a Gröbner basis for <F>
        """
        phi = self.phi
        sig = self.sig
        voo = self.voo
        poly = self.poly
        incremental_basis = self.incremental_basis
        interreduce_basis = self.interreduce_basis
        regular_s_interreduce_basis = self.regular_s_interreduce_basis

        if not F: return [], []
        self.__init__(F)

        Rules = self.Rules
        L = self.L
        R = self.R
        m = len(F)

        f = F[0]
        L[0] = (Signature(1, 0), f/f.lc(), unit_vec(R, 0, m)/f.lc())
        assert phi(0) == poly(0), f" [!] Something wrong:\n{voo(0)}\n{poly(0)}."
        Rules.append([])

        Gprev = set([0])
        B = [poly(0)]
        B_voo = [voo(0)]

        for i in range(1,m):
            if get_verbose() >= 0: print(f"Starting incremental basis up to {i}")
            f = F[i]
            L.append( (Signature(1,i), f/f.lc(), unit_vec(R, i, m)/f.lc()) )
            Gcurr = incremental_basis(i, B, B_voo, Gprev)
            for j in range(len(Gcurr)):
                if poly(j) == 1:
                    return [poly(j)], [voo(j)]
            Gcurr = regular_s_interreduce_basis(Gcurr)
            Gprev = Gcurr
            B = [poly(l) for l in Gcurr]
            B_voo = [voo(l) for l in Gcurr]
        Gcurr = interreduce_basis(Gcurr)
        B = [poly(l) for l in Gcurr]
        B_voo = [voo(l) for l in Gcurr]
        return B, B_voo

    def incremental_basis(self, i, B, B_voo, Gprev):
        """
        adapted from Justin Gash (p.49):

        'This is the portion of the algorithm that is called (m-1)
         times and at the end of each call a new Gröbner basis is
         produced. After the first call to this algorithm, the
         Gröbner basis for <f_0 , f_1> is returned. In general, after
         the k-th call to this algorithm, the Gröbner basis for
         <f_0,...,f_k>. This is why F5 is called an iterative
         algorithm.  The process used by this algorithm is similar to
         many Gröbner basis algorithms: it moves degree-by-degree; it
         generates a new set of S- polynomials S to consider; it
         reduces this new set S of S-polynomials by G_curr and
         G_prev.'
        """
        L = self.L
        critical_pair = self.critical_pair
        compute_spols = self.compute_spols
        reduction = self.reduction
        Rules = self.Rules

        curr_idx = len(L) - 1
        Gcurr = Gprev.union([curr_idx])
        Rules.append( [] )

        P = reduce(lambda x,y: x.union(y), [critical_pair(curr_idx, j, i, Gprev) for j in Gprev], set())
        while P:
            d = min(t.degree() for (t,k,u,l,v) in P)
            Pd = [(t,k,u,l,v) for (t,k,u,l,v) in P if t.degree() == d]
            if get_verbose() >= 0: print(f"Processing {len(Pd)} pair{'' if len(Pd)==1 else 's'} of degree {d}")
            if get_verbose() >= 2: [print(each) for each in Pd]
            P = P.difference(Pd)
            S = compute_spols(Pd)
            R = reduction(S, B, B_voo, Gprev, Gcurr)
            for k in R:
                P = reduce(lambda x,y: x.union(y), [critical_pair(j, k, i, Gprev) for j in Gcurr], P)
                Gcurr.add(k)
        if get_verbose() >= 2: print(f"Ended with {len(Gcurr)} polynomials")
        return Gcurr

    def regular_s_interreduce_basis(self, Gcurr):
        voo, sig, poly = self.voo, self.sig, self.poly
        regular_s_reduce = self.regular_s_reduce
        get_sig_from_voo = self.get_sig_from_voo
        phi = self.phi
        L = self.L

        G_red = set()
        if get_verbose() >= 1: print(f"Regular-s-interreducing {len(Gcurr)} polynomials.")
        for k in Gcurr:
            if not poly(k): continue
            if sig(k) in [sig(j) for j in G_red]: continue
            assert sig(k) == get_sig_from_voo(voo(k)), "Mismatching sig and voo."
            v = regular_s_reduce(k, [j for j in Gcurr if j != k])
            p = phi(v)
            if not p: continue
            v, p = v/p.lc(), p/p.lc()
            s = get_sig_from_voo(v)
            assert s == sig(k), f"Signature has changed during regular-s-reduction."
            if v == voo(k):
                G_red.add(k)
            else:
                assert s == get_sig_from_voo(v), "Mismatching sig and voo."
                L.append( (s, p, v) )
                G_red.add(len(L) - 1)
        return G_red

    def interreduce_basis(self, Gcurr):
        """
        Given a list of indices comprising a preliminary Gröbner basis,
        returns a list of indices of an interreduced basis spanning the same ideal.
        """
        voo, sig, poly = self.voo, self.sig, self.poly
        get_sig_from_voo = self.get_sig_from_voo
        voo_reduce = self.voo_reduce
        phi = self.phi
        L = self.L

        G_red = set()
        for k in Gcurr:
            if not poly(k): continue
            B = [poly(j) for j in Gcurr if j != k]
            Bv = [voo(j) for j in Gcurr if j != k]
            p, v = voo_reduce(k, B, Bv)
            assert poly(k).reduce(B) == p, f" [!] Buggy behavior in 'voo_reduce'!"
            if not p: continue
            if v == voo(k):
                G_red.add(k)
            else:
                s = get_sig_from_voo(v)
                L.append( (s, p, v) )
                G_red.add(len(L) - 1)
        return G_red

    def regular_s_reduce(self, k, G):
        """
        Given index k and list of indices G, returns true iff polynomial k is
        reducible by any element in G without changing signature of k.
        """
        voo, sig, poly = self.voo, self.sig, self.poly
        get_sig_from_voo = self.get_sig_from_voo
        phi = self.phi
        R = self.R

        v = voo(k)[:] # copy, don't modify
        r = 0
        while phi(v) != r:
            m = (phi(v) - r).lt()
            i = 0
            reduction_occured = False
            while i < len(G) and not reduction_occured:
                assert phi(voo(G[i])) == poly(G[i]), f"Mismatching poly and voo."
                assert sig(G[i]) == get_sig_from_voo(voo(G[i])), f"Mismatching sig and voo."
                h_i_lt = poly(G[i]).lt()
                if h_i_lt.divides(m):
                    sig_d = R(m / h_i_lt) * sig(G[i])
                    sig_d = Signature(sig_d[0]/sig_d[0].lc(), sig_d[1])
                    sig_v = get_sig_from_voo(v)
                    if sig_d < sig_v:
                        v -= R(m / h_i_lt) * voo(G[i])
                        reduction_occured = True
                        assert sig_v == get_sig_from_voo(v), f"Signature has changed!"
                    else:
                        i += 1
                else:
                    i += 1
            if not reduction_occured:
                r += m
        return v

    def get_sig_from_voo(self, voo):
        '''
        Given vector of origin, return its pot-signature.
        '''
        if voo.is_zero(): return Signature(0, 0)
        i = next((j for j, x in enumerate(voo[::-1]) if x), None)
        i = len(voo) - i - 1 # reverse above reversing: need rightmost non-zero index
        mon = max(voo[i].monomials())
        return Signature(mon, i)

    def critical_pair(self, k, l, i, Gprev):
        """
        adapted from Justin Gash (p.51):

        'It is the subroutine critical_pair that is responsible for
         imposing the F5 Criterion from Theorem 3.3.1. Note that in
         condition (3) of Theorem 3.3.1, it is required that all pairs
         (r_i, r_j) be normalized. The reader will recall from
         Definition 3.2.2 that a pair is normalized if:

         (1) S(k) = m_0*F_{e_0} is not top-reducible by <f_0, ..., f_{e_0}-1>

         (2) S(l) = m_1*F_{e_1} is not top-reducible by <f_0, ..., f_{e_1}-1>

         (3) S(m_0*k) > S(m_1*l)

         If these three conditions are not met in critical_pair (note
         that the third condition will always be met because
         cirtical_pair forces it to be met), the nominated critical
         pair is dropped and () is returned.

         Once we have collected the nominated critical pairs that pass
         the F5 criterion test of critical_pair, we send them to
         compute_spols.'
        """
        poly = self.poly
        sig = self.sig
        is_top_reducible = self.is_top_reducible
        is_rewritable = self.is_rewritable

        #print(f"crit_pair({k}, {l}, {i}, {Gprev})")
        #print(f"{self.L}")
        tk = poly(k).lt()
        tl = poly(l).lt()
        t = LCM(tk, tl)
        u0 = t//tk
        u1 = t//tl
        m0, e0 = sig(k)
        m1, e1 = sig(l)
        if e0 == e1 and u0*m0 == u1*m1:
            return set()

        # Stegers and Gash leave out the == i check, Faugere and Perry
        # have it. It is unclear for now, whether the check is
        # necessary.
        if e0 == i and is_top_reducible(u0*m0, Gprev):
            return set()
        if e1 == i and is_top_reducible(u1*m1, Gprev):
            return set()
        # This check was introduced by Stegers, it isn't strictly
        # necessary, see class F5SansRewriting below
        if is_rewritable(u0, k) or is_rewritable(u1, l):
            return set()
        if u0/u0.lc() * sig(k) < u1/u1.lc() * sig(l):
            u0, u1 = u1, u0
            k, l = l, k
        return set([(t,k,u0,l,u1)])

    def compute_spols(self, P):
        """
        adapted from Justin Gash (p.51):

        'Though at first glance this subroutine may look complicated,
         compute_spols essentially does one thing: form the new
         S-polynomials output from critical_pairs as admissible signed
         polynomials. We note that, because critical_pairs ensured
         that S(u*k) < S(v*l), we know that the signature of all new
         polynomials will always be of the form u_L*S(r_{i_L}) in
         compute_spols.'
        """
        voo, sig, poly = self.voo, self.sig, self.poly
        phi = self.phi
        get_sig_from_voo = self.get_sig_from_voo
        is_rewritable = self.is_rewritable
        syzygies = self.syzygies
        add_rule = self.add_rule

        L = self.L

        S = list()
        P = sorted(P, key=lambda x: x[0])
        for (t,k,u,l,v) in P:
            if not is_rewritable(u,k) and not is_rewritable(v,l):
                assert sig(k) == get_sig_from_voo(voo(k)), f"Mismatching sig and voo: index {k}."
                assert sig(l) == get_sig_from_voo(voo(l)), f"Mismatching sig and voo: index {l}."
                s = u*poly(k)-v*poly(l) # S-Polynomial
                s_voo = u*voo(k)-v*voo(l)
                if s != 0:
                    s_voo /= s.lc() # normalize
                    s /= s.lc()
                sig_k = u/u.lc() * sig(k)
                assert phi(s_voo) == s, "Mismatching voo and poly."
                assert sig_k == get_sig_from_voo(s_voo), "Mismatching sig and voo."
                L.append( (sig_k, s, s_voo) )
                add_rule(sig_k, len(L)-1)
                if s != 0:
                    S += [len(L)-1]
                else:
                    if get_verbose() >= 0: print(f"S-Polynomial reduced to zero! {k} and {l}")
                    syzygies += [len(L)-1]
        S = sorted(S, key=lambda x: sig(x))
        return S

    def reduction(self, S, B, B_voo, Gprev, Gcurr):
        """
        adapted from Justin Gash (p.54ff):

        'Let's begin our discussion by focusing our attention to the
         outer layer of the reduction subroutine(s): reduction. The
         signed polynomial with smallest signature, denoted h, is
         grabbed and removed from the todo list of polynomial to be
         reduced.

         It's normal form, with respect to the previous Gröbner basis,
         and other information is sent to the sub-subroutine
         top_reduction. If top_reduction determines that the signed
         polynomial can be reduced, then nothing will be added to
         completed and the reduced (still signed) version of h will be
         placed back into todo. If no top reduction is possible, h is
         made monic by K-multiplication and the resulting signed
         polynomial is placed in completed.

         This description of reduction seems very similar to other
         reduction routines from other algorithms. The difference lies
         in the phrase, "If top_reduction determines that the signed
         polynomial can be reduced ..."'
        """
        phi = self.phi
        sig, voo, poly = self.sig, self.voo, self.poly
        get_sig = lambda x : self.get_sig_from_voo(voo(x))
        top_reduction = self.top_reduction
        voo_reduce = self.voo_reduce
        L = self.L

        assert all([b == phi(bv) for b, bv in zip(B, B_voo)]), f"Basis B and VoO's don't match."
        assert all([phi(s) == poly(s) for s in S]), f"reduction: Inconsistent voo and poly in S."
        assert all([sig(s) == get_sig(s) for s in S]), f"reduction: Inconsistent voo and sig in S."

        to_do = S
        completed = set()
        while to_do:
            k, to_do = to_do[0], to_do[1:]
            assert phi(k) == poly(k), f"reduction: Mismatching voo and poly at index {k}."
            assert sig(k) == get_sig(k), f"reduction: Mismatching voo and sig at index {k}."
            if get_verbose() >= 2: print(f"Processing {k} – {sig(k)}, {poly(k)}")
            h = poly(k).reduce(B)
            pol, voo_h = voo_reduce(k, B, B_voo)
            assert h == pol, f"\nBuggy behavior in 'voo_reduce'."
            if get_verbose() >= 2 and h != poly(k): print(f"Reduced {poly(k)} to {h}")
            L[k] = (sig(k), h, voo_h)
            newly_completed, redo = top_reduction(k, Gprev, Gcurr.union(completed))
            completed = completed.union( newly_completed )
            if get_verbose() >= 2 and k in newly_completed: print(f"completed {k} lm {poly(k).lt()}")
            to_do += redo
            to_do.sort(key=lambda x: sig(x))
        return completed

    def voo_reduce(self, i, basis, basis_voo):
        """
        Perform complete reduction of polynomial with index i by basis,
        and keep track of how that alters it's vector of origin.
        Returns the fully reduced polynomial and corresponding VoO.
        """
        sig, poly, voo = self.sig, self.poly, self.voo
        get_sig_from_voo = self.get_sig_from_voo
        phi = self.phi

        assert all([phi(b_voo) == b for b, b_voo in zip(basis, basis_voo)]), f"basis mismatches voos."
        assert phi(voo(i)) == poly(i), f"voo_reduce: Mismatching voo and poly at index {i}."
        assert sig(i) == get_sig_from_voo(voo(i)), f"voo_reduce: Mismatching voo and sig at index {i}."
        p = poly(i)
        v = voo(i)
        reduced = True
        while reduced:
            reduced = False
            for b, b_voo in zip(basis, basis_voo):
                p_before_b = p
                reduced_by_b = True # In order to mimic built-in reduction, first completely reduce using b
                while reduced_by_b:
                    reduced_by_b = False
                    quo, rem = p.quo_rem(b.lt())
                    if quo:
                        p = p - quo*b
                        v = v - quo*b_voo
                        reduced = True
                        reduced_by_b = True
                assert p_before_b.reduce([b]) == p, "voo_reduce: reduction of p by b seems wrong"
        assert phi(v) == p, f"voo_reduce: Kept track of voo incorrectly."
        return p, v

    def top_reduction(self, k, Gprev, Gcurr):
        """
        adapted from Justin Gash (p.55ff):

        'We will go through top_reduction step-by-step. If the signed
         polynomial being examined has polynomial part 0, then there
         is no data left in that particular signed polynomial - an
         empty ordered pair is returned. Otherwise top_reduction calls
         upon another sub-subroutine find_reductor. Essentially, if
         find_reductor comes back negative, the current signed
         polynomial is made monic and returned to reduction to be
         placed in completed. If a top-reduction is deemed possible,
         then there are two possible cases: either the reduction will
         increase the signature of polynomial or it won't. In the
         latter case, the signature of r_{k_0} is maintained, the
         polynomial portion is top-reduced and the signed polynomial
         is returned to reduction to be added back into todo; this
         case corresponds to top-reduction in previous algorithms.

         In the former case, however, the signature will change. This
         is marked by adding a new polynomial r_N (our notation here
         describes N after N was incremented) with appropriate
         signature based upon the reductor, not S(r_{k_0}). A new rule
         is added (as I mentioned previously, this will be explained
         later) and then both r_{k_0} and r_N are sent back to
         reduction to be added back into todo. This is done because
         r_N has a different signature than r_{k_0} and r_{k_0} might
         still be reducible by another signed polynomial.
        """
        sig, poly, voo = self.sig, self.poly, self.voo
        get_sig_from_voo = self.get_sig_from_voo
        find_reductor = self.find_reductor
        syzygies = self.syzygies
        add_rule = self.add_rule
        phi = self.phi
        L = self.L

        if poly(k) == 0:
            if get_verbose() >= 0: print(f"Reduction of {k} to zero.")
            self.zero_reductions += 1
            self.syzygies += [k]
            return set(), set()
        p = poly(k)
        p_voo = voo(k)
        J = find_reductor(k, Gprev, Gcurr)
        if not J:
            L[k] = ( sig(k), p/p.lc(), p_voo/p.lc() )
            assert phi(k) == poly(k), f"top_reduction: Mismatching voo and poly at index {k}."
            assert sig(k) == get_sig_from_voo(voo(k)), f"top_reduction: Mismatching voo and sig at index {k}."
            return {k}, set()
        j = J.pop()
        q = poly(j)
        q_voo = voo(j)
        u = p.lt()//q.lt()
        p = p - u*q
        p_voo = p_voo - u*q_voo
        self.reductions += 1
        if p != 0:
            p_voo /= p.lc()
            p /= p.lc()
        sig_j = u/u.lc() * sig(j)
        # no need to add k to syzygies below: calling function “reduction” will redo k
        if sig_j < sig(k):
            L[k] = (sig(k), p, p_voo)
            assert phi(k) == poly(k), f"top_reduction: Mismatching voo and poly at index {k}."
            assert sig(k) == get_sig_from_voo(voo(k)), f"top_reduction: Mismatching voo and sig at index {k}."
            return set(), {k}
        else:
            assert p == phi(p_voo), f"top_reduction: Mismatching voo and poly at index {len(L)-1}."
            assert sig_j == get_sig_from_voo(p_voo), f"top_reduction: Mismatching voo and sig at index {len(L)-1}."
            L.append( (sig_j, p, p_voo) )
            add_rule(sig_j, len(L)-1)
            return set(), {k, len(L)-1}

    def find_reductor(self, k, Gprev, Gcurr):
        """
        adapted from Justin Gash (p.56ff):

        'For a previously added signed polynomial in G_curr to become
         a reductor of r_{k_0}, it must meet four requirements:

         (1) u = HT(r_{k_0})/HT(r_{k_j}) ∈ T
         (2) NF(u_{t_j}, G_curr) = u_{t_j}
         (3) not is_rewriteable(u, r_{k_j})
         (4) u_{t_j} F_{k_j} = S(r_{k_0})

         We will go through each requirement one-by-one.

         Requirement (1) is simply the normal top-reduction
         requirement. The only thing of note here is that, in testing
         for the top-reducibility, u is assigned a particular value to
         be used in subsequent tests.

         Requirement (2) is making sure that the signature of the
         reductor is normalized.  Recall that we only want signatures
         of our polynomials to be normalized - we are discarding
         non-normalized S-polynomials. If we ignored this condition
         and our re- ductor wound up having larger signature than
         S(r_{k_0}), then top_reduction would create a new signed
         polynomial with our reductor's non-normalized signature. (We
         might add that, if the reductor had smaller signature than
         S(r_{k_0}), it would be fine to reduce by it; however, F5
         doesn't miss anything by forgoing this opportunity because,
         by Lemma 3.2.1 (The Normalization Lemma), there will be
         another normalized reductor with the same head term and
         smaller signature.)

         Requirement (3) will be discussed when we discuss
         is_rewriteable. That discussion is approaching rapidly.

         Requirement (4) is a check that makes sure we don't reduce by
         something that has the same signature as r_{k_0} . Recall
         that we want all signed polynomials used during the run of F5
         to be admissible. If we reduced by a polynomial that has the
         same signature, we would be left with a new polynomial for
         which we would have no idea what the signature is. The act of
         reduction would have certainly lowered the signature, thus
         causing admissibility to be lost. (We will comment on this
         requirement later in subsection 3.5. With a little care, we
         can loosen this requirement.)
        """
        is_rewritable = self.is_rewritable
        is_top_reducible = self.is_top_reducible

        poly = self.poly
        sig = self.sig
        t = poly(k).lt()
        for j in Gcurr:
            tprime = poly(j).lt()
            if divides(tprime,t):
                u = t // tprime
                mj, ej = sig(j)
                if u/u.lc() * sig(j) != sig(k) and not is_rewritable(u, j) and not is_top_reducible(u*mj, Gprev):
                    return set([j])
        return set()

    def is_top_reducible(self, t, Gprev):
        """
        Note, that this function test traditional top reduction and
        not top_reduction as implemented in the function with the same
        name of this class.
        """
        R = self.R
        poly = self.poly
        return any([R.monomial_divides(poly(g).lm(), t) for g in Gprev])

    def add_rule(self, s, k):
        self.Rules[s[1]].append( (s[0],k) )

    def is_rewritable(self, u, k):
        j = self.find_rewriting(u, k)
        return j != k

    def find_rewriting(self, u, k):
        """
        adapted from Justin Gash (p.57):

        'find_rewriting gives us information to be used as an
         additional criterion for eliminating critical pairs. Proof of
         this fact is given in section 3.4.3. In short, we could
         remove all discussion of rules and find_rewriting and F5
         would work fine. (But it would work much more slowly.) So we
         will treat these final four subroutines as a separate module
         that works in conjunction with F5, but is not an official
         part of the F5 criteria per se.'
        """
        Rules = self.Rules
        mk, v = self.sig(k)
        for ctr in reversed(range(len(Rules[v]))):
            mj, j = Rules[v][ctr]
            if divides(mj, u * mk):
                return j
        return k

class F5R(F5):
    def basis(self, F):
        """
        F5's main routine. Computes a Gröbner basis for F.

        INPUT:
            F -- a list of polynomials

        OUTPUT:
            G -- a list of polynomials; a Gröbner basis for <F>
        """
        poly = self.poly
        incremental_basis = self.incremental_basis

        self.__init__(F)

        Rules = self.Rules
        L = self.L

        m = len(F)
        F = sorted(F, key=cmp_to_key(compare_by_degree))

        f0 = F[0]
        L[0] = (Signature(1, 0), f0*f0.lc()**(-1))
        Rules.append(list())

        Gprev = set([0])
        B = [f0]

        for i in range(1, m):
            if get_verbose() >= 1: print(f"Starting incremental basis up to {i}")
            f = F[i]
            L.append( (Signature(1,i), f*f.lc()**(-1)) )
            Gcurr = incremental_basis(i, B, Gprev)
            if any(poly(lambd) == 1 for lambd in Gcurr):
                return set(1)
            Gprev = Gcurr
            B = Ideal([poly(l) for l in Gprev]).interreduced_basis()
            #B = self.interreduce([poly(l) for l in Gprev])

        return B

    def interreduce(self, RF):
        """
        interreduce RF and count the number of reductions performed.

        INPUT:
            RF -- a list of polynomial
        """
        F = list(RF)
        for each in range(len(F)):
           F[each] = F[each]*F[each].lc()**(-1)
        i = 0
        while i < len(F):
           reduceme = F.pop(0)
           reduced = False
           for j in range(len(F)):
              quo, rem = self.divide(reduceme,F[j])
              reduceme = rem
              if (quo != 0) and (rem != 0):
                 reduceme = rem*rem.lc()**(-1)
                 j = -1
                 reduced = True
           if (reduceme != 0):
              F.append(reduceme)
           if reduced:
              i = -1
           i = i + 1
        return F

    def normal_form(self, f, B):
        """
        Compute the normal form of f w.r.t. B and count the number of
        reductions.

        INPUT:
            f -- a polynomial
            B -- a set of polynomials
        """
        remainder = f
        quotient = [0 for each in B]
        i = 0
        while i < len(B):
           quo, rem = self.divide(remainder, B[i])
           remainder = rem
           if quo != 0:
              i = -1
           i = i + 1
        return remainder

    def divide(self,dividend, divisor):
        """
        Divide dividend by divisor and count number of reductions.

        INPUT:
            dividend -- a polynomial
            divisor -- a polynomial
        """
        remainder = dividend
        quotient = 0
        mons = remainder.monomials()
        coeffs = remainder.coefficients()
        t = divisor.lm()
        c = divisor.lc()
        i = 0
        while (remainder != 0) and (i < len(mons)):
           if t.divides(mons[i]):
              self.reductions += 1
              quotient = quotient + (mons[i]/t*coeffs[i]/c).numerator()
              remainder = remainder - (mons[i]/t*coeffs[i]/c*divisor).numerator()
              mons = remainder.monomials()
              coeffs = remainder.coefficients()
           else:
              i = i + 1
        return quotient, remainder

## See Eder & Faugère (2014), Section 8.2 “F5C – Improved S-pair generation”
## To enhance with VoOs, keep track of all re-definitions of mapping φ.
class F5C(F5):
    def basis(self, F):
        """
        F5's main routine. Computes a Gröbner basis for F.

        INPUT:
            F -- a list of polynomials

        OUTPUT:
            G -- a list of polynomials; a Gröbner basis for <F>
        """
        incremental_basis = self.incremental_basis
        voo = self.voo
        poly = self.poly

        self.__init__(F)

        Rules = self.Rules
        R = self.R
        L = self.L

        m = len(F)

        f = F[0]
        L[0] = (Signature(1, 0), f/f.lc(), unit_vec(R, 0, m)/f.lc())
        Rules.append(list())

        Gprev = {0}
        B = [poly(0)]
        B_voo = [voo(0)]

        for i in range(1, m):
            if get_verbose() >= 1: print(f"Starting incremental basis up to {i}")
            f = F[i]
            L.append( (Signature(1, len(L)), f/f.lc(), unit_vec(R, len(L), m)/f.lc())  )
            Gcurr = incremental_basis(len(L)-1, B, B_voo, Gprev)
            for j in range(len(Gcurr)):
                if poly(j) == 1: return [poly(j)], [voo(j)]
            B = Ideal([poly(l) for l in Gcurr]).interreduced_basis()
            B_2, B_voo = self.interreduce_basis_voo([poly(l) for l in Gcurr], B_voo)
            assert all([b in B_2 for b in B]) and all([b in B for b in B_2]), "Buggy behavior in interreduce_basis_voo."
            if i != m-1:
                Gprev = self.setup_reduced_basis(B)
        return B, B_voo

    def setup_reduced_basis(self, B):
        """
        Update the global L and Rules to match the reduced basis B.

        OUTPUT:
            Gcurr -- index set for B
        """
        raise NotImplementedError("Need to keep track of the basis changes – See Eder & Faugère (2014), Sec. 8.2")
        add_rule = self.add_rule
        Rules = self.Rules
        m = len(self.F)
        L = self.L
        R = self.R

        # we don't want to replace L but modify it
        L[:] = [(Signature(1, i), f, "some vector") for i, f in enumerate(B)]
        Rules[:] = [[]] * len(B)
        Gcurr = set()

        for i,f in enumerate(B):
            Gcurr += [i]
            t = B[i].lt()
            for j in range(i+1, len(B)):
                fjlt = B[j].lt()
                u = LCM(t, fjlt)//fjlt
                add_rule( Signature(u, j), -1 )
        return Gcurr

class F4F5(F5C):
    """
    F4-Style F5

    Till Steger's calls this F4.5. We don't know how Jean-Charles
    Faugère calls it.
    """
    def reduction(self, S, B, Gprev, Gcurr):
        """
        INPUT:
            S -- a list of components of S-polynomials
            B -- ignored
            Gprev -- the previous Gröbner basis indexed in L
            Ccurr -- the Gröbner basis computed so far indexed in L
        """
        L = self.L
        add_rule = self.add_rule
        poly = self.poly

        S = self.symbolic_preprocessing(S, Gprev, Gcurr)
        St = self.gauss_elimination(S)

        Ret = []

        for k, (s,p,idx) in enumerate(St):
            if (s,p,idx) == S[k] and idx == -1:
                continue # ignore unchanged new polynomials
            if idx >= 0:
                L[idx] = L[idx][0], p # update p
                if p != 0:
                    Ret.append(idx)
            else:
                L.append( (s,p) ) # we have a new polynomial
                add_rule( s, len(L)-1 )
                if p != 0:
                    Ret.append(len(L)-1)
        return Ret

    def symbolic_preprocessing(self,S, Gprev, Gcurr):
        """
        Add polynomials to the set S such that all possible reductors
        for all elements in S are available.

        INPUT:
            S -- a list of components of S-polynomials
            Gprev -- the previous Gröbner basis indexed in L
            Ccurr -- the Gröbner basis computed so far indexed in L
        """
        poly = self.poly
        L = self.L
        find_reductor = self.find_reductor

        # We add a new marker for each polynomial which encodes
        # whether the polynomial was added by this routine or is an
        # original input polynomial.

        F = [L[k]+(k,) for k in S]
        Done = set()

        # the set of all monomials
        M = set([m for f in F for m in f[1].monomials()])
        while M != Done:
            m = M.difference(Done).pop()
            Done.add(m)
            t, g = find_reductor(m, Gprev, Gcurr)
            if t!=0:
                F.append( (t*g[0], t*g[1], -1) )
                M = M.union((t*g[1]).monomials())
        return sorted(F, key=lambda f: f[0]) # sort by signature

    def find_reductor(self, m, Gprev, Gcurr):
        r"""
        Find a reductor $g_i$ for $m$ in $G_{prev}$ and $G_{curr}$ subject
        to the following contraint.  is a
         * the leading monomial of $g_i$ divides $m$
         * $g_i ∈ G_{prev}$ is preferred over $g_i ∈ G_{curr}$
         * if $g_i ∈ G_{curr}$ then
           * $g_i$ is not rewritable
           * $g_i$ is not top reducible by $G_prev$

        INPUT:
            m -- a monomial
            Gprev -- the previous Gröbner basis indexed in L
            Ccurr -- the Gröbner basis computed so far indexed in L

        """
        is_rewritable = self.is_rewritable
        is_top_reducible = self.is_top_reducible
        sig = self.sig
        poly = self.poly

        L = self.L
        R = m.parent()

        for k in Gprev:
            if R.monomial_divides(poly(k).lm(),m):
                return  R.monomial_quotient(m,poly(k).lm()), L[k]
        for k in Gcurr:
            if R.monomial_divides(poly(k).lm(),m):
                t =  R.monomial_quotient(m,poly(k).lm())
                if is_rewritable(t, k):
                    continue
                if is_top_reducible(t/t.lc() * sig(k)[0], Gprev):
                    continue
                return t, L[k]
        return 0, -1


    def gauss_elimination(self, F1):
        """
        Perform permuted Gaussian elimination on F1.

        INPUT:
            F1 -- a list of tuples (sig, poly, idx)
        """
        F = [f[1] for f in F1]
        if len(F) == 0:
            return F
        A,v = Sequence(F).coefficient_matrix()
        self.zero_reductions += A.nrows()-A.rank()
        if get_verbose() >= 1: print(f"{A.nrows():>4} x {A.ncols():>4}, {A.rank():>4}, {A.nrows()-A.rank():>4}")
        nrows, ncols = A.nrows(), A.ncols()
        for c in range(ncols):
            for r in range(0,nrows):
                if A[r,c] != 0:
                    if any(A[r,i] for i in range(c)):
                        continue
                    a_inverse = ~A[r,c]
                    A.rescale_row(r, a_inverse, c)
                    for i in range(r+1,nrows):
                        if A[i,c] != 0:
                            if any(A[i,_] for _ in range(c)):
                                continue
                            minus_b = -A[i,c]
                            A.add_multiple_of_row(i, r, minus_b, c)
                    break
        F = (A*v).list()
        return [(F1[i][0],F[i],F1[i][2]) for i in range(len(F))]

class F5SansRewriting(F5):
    """
    A variant of F5 which does not use the rewriting rule. This is
    motivated by the following observation by Justin Gash (p.57):

    'Rewritten gives us information to be used as an additional
     criterion for eliminating critical pairs. Proof of this fact is
     given in section 3.4.3. In short, we could remove all discussion
     of rules and Rewritten and F5 would work fine. (But it would work
     much more slowly.) So we will treat these final four subroutines
     as a separate module that works in conjunction with F5, but is
     not an official part of the F5 criteria per se.'
    """
    def critical_pair(self, k, l, i, Gprev):
        """
        Like the function of base class F5, without the
        rewriteability checks.
        """
        poly = self.poly
        sig = self.sig
        is_top_reducible = self.is_top_reducible

        tk = poly(k).lt()
        tl = poly(l).lt()
        t = LCM(tk, tl)
        u0 = t//tk
        u1 = t//tl
        m0, e0 = sig(k)
        m1, e1 = sig(l)
        if e0 == e1 and u0*m0 == u1*m1:
            return set()

        # Stegers and Gash leave out the == i check, Faugere and Perry
        # have it. It is unclear for now, whether the check is
        # necessary.
        if e0 == i and is_top_reducible(u0*m0, Gprev):
            return set()
        if e1 == i and is_top_reducible(u1*m1, Gprev):
            return set()
        if u0/u0.lc() * sig(k) < u1/u1.lc() * sig(l):
            u0, u1 = u1, u0
            k, l = l, k
        return set([(t,k,u0,l,u1)])

    def compute_spols(self, P):
        """
        Like the function of base class F5, without the
        rewriteability checks.
        """
        syzygies = self.syzygies
        sig = self.sig
        poly = self.poly
        voo = self.voo
        L = self.L

        S = list()
        P = sorted(P, key=lambda x: x[0])
        for (t,k,u,l,v) in P:
            s = u*poly(k)-v*poly(l) # S-Polynomial
            s_voo = u*voo(k)-v*voo(l) # S-Polynomial
            if s != 0:
                s_voo /= s.lc()
                s /= s.lc() # normalize
            L.append( (u/u.lc() * sig(k), s, s_voo) )
            if s != 0:
                S += [len(L)-1]
            else:
                syzygies += [len(L)-1]
        S = sorted(S, key=lambda x: sig(x))
        return S

    def top_reduction(self, k, Gprev, Gcurr):
        """
        Like the function of base class F5, without the
        adding of rules.
        """
        find_reductor = self.find_reductor
        syzygies = self.syzygies
        sig = self.sig
        poly = self.poly
        voo = self.voo
        phi = self.phi
        L = self.L

        if poly(k) == 0:
            if get_verbose() >= 0: print(f"Reduction of {k} to zero.")
            self.zero_reductions += 1
            self.syzygies += [k]
            return set(), set()
        p = poly(k)
        p_voo = voo(k)
        J = find_reductor(k, Gprev, Gcurr)
        if not J:
            L[k] = ( sig(k), p/p.lc() , p_voo/p.lc() )
            assert phi(k) == poly(k), f" [!] Something wrong:\n{voo(k)}\n{poly(k)}."
            return set([k]), set()
        j = J.pop()
        q = poly(j)
        q_voo = voo(j)
        u = p.lt()//q.lt()
        p = p - u*q
        p_voo = p_voo - u*q_voo
        self.reductions += 1
        if p != 0:
            p_voo /= p.lc()
            p /= p.lc()
        # no need to add k to syzygies below: calling function “reduction” will redo k
        if u/u.lc() * sig(j) < sig(k):
            L[k] = (sig(k), p, p_voo)
            assert phi(k) == poly(k), f" [!] Something wrong:\n{voo(k)}\n{poly(k)}."
            return set(), set([k])
        else:
            L.append((u/u.lc() * sig(j), p, p_voo))
            assert phi(-1) == poly(-1), f" [!] Something wrong:\n{voo(-1)}\n{poly(-1)}."
            return set(), set([k, len(L)-1])

    def find_reductor(self, k, Gprev, Gcurr):
        """
        Like the function of base class F5, without the
        rewriteability check.
        """
        is_top_reducible = self.is_top_reducible

        poly = self.poly
        sig = self.sig
        t = poly(k).lt()
        for j in Gcurr:
            tprime = poly(j).lt()
            if divides(tprime,t):
                u = t // tprime
                mj, ej = sig(j)
                if u/u.lc() * sig(j) != sig(k) and not is_top_reducible(u*mj, Gprev):
                    return set([j])
        return set()

from collections import UserList

class Signature(UserList):
    def __init__(self, multiplier, index):
        """
        Create a new signature from the mulitplier and the index.
        """
        UserList.__init__(self, (multiplier, index))

    def __lt__(self, other):
        if self[1] != other[1]:
            return self[1] < other[1]
        return self[0] < other[0]

    def __gt__(self, other):
        if self[1] != other[1]:
            return self[1] > other[1]
        return self[0] > other[0]

    def __eq__(self, other):
        return self[0] == other[0] and self[1] == other[1]

    def __neq__(self, other):
        return self[0] != other[0] or self[1] != other[1]

    def __rmul__(self, other):
        if isinstance(self, Signature):
            return Signature(other * self[0], self[1])
        raise TypeError

    def __hash__(self):
        return hash(self[0]) + hash(self[1])
