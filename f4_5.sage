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
LM = lambda f: f.lm()
LT = lambda f: f.lt()

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

    def phi(self, i):
        """
        Maps vector of origin at index 'i' to its polynomial under
        the input system. Always results in the polynomial at 'i'.
        Retuns the polynomial, allowing to check consistency of voos.
        """
        return self.voo(i) * vector(self.R, self.F)

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

        self.__init__(F)

        Rules = self.Rules
        L = self.L

        m = len(F)

        f = F[0]
        L[0] = (Signature(1, 0), f/f.lc(), unit_vec(R, 0, m)/f.lc())
        if phi(0) != poly(0): print(f" [!] Something wrong:\n{voo(0)}\n{poly(0)}.")
        Rules.append([])

        Gprev = set([0])
        B = [poly(0)]
        B_voo = [voo(0)]

        for i in range(1,m):
            if get_verbose() >= 0: print(f"Increment {i}")
            f = F[i]
            L.append( (Signature(1,i), f/f.lc(), unit_vec(R, i, m)/f.lc()) )
            Gcurr = incremental_basis(i, B, B_voo, Gprev)
            for j in range(len(Gcurr)):
                if poly(j) == 1:
                    return [poly(j)], [voo(j)]
            Gprev = Gcurr
            B = [poly(l) for l in Gprev]
            B_voo = [voo(l) for l in Gprev]
        return [poly(l) for l in Gprev], [voo(l) for l in Gprev]

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
        Rules.append( list() )

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
        # necessary
        if is_rewritable(u0, k) or is_rewritable(u1, l):
            return set()
        if u0 * sig(k) < u1 * sig(l):
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
        voo = self.voo
        poly = self.poly
        sig = self.sig
        is_rewritable = self.is_rewritable
        syzygies = self.syzygies
        add_rule = self.add_rule

        L = self.L

        S = list()
        P = sorted(P, key=lambda x: x[0])
        for (t,k,u,l,v) in P:
            if not is_rewritable(u,k) and not is_rewritable(v,l):
                s = u*poly(k)-v*poly(l) # S-Polynomial
                s_voo = u*voo(k)-v*voo(l)
                if s != 0:
                    s_voo /= s.lc() # normalize
                    s /= s.lc()
                L.append( (u * sig(k), s, s_voo) )
                add_rule(u * sig(k), len(L)-1)
                if s != 0:
                    S += [len(L)-1]
                else:
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
        F = self.F
        L = self.L
        phi = self.phi
        sig = self.sig
        voo = self.voo
        poly = self.poly
        top_reduction = self.top_reduction
        voo_reduce = self.voo_reduce

        to_do = S
        completed = set()
        while to_do:
            k, to_do = to_do[0], to_do[1:]
            if get_verbose() >= 2: print(f"Processing {k} – {sig(k)}, {poly(k)}")
            h = poly(k).reduce(B)
            pol, voo_h = voo_reduce(k, B, B_voo)
            assert h == pol, f"\nBuggy behavior in 'voo_reduce':\n    {h}\n        !=\n    {pol}"
            if get_verbose() >= 1 and pol != poly(k): print(f"Reduced {poly(k)} to {pol}")
            L[k] = (sig(k), pol, voo_h)
            if phi(k) != poly(k): print(f" [!] Something wrong:\n{voo(k)}\n{poly(k)}.")
            newly_completed, redo = top_reduction(k, Gprev, Gcurr.union(completed))
            completed = completed.union( newly_completed )
            if get_verbose() >= 2 and k in newly_completed: print(f"completed {k} lm {poly(k).lt()}")
            to_do += redo
            to_do.sort(key=lambda x: sig(x))
        return completed

    def voo_reduce(self, i, base, base_voo):
        """
        Perform complete reduction of polynomial with index i by base,
        and keep track of how that alters it's vector of origin.
        Returns the fully reduced polynomial and corresponding VoO.
        """
        pol = self.poly(i)
        voo = self.voo(i)
        reduced = True
        while reduced:
            reduced = False
            for b, b_voo in zip(base, base_voo):
                quo, rem = pol.quo_rem(b.lt())
                if quo:
                    pol = pol - quo*b # pol == rem
                    voo = voo - quo*b_voo
                    reduced = True
                # assert rem == pol - quo*b
                # if b.lt().divides(pol.lt()):
                #     div = pol.lt()//b.lt()
                #     pol -= div*b
                #     voo -= div*b_voo
                #     reduced = True
        return pol, voo

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
        find_reductor = self.find_reductor
        syzygies = self.syzygies
        add_rule = self.add_rule
        phi = self.phi
        voo = self.voo
        poly = self.poly
        sig = self.sig
        L = self.L
        F = self.F

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
            if phi(k) != poly(k): print(f" [!] Something wrong:\n{voo(k)}\n{poly(k)}.")
            return set([k]),set()
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
        if u * sig(j) < sig(k):
            L[k] = (sig(k), p, p_voo)
            if phi(k) != poly(k): print(f" [!] Something wrong:\n{voo(k)}\n{poly(k)}.")
            return set(), set([k])
        else:
            L.append( (u * sig(j), p, p_voo) )
            if phi(-1) != poly(-1): print(f" [!] Something wrong:\n{voo(-1)}\n{poly(-1)}.")
            add_rule(u * sig(j), len(L)-1)
            return set(), set([k, len(L)-1])

    def find_reductor(self, k, Gprev, Gcurr):
        """
        adapted from Justin Gash (p.56ff):

        'For a previously added signed polynomial in G_curr to become
         a reductor of r_{k_0}, it must meet four requirements:

         (1) u = HT(r_{k_0})/HT(r_{k_j}) \in T
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
                if u * sig(j) != sig(k) and not is_rewritable(u, j) \
                        and not is_top_reducible(u*mj, Gprev):
                    return set([j])
        return set()

    def is_top_reducible(self, t, l):
        """
        Note, that this function test traditional top reduction and
        not top_reduction as implemented in the function with the same
        name of this class.
        """
        R = self.R
        poly = self.poly
        for g in l:
            if R.monomial_divides(poly(g).lm(),t):
                return True
        return False

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
            print(f"Increment {i}")
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
        poly = self.poly

        self.__init__(F)

        Rules = self.Rules
        L = self.L

        m = len(F)
        F = sorted(F, key=cmp_to_key(compare_by_degree))

        f0 = F[0]
        L[0] = (Signature(1, 0), f0*f0.lc()**(-1))
        Rules.append(list())

        Gprev = set([0])
        B = set([f0])

        for i in range(1, m):
            print(f"Increment {i}")
            f = F[i]
            L.append( (Signature(1,len(L)), f*f.lc()**(-1)) )
            Gcurr = incremental_basis(len(L)-1, B, Gprev)
            if any(poly(lambd) == 1 for lambd in Gcurr):
                return set(1)
            B = Ideal([poly(l) for l in Gcurr]).interreduced_basis()
            #B = self.interreduce([poly(l) for l in Gcurr])
            if i != m-1:
                Gprev = self.setup_reduced_basis(B)
        return B

    def setup_reduced_basis(self, B):
        """
        Update the global L and Rules to match the reduced basis B.

        OUTPUT:
            Gcurr -- index set for B
        """
        add_rule = self.add_rule
        L = self.L
        Rules = self.Rules

        # we don't want to replace L but modify it
        L[:] = [(Signature(1,i), f) for i,f in enumerate(B)]
        Rules[:] = [[] for _ in range(len(B))]
        Gcurr = set()

        for i,f in enumerate(B):
            Gcurr.add( i )
            t = B[i].lt()
            for j in range(i+1, len(B)):
                fjlt = B[j].lt()
                u = LCM(t, fjlt)//fjlt
                add_rule( Signature(u, j), -1 )
        return Gcurr


class F4F5(F5C):
#class F4F5(F5):
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
         * $g_i \in G_{prev}$ is preferred over $g_i \in G_{curr}$
         * if $g_i \in G_{curr}$ then
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
                if is_top_reducible(t * sig(k)[0], Gprev):
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
        print(f"{A.nrows():>4} x {A.ncols():>4}, {A.rank():>4}, {A.nrows()-A.rank():>4}")
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
        if u0 * sig(k) < u1 * sig(l):
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
            L.append( (u * sig(k), s, s_voo) )
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
            if phi(k) != poly(k): print(f" [!] Something wrong:\n{voo(k)}\n{poly(k)}.")
            return set([k]), set()
        j = J.pop()
        q = poly(j)
        q_voo = voo(j)
        u = p.lt()//q.lt()
        print(f"\n{k}, {j}")
        p = p - u*q
        p_voo = p_voo - u*q_voo
        self.reductions += 1
        if p != 0:
            p_voo /= p.lc()
            p /= p.lc()
        # no need to add k to syzygies below: calling function “reduction” will redo k
        if u * sig(j) < sig(k):
            L[k] = (sig(k), p, p_voo)
            if phi(k) != poly(k): print(f" [!] Something wrong:\n{voo(k)}\n{poly(k)}.")
            return set(), set([k])
        else:
            L.append((u * sig(j), p, p_voo))
            if phi(-1) != poly(-1): print(f" [!] Something wrong:\n{voo(-1)}\n{poly(-1)}.")
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
                if u * sig(j) != sig(k) and not is_top_reducible(u*mj, Gprev):
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

f5  = F5()
f5r = F5R()
f5c = F5C()
ff  = F4F5()


if __name__ == "__main__":

    set_verbose(1)
    R = PolynomialRing(GF(283), 'x', 6, order="lex")
    R.inject_variables()
    polys = [
        # x0 + x1,
        # x0*x1,

        # x0 - x2,
        # x0^2 + x2^2 - x1,
        # x0^2 + x1^2 + x2^2 - 1,

        # x0^5 + x1^4 + x2^3 - 1,
        # x0^3 + x1^2 + x2^2 - 1,

        # Has one syzygy vector with signature (x1*x3, 2)
        # x0*x1 - x1*x3,
        # x0^2  - x2*x3,
        # x2^3  - x3^3,

        # p = 283
        -x0 + 15*x1 + x2,
        -x3 + 15*x4 + x5,
        x1^19 + 112*x1^18 - 18*x1^17 - 10*x1^16 - 67*x1^15 - 88*x1^14 - 25*x1^13 - 70*x1^12 + 40*x1^11 + 24*x1^10 - 98*x1^9 - 69*x1^8 + 20*x1^7 + 61*x1^6 + 111*x1^5 + 105*x1^4 - 110*x1^3 - 35*x1^2 + 116*x1,
        95*x1^18 + 30*x1^17 + 67*x1^16 - 109*x1^15 - 72*x1^14 + 5*x1^13 + 13*x1^12 + 20*x1^11 + 127*x1^10 + 139*x1^9 - 131*x1^8 - 5*x1^7 - 129*x1^6 + 10*x1^5 - 78*x1^4 - 14*x1^3 + x1^2 + 28*x1 - x4 + 9,
        x2^15 - 105*x2^14 - 89*x2^13 - 127*x2^12 + 119*x2^11 - 140*x2^10 - 49*x2^9 + 91*x2^8 + 123*x2^7 + 53*x2^6 + 63*x2^5 - 60*x2^4 - 137*x2^3 - 9*x2^2 - 17*x2,
        103*x2^14 - 117*x2^13 + 14*x2^12 - 26*x2^11 - 121*x2^10 - 78*x2^9 + 69*x2^8 + 5*x2^7 - 74*x2^6 - 86*x2^5 - 40*x2^4 - 23*x2^3 - 106*x2^2 - 89*x2 - x5 + 9,

        # p = 347
        # -x0 + 12*x1 + x2,
        # -x3 + 12*x4 + x5,
        # x1^29 - 59*x1^28 + 139*x1^27 - 123*x1^26 + 40*x1^25 - 47*x1^24 + 168*x1^23 - 25*x1^22 + 29*x1^21 - 37*x1^20 + 77*x1^19 - 145*x1^18 - 173*x1^17 - 134*x1^16 + 114*x1^15 - 115*x1^14 + 57*x1^13 - 82*x1^12 + 164*x1^11 + 142*x1^10 - 166*x1^9 + 139*x1^8 + 10*x1^7 - 42*x1^6 + 127*x1^5 + 85*x1^4 - 115*x1^3 + 26*x1^2 - 55*x1,
        # 130*x1^28 + 67*x1^27 + 123*x1^26 + 105*x1^25 + 93*x1^24 + 61*x1^23 + 148*x1^22 - 51*x1^21 + 64*x1^20 - 43*x1^19 + 163*x1^18 + 100*x1^17 + 165*x1^16 - 107*x1^15 - 94*x1^14 + 115*x1^13 + 89*x1^12 + 48*x1^11 - 123*x1^10 + 121*x1^9 - 65*x1^8 + 104*x1^7 - 59*x1^6 - 26*x1^5 - 11*x1^4 - 33*x1^3 + 25*x1^2 - 70*x1 - x4 + 6,
        # x2^12 - 66*x2^11 - 157*x2^10 - 52*x2^9 + 13*x2^8 - 11*x2^7 + 161*x2^6 + 161*x2^5 + 137*x2^4 - 42*x2^3 - 143*x2^2 - 2*x2,
        # -75*x2^11 - 83*x2^10 + 162*x2^9 - 67*x2^8 + 141*x2^7 + 17*x2^6 - 89*x2^5 + 110*x2^4 - 151*x2^3 + 43*x2^2 - 10*x2 - x5 + 6,

        # p = 911
        # -x0 + 24*x1 + x2,
        # -x3 + 24*x4 + x5,
        # x1^38 + 208*x1^37 - 365*x1^36 + 437*x1^35 + 98*x1^34 + 343*x1^33 + 62*x1^32 + 434*x1^31 + 244*x1^30 + 140*x1^29 - 348*x1^28 + 239*x1^27 + 52*x1^26 - 368*x1^25 + 391*x1^24 - 433*x1^23 - 435*x1^22 + 387*x1^21 + 64*x1^20 - 275*x1^19 + 39*x1^18 + 289*x1^17 - 281*x1^16 - 419*x1^15 - 178*x1^14 - 18*x1^13 - 408*x1^12 + 183*x1^11 + 441*x1^10 + 16*x1^9 - 126*x1^8 - 58*x1^7 + 109*x1^6 - 129*x1^5 - 17*x1^4 - 384*x1^3 + 437*x1^2 - 372*x1,
        # 441*x1^37 - 131*x1^36 + 43*x1^35 + 264*x1^34 - 61*x1^33 - 402*x1^32 + 273*x1^31 - x1^30 + 316*x1^29 - 441*x1^28 + 48*x1^27 + 452*x1^26 - 167*x1^25 - 230*x1^24 - 149*x1^23 + 114*x1^22 + 403*x1^21 + 381*x1^20 - 114*x1^19 - 264*x1^18 - 397*x1^17 - 278*x1^16 + 134*x1^15 + 93*x1^14 - 292*x1^13 - 204*x1^12 + 309*x1^11 + 182*x1^10 + 427*x1^9 + 330*x1^8 + 448*x1^7 - 270*x1^6 - 192*x1^5 + 214*x1^4 + 285*x1^3 - 377*x1^2 - 275*x1 - x4 + 5,
        # x2^24 - 276*x2^23 + 397*x2^22 - 267*x2^21 + 13*x2^20 - 447*x2^19 + 219*x2^18 - 206*x2^17 + 119*x2^16 + 174*x2^15 + 225*x2^14 + 283*x2^13 - 14*x2^12 - 129*x2^11 - 209*x2^10 + 335*x2^9 - 359*x2^8 - 439*x2^7 - 276*x2^6 + 358*x2^5 - 455*x2^4 - 304*x2^3 + 181*x2^2 + 165*x2,
        # 247*x2^23 + 301*x2^22 + 7*x2^21 - 218*x2^20 - 289*x2^19 - 257*x2^18 - 439*x2^17 + 424*x2^16 + 345*x2^15 + 329*x2^14 - 24*x2^13 + 38*x2^12 + 318*x2^11 - 230*x2^10 - 170*x2^9 + 15*x2^8 - 304*x2^7 - 389*x2^6 - 141*x2^5 + 59*x2^4 - 14*x2^3 + 274*x2^2 + 119*x2 - x5 + 5,
    ]
    if get_verbose() >= 2: print(f"––––––––––––\n Input:\n{polys}\n––––––––––––")
    gb, voos = f5(polys, homogenize=False)

    products = [voo * vector(R, polys) for voo in voos]
    same = [products[i] == gb[i] for i in range(len(gb))]

    if get_verbose() >= 0: print(f"––––––––––––\n is Gröbner basis: {Ideal(gb).basis_is_groebner()}")
    if get_verbose() >= 2: print(f"––––––––––––\n GB:\n{gb}")
    if get_verbose() >= 2: print(f"––––––––––––\n VoOs:\n{voos}")
    if get_verbose() >= 2:
        print(f"––––––––––––\n Products of Vectors of Origin and input F:")
        col_1_len = max([len(str(voo)) for voo in voos])
        col_2_len = max([len(str(pr)) for pr in products])
        [print(f"{'   ' if same[i] else '[!]'} {str(voos[i]):>{col_1_len}} * F = {str(products[i]):{col_2_len}} {'==' if same[i] else '!='} {gb[i]}") for i in range(len(gb))]
    if get_verbose() >= 1 and all(same): print(f"––––––––––––\n \\o/ VoOs are all good! \\o/")
    assert all(same), "Some Vector of Origin is misbehaving…"
    if get_verbose() >= 1:
        involvement_buckets = [0] * 2**(len(voos[0]))
        for voo in voos:
            idx = sum([2**i for i, x in enumerate(voo) if x])
            involvement_buckets[idx] += 1
        print(f"––––––––––––\n Involvement:         {involvement_buckets}")
    if get_verbose() >= 1: print(f"––––––––––––\n Zero reductions:     {f5.zero_reductions}")
    if get_verbose() >= 1:
        print(f"––––––––––––\n Non-Koszul Syzygies: {'None' if not f5.syzygies else ''}")
        for i in f5.syzygies:
            print(f"                      {f5.voo(i)}")
    with open("./voo_file.txt", 'w') as voo_file:
        for voo in voos:
            voo_file.write(f"{voo}\n")
