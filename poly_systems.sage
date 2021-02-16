load('f4_5.sage')

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

f5  = F5() # or: FR5, F5C, F4F5

R = PolynomialRing(GF(283), 'x', 6)
R.inject_variables()

polys = [ # contains a bunch of test-systems, each seperated by a new line
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

    # Definitely has more than one non-Koszul syzygy vector
    # x0*x1 - x1*x3,
    # x0^2  - x2*x3,
    # x1*x2^3  - x1*x3^3,
    # x2*x2^3  - x2*x3^3,
    # x2*x2^3  - x2*x3^3,

    # Eder & Faugère, example 4.3
    x0*x1*x2 - x2^2*x3,
    x0^2*x1 - x1^3,
    x1^3 - x2*x3^2 - x3^3,


    ## Bar
    # p = 283
    # -x0 + 15*x1 + x2,
    # -x3 + 15*x4 + x5,
    # x1^19 + 112*x1^18 - 18*x1^17 - 10*x1^16 - 67*x1^15 - 88*x1^14 - 25*x1^13 - 70*x1^12 + 40*x1^11 + 24*x1^10 - 98*x1^9 - 69*x1^8 + 20*x1^7 + 61*x1^6 + 111*x1^5 + 105*x1^4 - 110*x1^3 - 35*x1^2 + 116*x1,
    # 95*x1^18 + 30*x1^17 + 67*x1^16 - 109*x1^15 - 72*x1^14 + 5*x1^13 + 13*x1^12 + 20*x1^11 + 127*x1^10 + 139*x1^9 - 131*x1^8 - 5*x1^7 - 129*x1^6 + 10*x1^5 - 78*x1^4 - 14*x1^3 + x1^2 + 28*x1 - x4 + 9,
    # x2^15 - 105*x2^14 - 89*x2^13 - 127*x2^12 + 119*x2^11 - 140*x2^10 - 49*x2^9 + 91*x2^8 + 123*x2^7 + 53*x2^6 + 63*x2^5 - 60*x2^4 - 137*x2^3 - 9*x2^2 - 17*x2,
    # 103*x2^14 - 117*x2^13 + 14*x2^12 - 26*x2^11 - 121*x2^10 - 78*x2^9 + 69*x2^8 + 5*x2^7 - 74*x2^6 - 86*x2^5 - 40*x2^4 - 23*x2^3 - 106*x2^2 - 89*x2 - x5 + 9,

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


    ## Rescue XLIX, p = 101
    # 1 round
    # x2,
    # x0^3 + 14*x1^3 + 39*x2^3 + 45*x3^3 + 32*x3^2 + 42*x3 - 35,
    # x0^3 - 16*x1^3 + 23*x2^3 - x3^3 + 15*x3^2 + 26*x3 + 31,

    # 2 rounds
    # x2,
    # x0^3 + 14*x1^3 + 39*x2^3 - 29*x3^3 + 9*x3^2*x4 + 13*x3*x4^2 + 10*x4^3 - 26*x3^2*x5 - 19*x3*x4*x5 + 48*x4^2*x5 + 5*x3*x5^2 - 4*x4*x5^2 + 45*x5^3 - 49*x3^2 + 38*x3*x4 + 5*x4^2 - 20*x3*x5 + 16*x4*x5 + 33*x5^2 + 20*x3 - 16*x4 + 35*x5 - 22,
    # x0^3 - 16*x1^3 + 23*x2^3 + 31*x3^3 + 10*x3^2*x4 + 38*x3*x4^2 + x4^3 - 10*x3^2*x5 + 25*x3*x4*x5 - 3*x4^2*x5 + 38*x3*x5^2 + 3*x4*x5^2 - x5^3 - 43*x3^2 - 44*x3*x4 - 23*x4^2 + 44*x3*x5 + 46*x4*x5 - 23*x5^2 + 34*x3 + 8*x4 - 8*x5 + 9,
    # x0^3 + 10*x1^3 + 2*x2^3 + 8*x3^3 - 2*x3^2*x4 + 17*x3*x4^2 + 36*x4^3 + 41*x3^2*x5 + 10*x3*x4*x5 + 8*x4^2*x5 + 49*x3*x5^2 + 38*x4*x5^2 - 24*x5^3 - 19*x3^2 + 20*x3*x4 + 16*x4^2 - 6*x3*x5 - 50*x4*x5 - 43*x5^2 - 6*x3 - 50*x4 + 15*x5 + 8,
    # x3^3 + 14*x4^3 + 39*x5^3 + 45*x6^3 + 47*x6^2 - 39*x6 - 43,
    # x3^3 - 16*x4^3 + 23*x5^3 - x6^3 + 43*x6^2 - 44*x6 - 45,

    # 3 rounds
    # x2,
    # x0^3 + 14*x1^3 + 39*x2^3 - 29*x3^3 + 9*x3^2*x4 + 13*x3*x4^2 + 10*x4^3 - 26*x3^2*x5 - 19*x3*x4*x5 + 48*x4^2*x5 + 5*x3*x5^2 - 4*x4*x5^2 + 45*x5^3 - 49*x3^2 + 38*x3*x4 + 5*x4^2 - 20*x3*x5 + 16*x4*x5 + 33*x5^2 + 20*x3 - 16*x4 + 35*x5 - 22,
    # x0^3 - 16*x1^3 + 23*x2^3 + 31*x3^3 + 10*x3^2*x4 + 38*x3*x4^2 + x4^3 - 10*x3^2*x5 + 25*x3*x4*x5 - 3*x4^2*x5 + 38*x3*x5^2 + 3*x4*x5^2 - x5^3 - 43*x3^2 - 44*x3*x4 - 23*x4^2 + 44*x3*x5 + 46*x4*x5 - 23*x5^2 + 34*x3 + 8*x4 - 8*x5 + 9,
    # x0^3 + 10*x1^3 + 2*x2^3 + 8*x3^3 - 2*x3^2*x4 + 17*x3*x4^2 + 36*x4^3 + 41*x3^2*x5 + 10*x3*x4*x5 + 8*x4^2*x5 + 49*x3*x5^2 + 38*x4*x5^2 - 24*x5^3 - 19*x3^2 + 20*x3*x4 + 16*x4^2 - 6*x3*x5 - 50*x4*x5 - 43*x5^2 - 6*x3 - 50*x4 + 15*x5 + 8,
    # x3^3 + 14*x4^3 + 39*x5^3 - 29*x6^3 + 9*x6^2*x7 + 13*x6*x7^2 + 10*x7^3 - 26*x6^2*x8 - 19*x6*x7*x8 + 48*x7^2*x8 + 5*x6*x8^2 - 4*x7*x8^2 + 45*x8^3 + 41*x6^2 - 5*x6*x7 + 2*x7^2 - 8*x6*x8 - 34*x7*x8 - 7*x8^2 - 17*x6 - 47*x7 + 46*x8 - 47,
    # x3^3 - 16*x4^3 + 23*x5^3 + 31*x6^3 + 10*x6^2*x7 + 38*x6*x7^2 + x7^3 - 10*x6^2*x8 + 25*x6*x7*x8 - 3*x7^2*x8 + 38*x6*x8^2 + 3*x7*x8^2 - x8^3 - 21*x6^2 + 2*x6*x7 + 24*x7^2 - 2*x6*x8 - 48*x7*x8 + 24*x8^2 + 8*x6 - 10*x7 + 10*x8 - 35,
    # x3^3 + 10*x4^3 + 2*x5^3 + 8*x6^3 - 2*x6^2*x7 + 17*x6*x7^2 + 36*x7^3 + 41*x6^2*x8 + 10*x6*x7*x8 + 8*x7^2*x8 + 49*x6*x8^2 + 38*x7*x8^2 - 24*x8^3 - 26*x6^2 + 38*x6*x7 - 10*x7^2 + 29*x6*x8 + 6*x7*x8 - 11*x8^2 + 45*x6 - 29*x7 + 39*x8 - 3,
    # x6^3 + 14*x7^3 + 39*x8^3 + 45*x9^3 + 16*x9^2 - 40*x9 + 33,
    # x6^3 - 16*x7^3 + 23*x8^3 - x9^3 + 29*x9^2 - 11*x9 - 42,


    ## bar pow bar
    # p = 101
    # -x0 + 6*x1 + x2,
    # -x3 + 6*x4 + x5,
    # x1^17 - 35*x1^16 + 16*x1^15 + 25*x1^14 + 13*x1^13 - 48*x1^12 + 46*x1^11 - 29*x1^10 - 17*x1^9 + 25*x1^8 + 18*x1^7 + 22*x1^6 - 25*x1^5 - 40*x1^4 + 13*x1^3 + 40*x1^2 - 25*x1,
    # -48*x1^16 - 46*x1^15 - 40*x1^14 + 29*x1^13 + 3*x1^12 - 8*x1^11 + 47*x1^10 - 27*x1^9 - 19*x1^8 - 16*x1^7 + 11*x1^6 - 12*x1^5 - 6*x1^4 + 24*x1^3 - 24*x1^2 + 32*x1 - x4,
    # x2^6 - 15*x2^5 - 16*x2^4 - 23*x2^3 - 29*x2^2 - 19*x2,
    # -17*x2^5 - 23*x2^4 + 40*x2^3 - 6*x2^2 + 7*x2 - x5,
    # x3^2 - x6,
    # -x6 + 6*x7 + x8,
    # -x9 + 6*x10 + x11,
    # x7^17 - 35*x7^16 + 16*x7^15 + 25*x7^14 + 13*x7^13 - 48*x7^12 + 46*x7^11 - 29*x7^10 - 17*x7^9 + 25*x7^8 + 18*x7^7 + 22*x7^6 - 25*x7^5 - 40*x7^4 + 13*x7^3 + 40*x7^2 - 25*x7,
    # -48*x7^16 - 46*x7^15 - 40*x7^14 + 29*x7^13 + 3*x7^12 - 8*x7^11 + 47*x7^10 - 27*x7^9 - 19*x7^8 - 16*x7^7 + 11*x7^6 - 12*x7^5 - 6*x7^4 + 24*x7^3 - 24*x7^2 + 32*x7 - x10,
    # x8^6 - 15*x8^5 - 16*x8^4 - 23*x8^3 - 29*x8^2 - 19*x8,
    # -17*x8^5 - 23*x8^4 + 40*x8^3 - 6*x8^2 + 7*x8 - x11,

    # p = 1823
    # -x0 + 38*x1 + x2,
    # -x3 + 38*x4 + x5,
    # x1^48 + 695*x1^47 + 335*x1^46 - 686*x1^45 - 675*x1^44 - 313*x1^43 + 10*x1^42 - 229*x1^41 - 875*x1^40 + 615*x1^39 - 268*x1^38 + 427*x1^37 + 167*x1^36 - 45*x1^35 - 762*x1^34 + 519*x1^33 - 683*x1^32 + 130*x1^31 + 726*x1^30 + 101*x1^29 + 878*x1^28 + 119*x1^27 - 458*x1^26 + 393*x1^25 - 461*x1^24 + 82*x1^23 - 498*x1^22 - 701*x1^21 - 384*x1^20 + 769*x1^19 + 61*x1^18 + 583*x1^17 + 8*x1^16 + 694*x1^15 + 172*x1^14 - 784*x1^13 - 200*x1^12 - 818*x1^11 + 206*x1^10 + 645*x1^9 - 715*x1^8 + 700*x1^7 + 59*x1^6 + 705*x1^5 - 787*x1^4 + 868*x1^3 - 124*x1^2 - 202*x1,
    # 820*x1^47 + 733*x1^46 - 524*x1^45 - 627*x1^44 - 392*x1^43 - 33*x1^42 - 560*x1^41 - 390*x1^40 + 407*x1^39 + 706*x1^38 + 810*x1^37 - 436*x1^36 - 884*x1^35 + 438*x1^34 + 871*x1^33 + 601*x1^32 - 22*x1^31 - 233*x1^30 + 875*x1^29 - 103*x1^28 - 200*x1^27 - 601*x1^26 - 661*x1^25 + 196*x1^24 + 407*x1^23 + 539*x1^22 - 704*x1^21 - 29*x1^20 + 844*x1^19 + 363*x1^18 + 51*x1^17 + 56*x1^16 + 424*x1^15 - 429*x1^14 - 174*x1^13 - 532*x1^12 + 271*x1^11 - 808*x1^10 - 232*x1^9 - 638*x1^8 - 467*x1^7 - 911*x1^6 - 265*x1^5 - 733*x1^4 - 45*x1^3 - 816*x1^2 - 632*x1 - x4 + 29,
    # x2^38 - 703*x2^37 - 496*x2^36 - 250*x2^35 - 867*x2^34 - 885*x2^33 + 627*x2^32 - 773*x2^31 - 878*x2^30 + 571*x2^29 + 796*x2^28 - 190*x2^27 + 78*x2^26 + 683*x2^25 - 535*x2^24 + 688*x2^23 + 33*x2^22 - 623*x2^21 + 572*x2^20 + 164*x2^19 + 659*x2^18 + 49*x2^17 + 44*x2^16 - 288*x2^15 - 769*x2^14 + 596*x2^13 - 664*x2^12 + 269*x2^11 + 698*x2^10 - 349*x2^9 - 246*x2^8 - 637*x2^7 + 406*x2^6 + 330*x2^5 + 52*x2^4 - 486*x2^3 - 816*x2^2 - 507*x2,
    # -808*x2^37 - 315*x2^36 + 322*x2^35 - 71*x2^34 - 645*x2^33 - 902*x2^32 + 523*x2^31 - 838*x2^30 + 607*x2^29 - 115*x2^28 - 703*x2^27 - 122*x2^26 - 366*x2^25 - 202*x2^24 + 684*x2^23 + 17*x2^22 + 854*x2^21 + 449*x2^20 + 98*x2^19 - 785*x2^18 + 606*x2^17 + 534*x2^16 + 490*x2^15 + 587*x2^14 - 32*x2^13 + 88*x2^12 - 847*x2^11 - 210*x2^10 - 112*x2^9 + 138*x2^8 - 33*x2^7 - 429*x2^6 - 165*x2^5 - 89*x2^4 + 34*x2^3 - 803*x2^2 + 715*x2 - x5 + 29,
    # x3^2 - x6,
    # -x6 + 38*x7 + x8,
    # -x9 + 38*x10 + x11,
    # x7^48 + 695*x7^47 + 335*x7^46 - 686*x7^45 - 675*x7^44 - 313*x7^43 + 10*x7^42 - 229*x7^41 - 875*x7^40 + 615*x7^39 - 268*x7^38 + 427*x7^37 + 167*x7^36 - 45*x7^35 - 762*x7^34 + 519*x7^33 - 683*x7^32 + 130*x7^31 + 726*x7^30 + 101*x7^29 + 878*x7^28 + 119*x7^27 - 458*x7^26 + 393*x7^25 - 461*x7^24 + 82*x7^23 - 498*x7^22 - 701*x7^21 - 384*x7^20 + 769*x7^19 + 61*x7^18 + 583*x7^17 + 8*x7^16 + 694*x7^15 + 172*x7^14 - 784*x7^13 - 200*x7^12 - 818*x7^11 + 206*x7^10 + 645*x7^9 - 715*x7^8 + 700*x7^7 + 59*x7^6 + 705*x7^5 - 787*x7^4 + 868*x7^3 - 124*x7^2 - 202*x7,
    # 820*x7^47 + 733*x7^46 - 524*x7^45 - 627*x7^44 - 392*x7^43 - 33*x7^42 - 560*x7^41 - 390*x7^40 + 407*x7^39 + 706*x7^38 + 810*x7^37 - 436*x7^36 - 884*x7^35 + 438*x7^34 + 871*x7^33 + 601*x7^32 - 22*x7^31 - 233*x7^30 + 875*x7^29 - 103*x7^28 - 200*x7^27 - 601*x7^26 - 661*x7^25 + 196*x7^24 + 407*x7^23 + 539*x7^22 - 704*x7^21 - 29*x7^20 + 844*x7^19 + 363*x7^18 + 51*x7^17 + 56*x7^16 + 424*x7^15 - 429*x7^14 - 174*x7^13 - 532*x7^12 + 271*x7^11 - 808*x7^10 - 232*x7^9 - 638*x7^8 - 467*x7^7 - 911*x7^6 - 265*x7^5 - 733*x7^4 - 45*x7^3 - 816*x7^2 - 632*x7 - x10 + 29,
    # x8^38 - 703*x8^37 - 496*x8^36 - 250*x8^35 - 867*x8^34 - 885*x8^33 + 627*x8^32 - 773*x8^31 - 878*x8^30 + 571*x8^29 + 796*x8^28 - 190*x8^27 + 78*x8^26 + 683*x8^25 - 535*x8^24 + 688*x8^23 + 33*x8^22 - 623*x8^21 + 572*x8^20 + 164*x8^19 + 659*x8^18 + 49*x8^17 + 44*x8^16 - 288*x8^15 - 769*x8^14 + 596*x8^13 - 664*x8^12 + 269*x8^11 + 698*x8^10 - 349*x8^9 - 246*x8^8 - 637*x8^7 + 406*x8^6 + 330*x8^5 + 52*x8^4 - 486*x8^3 - 816*x8^2 - 507*x8,
    # -808*x8^37 - 315*x8^36 + 322*x8^35 - 71*x8^34 - 645*x8^33 - 902*x8^32 + 523*x8^31 - 838*x8^30 + 607*x8^29 - 115*x8^28 - 703*x8^27 - 122*x8^26 - 366*x8^25 - 202*x8^24 + 684*x8^23 + 17*x8^22 + 854*x8^21 + 449*x8^20 + 98*x8^19 - 785*x8^18 + 606*x8^17 + 534*x8^16 + 490*x8^15 + 587*x8^14 - 32*x8^13 + 88*x8^12 - 847*x8^11 - 210*x8^10 - 112*x8^9 + 138*x8^8 - 33*x8^7 - 429*x8^6 - 165*x8^5 - 89*x8^4 + 34*x8^3 - 803*x8^2 + 715*x8 - x11 + 29,


    ## Systems from Stegers PhD thesis
    ## From Faugère's paper
    # x0^2*x1-x2^2*x3,
    # x0*x2^2-x1^2*x3,
    # x0*x1^3*x3-x2^4*x3,
    # x2^6*x3-x1^5*x3^2,
    # x1*x2^3-x0^2*x3^2,
    # x1^3*x2*x3-x0^3*x3^2,
    # x2^5*x3-x0^4*x3^2,
    # x1^5*x3^2-x0^4*x2*x3^2,
    # -x0^5*x3^2+x1^2*x2^3*x3^2,
    # x1^6*x3^2-x0*x1^2*x2*x3^4,

    ## Cyclic 3
    # x0+x1+x2,
    # x0*x1+x1*x2+x2*x0,
    # x0*x1*x2-1,

    ## Cyclic 4
    # x0+x1+x2+x3,
    # x0*x1+x1*x2+x2*x3+x3*x0,
    # x0*x1*x2+x1*x2*x3+x2*x3*x0+x3*x0*x1,
    # x0*x1*x2*x3-1,

    ## Cyclic 5
    # x0+x1+x2+x3+x4,
    # x0*x1+x1*x2+x2*x3+x3*x4+x4*x0,
    # x0*x1*x2+x1*x2*x3+x2*x3*x4+x3*x4*x0+x4*x0*x1,
    # x0*x1*x2*x3+x1*x2*x3*x4+x2*x3*x4*x0+x3*x4*x0*x1+x4*x0*x1*x2,
    # x0*x1*x2*x3*x4-1,

    ## cyclic 6
    # x0*x1*x2*x3*x4*x5-1,
    # x0*x1*x2*x3*x4+x0*x1*x2*x3*x5+x0*x1*x2*x4*x5+x0*x1*x3*x4*x5+x0*x2*x3*x4*x5+x1*x2*x3*x4*x5,
    # x0*x1*x2*x3+x0*x1*x2*x5+x0*x1*x4*x5+x0*x3*x4*x5+x1*x2*x3*x4+x2*x3*x4*x5,
    # x0*x1*x2+x0*x1*x5+x0*x4*x5+x1*x2*x3+x2*x3*x4+x3*x4*x5,
    # x0*x1+x0*x5+x1*x2+x2*x3+x3*x4+x4*x5,
    # x0+x1+x2+x3+x4+x5,

    ## cyclic 7
    # x0*x1*x2*x3*x4*x5*x6-1,
    # x0*x1*x2*x3*x4*x5+x0*x1*x2*x3*x4*x6+x0*x1*x2*x3*x5*x6+x0*x1*x2*x4*x5*x6+x0*x1*x3*x4*x5*x6+x0*x2*x3*x4*x5*x6+x1*x2*x3*x4*x5*x6,
    # x0*x1*x2*x3*x4+x0*x1*x2*x3*x6+x0*x1*x2*x5*x6+x0*x1*x4*x5*x6+x0*x3*x4*x5*x6+x1*x2*x3*x4*x5+x2*x3*x4*x5*x6,
    # x0*x1*x2*x3+x0*x1*x2*x6+x0*x1*x5*x6+x0*x4*x5*x6+x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6,
    # x0*x1*x2+x0*x1*x6+x0*x5*x6+x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6,
    # x0*x1+x0*x6+x1*x2+x2*x3+x3*x4+x4*x5+x5*x6,
    # x0+x1+x2+x3+x4+x5+x6,

    ## cyclic 8
    # x0*x1*x2*x3*x4*x5*x6*x7-1,
    # x0*x1*x2*x3*x4*x5*x6+x0*x1*x2*x3*x4*x5*x7+x0*x1*x2*x3*x4*x6*x7+x0*x1*x2*x3*x5*x6*x7+x0*x1*x2*x4*x5*x6*x7+x0*x1*x3*x4*x5*x6*x7+x0*x2*x3*x4*x5*x6*x7+x1*x2*x3*x4*x5*x6*x7,
    # x0*x1*x2*x3*x4*x5+x1*x2*x3*x4*x5*x6+x0*x1*x2*x3*x4*x7+x0*x1*x2*x3*x6*x7+x0*x1*x2*x5*x6*x7+x0*x1*x4*x5*x6*x7+x0*x3*x4*x5*x6*x7+x2*x3*x4*x5*x6*x7,
    # x0*x1*x2*x3*x4+x1*x2*x3*x4*x5+x2*x3*x4*x5*x6+x0*x1*x2*x3*x7+x0*x1*x2*x6*x7+x0*x1*x5*x6*x7+x0*x4*x5*x6*x7+x3*x4*x5*x6*x7,
    # x0*x1*x2*x3+x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6+x0*x1*x2*x7+x0*x1*x6*x7+x0*x5*x6*x7+x4*x5*x6*x7,
    # x0*x1*x2+x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6+x0*x1*x7+x0*x6*x7+x5*x6*x7,
    # x0*x1+x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x0*x7+x6*x7,
    # x0+x1+x2+x3+x4+x5+x6+x7,

    ## Katsura 3
    # x1+2*x2+2*x3+2*x4-x0,
    # x2^2+2*x1*x3+2*x2*x4-x0*x3,
    # 2*x1*x2+2*x2*x3+2*x3*x4-x0*x2,
    # x1^2+2*x2^2+2*x3^2+2*x4^2-x0*x1,

    ## Katsura 5
    # 2*x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+x5^2-x5,
    # x0*x1+x1*x2+2*x2*x3+2*x3*x4+2*x4*x5-x4,
    # 2*x0*x2+2*x1*x3+2*x2*x4+x4^2+2*x3*x5-x3,
    # 2*x0*x3+2*x1*x4+2*x3*x4+2*x2*x5-x2,
    # x3^2+2*x0*x5+2*x1*x5+2*x2*x5-x1,
    # 2*x0+2*x1+2*x2+2*x3+2*x4+x5-1,
]

set_verbose(2)
if get_verbose() >= 2: print(f"––––––––––––\n Input:\n{polys}\n––––––––––––")
gb, voos = f5(polys, homogenize=False)

products = [f5.phi(voo) for voo in voos]
same = [pr == gb_elem for pr, gb_elem in zip(products, gb)]

if get_verbose() >= 0: print(f"––––––––––––\n is Gröbner basis: {Ideal(gb).basis_is_groebner()}")
if get_verbose() >= 1:
    gb_builtin = Ideal(polys).groebner_basis()
    gb_in_gb_builtin = all([b in gb_builtin for b in gb])
    gb_builtin_in_gb = all([b in gb for b in gb_builtin])
    gb_in_ideal_gb_builtin = all([b in Ideal(gb_builtin) for b in gb])
    gb_builtin_in_ideal_gb = all([b in Ideal(gb) for b in gb_builtin])
    span_same_ideal = Ideal(gb) == Ideal(gb_builtin)
    print(f"––––––––––––")
    if gb_in_gb_builtin and gb_builtin_in_gb:
        print(f" \\o/ Correctly computed reduced GB! \\o/")
if get_verbose() >= 2:
    print(f"––––––––––––")
    print(f"              len: {len(gb)}")
    print(f" len sage reduced: {len(gb_builtin)}")
    print(f"gb ⊆ <gb_builtin>: {gb_in_ideal_gb_builtin}")
    print(f"gb_builtin ⊆ <gb>: {gb_builtin_in_ideal_gb}")
    print(f"  span same ideal: {span_same_ideal}")
    print(f"  gb ⊆ gb_builtin: {gb_in_gb_builtin}")
    print(f"  gb_builtin ⊆ gb: {gb_builtin_in_gb}")
    print(f"––––––––––––")
    print(f" GB\\GB_builtin:")
    [print(f"                      {f5.get_sig_from_voo(voos[i])} : {gb[i]}") for i in range(len(gb)) if not gb[i] in gb_builtin]
    print(f" GB∩GB_builtin:")
    [print(f"                      {f5.get_sig_from_voo(voos[i])} : {gb[i]}") for i in range(len(gb)) if gb[i] in gb_builtin]
    print(f" GB_builtin\\GB:")
    [print(f"                      [?, ?] : {b}") for b in gb_builtin if not b in gb]
if get_verbose() >= 2:
    print(f"––––––––––––\n GB:\n{gb}")
    print(f"––––––––––––\n VoOs:\n{voos}")
    print(f"––––––––––––\n Products of Vectors of Origin and input F:")
    col_1_len = max([len(str(voo)) for voo in voos])
    col_2_len = max([len(str(pr)) for pr in products])
    [print(f"{'   ' if same[i] else '[!]'} {str(voos[i]):>{col_1_len}} * F = {str(products[i]):{col_2_len}} {'==' if same[i] else '!='} {gb[i]}") for i in range(len(gb))]
if get_verbose() >= 1 and all(same): print(f"––––––––––––\n \\o/ VoOs are all good! \\o/")
assert all(same), "Some Vector of Origin is misbehaving…"
vods = [] # vectors of destination – how to go from the g_i to the f_j.
for p in polys:
    quotients, rem = polynomial_division(p, gb)
    assert rem == 0, f"p is not in the ideal spanned by GB, or GB is not a Gröbner basis."
    vods += [vector(quotients)]
assert matrix(vods) * vector(gb) == vector(polys), f"Some vector of destination is misbehaving…"
if get_verbose() >= 1:
    len_voo = len(voos[0])
    invoolvement_buckets = [0] * (len_voo + 1)
    bucket_ids_voo = []
    for voo in voos:
        assert any(voo), "VoO is zero-vector: Interreduction of GB has failed."
        invlv_id = sum([1 for x in voo if x]) # number of non-zero entries of VoO
        invoolvement_buckets[invlv_id] += 1
        bucket_ids_voo += [invlv_id]
    print(f"–––––––––––– Involvement of VoOs:\n invoolvement:         {invoolvement_buckets}")
    print(f" mean of invoolvement: {n(mean(bucket_ids_voo), digits=3)}") # weighted mean
    print(f" invoolvement metric:  {n(mean(bucket_ids_voo) / len_voo, digits=3)}") # normalize
    len_vod = len(vods[0])
    invodlvement_buckets = [0] * (len_vod + 1)
    bucket_ids_vod = []
    for vod in vods:
        invlv_id = sum([1 for x in vod if x]) # number of non-zero entries of VoD
        invodlvement_buckets[invlv_id] += 1
        bucket_ids_vod += [invlv_id]
    print(f"–––––––––––– Involvement of VoDs:\n invodlvement:         {invodlvement_buckets}")
    print(f" mean of invodlvement: {n(mean(bucket_ids_vod), digits=3)}") # weighted mean
    print(f" invodlvement metric:  {n(mean(bucket_ids_vod) / len_vod, digits=3)}") # normalize
if get_verbose() >= 1: print(f"––––––––––––\n Zero reductions:     {f5.zero_reductions}")
if get_verbose() >= 1:
    print(f"––––––––––––\n Non-Koszul Syzygies: {'None' if not f5.syzygies else ''}")
    for i in f5.syzygies:
        print(f"                      {f5.voo(i)}")
with open("./voo_file.txt", 'w') as voo_file:
    for voo in voos:
        voo_file.write(f"{voo}\n")
with open("./vod_file.txt", 'w') as vod_file:
    for vod in vods:
        vod_file.write(f"{vod}\n")
