load('f4_5.sage')
load('common.sage')

def print_gb_analytics(polys, verbosity=2, write_to_disk=True, return_gb=False):
    f5 = F5() # or: FR5, F5C, F4F5
    gb, voos = f5(polys, homogenize=False)

    products = [f5.phi(voo) for voo in voos]
    same = [pr == gb_elem for pr, gb_elem in zip(products, gb)]

    if verbosity >= 2: print(f"––––––––––––\n Input:\n{f5.F}\n––––––––––––")
    if verbosity >= 0:
        deg_macaulay = 1 + sum([p.degree() - 1 for p in f5.F])
        deg_max = max([b.degree() for b in gb])
        print(f"––––––––––––")
        print(f" is Gröbner basis:       {Ideal(gb).basis_is_groebner()}")
        print(f" degree of regularity:   {f5.dreg}")
        print(f" macaulay bound:         {deg_macaulay}")
        print(f" highest degree in gb:   {deg_max}")
    if verbosity >= 3:
        gb_builtin = Ideal(polys).groebner_basis()
        gb_in_gb_builtin = all([b in gb_builtin for b in gb])
        gb_builtin_in_gb = all([b in gb for b in gb_builtin])
        if gb_in_gb_builtin and gb_builtin_in_gb:
            print(f"––––––––––––\n  \\o/ Correctly computed reduced GB! \\o/")
        else:
            print(f"––––––––––––\n  :-( Bugs happend when computing GB :-(")
        gb_in_ideal_gb_builtin = all([b in Ideal(gb_builtin) for b in gb])
        gb_builtin_in_ideal_gb = all([b in Ideal(gb) for b in gb_builtin])
        span_same_ideal = Ideal(gb) == Ideal(gb_builtin)
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
    if verbosity >= 2:
        print(f"–––––––––––– More details about the VoO's.")
        print(f" voo degrees – signature – gb_elem lead term")
        for voo, b in zip(voos, gb):
            print(f"[", end="")
            for v in voo:
                print(f"{v.degree() if v.degree() >= 0 else '   ':>3}", end=" ")
            print(f"] – {f5.get_sig_from_voo(voo)} – {b.lt()}")
    if verbosity >= 3:
        print(f"––––––––––––\n GB:\n{gb}")
        print(f"––––––––––––\n VoOs:\n{voos}")
        print(f"––––––––––––\n Products of Vectors of Origin and input F:")
        col_1_len = max([len(str(voo)) for voo in voos])
        col_2_len = max([len(str(pr)) for pr in products])
        [print(f"{'   ' if same[i] else '[!]'} {str(voos[i]):>{col_1_len}} * F = {str(products[i]):{col_2_len}} {'==' if same[i] else '!='} {gb[i]}") for i in range(len(gb))]
    if verbosity >= 1 and all(same): print(f"––––––––––––\n \\o/ VoOs are all good! \\o/")
    assert all(same), "Some Vector of Origin is misbehaving…"
    vods = [] # vectors of destination – how to go from the g_i to the f_j.
    for p in f5.F:
        quotients, rem = polynomial_division(p, gb)
        vods += [vector(quotients)]
        assert not rem, f"p is not in the ideal spanned by GB, or GB is not a Gröbner basis."
    assert matrix(vods) * vector(gb) == vector(f5.F), f"Some vector of destination is misbehaving…"
    ## Involvement metric. Subject to change.
    if verbosity >= 1:
        len_voo = len(voos[0])
        invoolvement_buckets = [0] * len_voo
        bucket_ids_voo = []
        for voo in voos:
            assert any(voo), "VoO is zero-vector: Interreduction of GB has failed."
            invlv_id = sum([1 for x in voo if x]) - 1
            invoolvement_buckets[invlv_id] += 1
            bucket_ids_voo += [invlv_id]
        print(f"–––––––––––– Involvement of VoOs:\n invoolvement:         {invoolvement_buckets}")
        print(f" mean of invoolvement: {n(mean(bucket_ids_voo), digits=3)}") # weighted mean
        print(f" invoolvement metric:  ", end='')
        print( n(mean(bucket_ids_voo) / (len_voo - 1), digits=3) if len_voo > 1 else 'NaN') # normalize
        len_vod = len(vods[0])
        invodlvement_buckets = [0] * len_vod
        bucket_ids_vod = []
        for vod in vods:
            invlv_id = sum([1 for x in vod if x]) - 1
            invodlvement_buckets[invlv_id] += 1
            bucket_ids_vod += [invlv_id]
        print(f"–––––––––––– Involvement of VoDs:\n invodlvement:         {invodlvement_buckets}")
        print(f" mean of invodlvement: {n(mean(bucket_ids_vod), digits=3)}") # weighted mean
        print(f" invodlvement metric:  ", end='')
        print( 1-n(mean(bucket_ids_vod) / (len_vod - 1), digits=3) if len_vod > 1 else 'NaN') # normalize
    if verbosity >= 1: print(f"––––––––––––\n Zero reductions:     {f5.zero_reductions}")
    if verbosity >= 1:
        print(f"––––––––––––\n Non-Koszul Syzygies: {'None' if not f5.syzygies else ''}")
        for i in f5.syzygies:
            print(f"                      {f5.voo(i)}")
    if write_to_disk:
        with open("./voo_file.txt", 'w') as voo_file:
            for voo in voos:
                voo_file.write(f"{voo}\n")
        with open("./vod_file.txt", 'w') as vod_file:
            for vod in vods:
                vod_file.write(f"{vod}\n")
    if return_gb:
        return gb
