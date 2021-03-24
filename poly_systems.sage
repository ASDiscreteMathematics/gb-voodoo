class PolynomialSystem():
    def __init__(self, field_size, var_name='x', mon_order="degrevlex"):
        self.field = GF(field_size)
        self.var_name = var_name
        self.mon_order = mon_order

    def super_simple(self):
        ring = PolynomialRing(self.field, self.var_name, 2, order=self.mon_order)
        x = ring.gens()
        polys = [x[0] + x[1],
                 x[0]*x[1]]
        return polys

    def perry_intro(self):
        ring = PolynomialRing(self.field, self.var_name, 2, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^2 - 1,
                 x[0] - x[1],
                 x[0]*x[1] - 1,]
        return polys

    def jfs_fglm(self): # Ferdinand's FGLM example
        ring = PolynomialRing(self.field, self.var_name, 2, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^3 - x[0]^2 - x[1]^6 + x[1]^5,
                 x[0]^2*x[1] + x[1]^7 - x[1]^6,
                 x[0]*x[1]^4 + x[1]^6,
                 x[1]^8]
        return polys

    def sys_a(self):
        ring = PolynomialRing(self.field, self.var_name, 3, order=self.mon_order)
        x = ring.gens()
        polys = [x[0] - x[2],
                 x[0]^2 + x[2]^2 - x[1],
                 x[0]^2 + x[1]^2 + x[2]^2 - 1]
        return polys

    def sys_b(self):
        ring = PolynomialRing(self.field, self.var_name, 3, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^5 + x[1]^4 + x[2]^3 - 1,
                 x[0]^3 + x[1]^2 + x[2]^2 - 1]
        return polys

    def sys_c(self):
        ring = PolynomialRing(self.field, self.var_name, 3, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[2] + 3*x[1],
                 x[1] + x[2] + 2,
                 x[0]*x[1] + x[1]^2]
        return polys

    def iva_ex(self): # Has one syzygy vector with signature (x1*x3, 2)
        ring = PolynomialRing(self.field, self.var_name, 4, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[1] - x[1]*x[3],
                 x[0]^2  - x[2]*x[3],
                 x[2]^3  - x[3]^3]
        return polys

    def iva_ex_redundant(self):
        ring = PolynomialRing(self.field, self.var_name, 4, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[1] - x[1]*x[3],
                 x[0]^2  - x[2]*x[3],
                 x[1]*x[2]^3  - x[1]*x[3]^3,
                 x[2]*x[2]^3  - x[2]*x[3]^3,
                 x[2]*x[2]^3  - x[2]*x[3]^3]
        return polys

    def eder_faugere_ex4_3(self):
        ring = PolynomialRing(self.field, self.var_name, 4, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[1]*x[2] - x[2]^2*x[3],
                 x[0]^2*x[1] - x[1]^3,
                 x[1]^3 - x[2]*x[3]^2 - x[3]^3]
        return polys

    def dreg_voo_diff(self):
        ring = PolynomialRing(self.field, self.var_name, 3, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^4*x[1]*x[2] + x[0]^2,
                 x[1]^4*x[2]*x[0] + x[1]^2 + 1,
                 x[2]^4*x[0]*x[1] + x[2]^2]
        return polys

    def dreg_voo_diff_2(self):
        ring = PolynomialRing(self.field, self.var_name, 2, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^100*x[1] + 1,
                 x[0]*x[1]^100]
        return polys

    def dreg_voo_diff_extreme(self):
        ring = PolynomialRing(self.field, self.var_name, 5, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^20*x[1] + x[0]^10,
                 x[1]^20*x[2] + x[1]^10,
                 x[2]^20*x[3] + x[2]^10,
                 x[3]^20*x[4] + x[3]^10,
                 x[4]^20*x[0] + x[4]^10]
        return polys

    def high_dreg(self, degree=7, num_vars=2):
        ring = PolynomialRing(self.field, self.var_name, num_vars, order=self.mon_order)
        x = ring.gens()
        polys = [x[1]^degree + x[0]]
        polys += x[1:]
        return polys

    def low_involv(self):
        ring = PolynomialRing(self.field, self.var_name, 4, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^10*x[1] + x[0],
                 x[1]^10*x[0] + x[1],
                 x[2]^6*x[3] + x[2]^3,
                 x[3]^6*x[2] + x[3]^3,
                 x[0]*x[2]]
                 # x[4]^6*x[5] + x[4]^3,
                 # x[5]^6*x[4] + x[5]^3,
                 # x[2]*x[4],
        return polys

    def faugere_paper(self):
        ring = PolynomialRing(self.field, self.var_name, 4, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]^2*x[1]-x[2]^2*x[3],
                 x[0]*x[2]^2-x[1]^2*x[3],
                 x[0]*x[1]^3*x[3]-x[2]^4*x[3],
                 x[2]^6*x[3]-x[1]^5*x[3]^2,
                 x[1]*x[2]^3-x[0]^2*x[3]^2,
                 x[1]^3*x[2]*x[3]-x[0]^3*x[3]^2,
                 x[2]^5*x[3]-x[0]^4*x[3]^2,
                 x[1]^5*x[3]^2-x[0]^4*x[2]*x[3]^2,
                 -x[0]^5*x[3]^2+x[1]^2*x[2]^3*x[3]^2,
                 x[1]^6*x[3]^2-x[0]*x[1]^2*x[2]*x[3]^4]
        return polys

    def cyclic_3(self):
        ring = PolynomialRing(self.field, self.var_name, 3, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]+x[1]+x[2],
                 x[0]*x[1]+x[1]*x[2]+x[2]*x[0],
                 x[0]*x[1]*x[2]-1]
        return polys

    def cyclic_4(self):
        ring = PolynomialRing(self.field, self.var_name, 4, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]+x[1]+x[2]+x[3],
                 x[0]*x[1]+x[1]*x[2]+x[2]*x[3]+x[3]*x[0],
                 x[0]*x[1]*x[2]+x[1]*x[2]*x[3]+x[2]*x[3]*x[0]+x[3]*x[0]*x[1],
                 x[0]*x[1]*x[2]*x[3]-1]
        return polys

    def cyclic_5(self):
        ring = PolynomialRing(self.field, self.var_name, 5, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]+x[1]+x[2]+x[3]+x[4],
                 x[0]*x[1]+x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[4]*x[0],
                 x[0]*x[1]*x[2]+x[1]*x[2]*x[3]+x[2]*x[3]*x[4]+x[3]*x[4]*x[0]+x[4]*x[0]*x[1],
                 x[0]*x[1]*x[2]*x[3]+x[1]*x[2]*x[3]*x[4]+x[2]*x[3]*x[4]*x[0]+x[3]*x[4]*x[0]*x[1]+x[4]*x[0]*x[1]*x[2],
                 x[0]*x[1]*x[2]*x[3]*x[4]-1]
        return polys

    def cyclic_6(self):
        ring = PolynomialRing(self.field, self.var_name, 6, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[1]*x[2]*x[3]*x[4]*x[5]-1,
                 x[0]*x[1]*x[2]*x[3]*x[4]+x[0]*x[1]*x[2]*x[3]*x[5]+x[0]*x[1]*x[2]*x[4]*x[5]+x[0]*x[1]*x[3]*x[4]*x[5]+x[0]*x[2]*x[3]*x[4]*x[5]+x[1]*x[2]*x[3]*x[4]*x[5],
                 x[0]*x[1]*x[2]*x[3]+x[0]*x[1]*x[2]*x[5]+x[0]*x[1]*x[4]*x[5]+x[0]*x[3]*x[4]*x[5]+x[1]*x[2]*x[3]*x[4]+x[2]*x[3]*x[4]*x[5],
                 x[0]*x[1]*x[2]+x[0]*x[1]*x[5]+x[0]*x[4]*x[5]+x[1]*x[2]*x[3]+x[2]*x[3]*x[4]+x[3]*x[4]*x[5],
                 x[0]*x[1]+x[0]*x[5]+x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[4]*x[5],
                 x[0]+x[1]+x[2]+x[3]+x[4]+x[5]]
        return polys

    def cyclic_7(self):
        ring = PolynomialRing(self.field, self.var_name, 7, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]-1,
                 x[0]*x[1]*x[2]*x[3]*x[4]*x[5]+x[0]*x[1]*x[2]*x[3]*x[4]*x[6]+x[0]*x[1]*x[2]*x[3]*x[5]*x[6]+x[0]*x[1]*x[2]*x[4]*x[5]*x[6]+x[0]*x[1]*x[3]*x[4]*x[5]*x[6]+x[0]*x[2]*x[3]*x[4]*x[5]*x[6]+x[1]*x[2]*x[3]*x[4]*x[5]*x[6],
                 x[0]*x[1]*x[2]*x[3]*x[4]+x[0]*x[1]*x[2]*x[3]*x[6]+x[0]*x[1]*x[2]*x[5]*x[6]+x[0]*x[1]*x[4]*x[5]*x[6]+x[0]*x[3]*x[4]*x[5]*x[6]+x[1]*x[2]*x[3]*x[4]*x[5]+x[2]*x[3]*x[4]*x[5]*x[6],
                 x[0]*x[1]*x[2]*x[3]+x[0]*x[1]*x[2]*x[6]+x[0]*x[1]*x[5]*x[6]+x[0]*x[4]*x[5]*x[6]+x[1]*x[2]*x[3]*x[4]+x[2]*x[3]*x[4]*x[5]+x[3]*x[4]*x[5]*x[6],
                 x[0]*x[1]*x[2]+x[0]*x[1]*x[6]+x[0]*x[5]*x[6]+x[1]*x[2]*x[3]+x[2]*x[3]*x[4]+x[3]*x[4]*x[5]+x[4]*x[5]*x[6],
                 x[0]*x[1]+x[0]*x[6]+x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[4]*x[5]+x[5]*x[6],
                 x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]]
        return polys

    def cyclic_8(self):
        ring = PolynomialRing(self.field, self.var_name, 8, order=self.mon_order)
        x = ring.gens()
        polys = [x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]-1,
                 x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]+x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[7]+x[0]*x[1]*x[2]*x[3]*x[4]*x[6]*x[7]+x[0]*x[1]*x[2]*x[3]*x[5]*x[6]*x[7]+x[0]*x[1]*x[2]*x[4]*x[5]*x[6]*x[7]+x[0]*x[1]*x[3]*x[4]*x[5]*x[6]*x[7]+x[0]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]+x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7],
                 x[0]*x[1]*x[2]*x[3]*x[4]*x[5]+x[1]*x[2]*x[3]*x[4]*x[5]*x[6]+x[0]*x[1]*x[2]*x[3]*x[4]*x[7]+x[0]*x[1]*x[2]*x[3]*x[6]*x[7]+x[0]*x[1]*x[2]*x[5]*x[6]*x[7]+x[0]*x[1]*x[4]*x[5]*x[6]*x[7]+x[0]*x[3]*x[4]*x[5]*x[6]*x[7]+x[2]*x[3]*x[4]*x[5]*x[6]*x[7],
                 x[0]*x[1]*x[2]*x[3]*x[4]+x[1]*x[2]*x[3]*x[4]*x[5]+x[2]*x[3]*x[4]*x[5]*x[6]+x[0]*x[1]*x[2]*x[3]*x[7]+x[0]*x[1]*x[2]*x[6]*x[7]+x[0]*x[1]*x[5]*x[6]*x[7]+x[0]*x[4]*x[5]*x[6]*x[7]+x[3]*x[4]*x[5]*x[6]*x[7],
                 x[0]*x[1]*x[2]*x[3]+x[1]*x[2]*x[3]*x[4]+x[2]*x[3]*x[4]*x[5]+x[3]*x[4]*x[5]*x[6]+x[0]*x[1]*x[2]*x[7]+x[0]*x[1]*x[6]*x[7]+x[0]*x[5]*x[6]*x[7]+x[4]*x[5]*x[6]*x[7],
                 x[0]*x[1]*x[2]+x[1]*x[2]*x[3]+x[2]*x[3]*x[4]+x[3]*x[4]*x[5]+x[4]*x[5]*x[6]+x[0]*x[1]*x[7]+x[0]*x[6]*x[7]+x[5]*x[6]*x[7],
                 x[0]*x[1]+x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[4]*x[5]+x[5]*x[6]+x[0]*x[7]+x[6]*x[7],
                 x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]]
        return polys

    def katsura_3(self):
        ring = PolynomialRing(self.field, self.var_name, 5, order=self.mon_order)
        x = ring.gens()
        polys = [x[1]+2*x[2]+2*x[3]+2*x[4]-x[0],
                 x[2]^2+2*x[1]*x[3]+2*x[2]*x[4]-x[0]*x[3],
                 2*x[1]*x[2]+2*x[2]*x[3]+2*x[3]*x[4]-x[0]*x[2],
                 x[1]^2+2*x[2]^2+2*x[3]^2+2*x[4]^2-x[0]*x[1]]
        return polys

    def katsura_5(self):
        ring = PolynomialRing(self.field, self.var_name, 6, order=self.mon_order)
        x = ring.gens()
        polys = [2*x[0]^2+2*x[1]^2+2*x[2]^2+2*x[3]^2+2*x[4]^2+x[5]^2-x[5],
                 x[0]*x[1]+x[1]*x[2]+2*x[2]*x[3]+2*x[3]*x[4]+2*x[4]*x[5]-x[4],
                 2*x[0]*x[2]+2*x[1]*x[3]+2*x[2]*x[4]+x[4]^2+2*x[3]*x[5]-x[3],
                 2*x[0]*x[3]+2*x[1]*x[4]+2*x[3]*x[4]+2*x[2]*x[5]-x[2],
                 x[3]^2+2*x[0]*x[5]+2*x[1]*x[5]+2*x[2]*x[5]-x[1],
                 2*x[0]+2*x[1]+2*x[2]+2*x[3]+2*x[4]+x[5]-1]
        return polys

    def dreg_is_off(self):
        ring = PolynomialRing(self.field, self.var_name, 3, order=self.mon_order)
        x = ring.gens()
        polys = [x[1]^4 - x[0]*x[2]^3,
                 x[0]^2*x[1]^2 + x[1]^3*x[2],
                 x[0]^3*x[1] - x[0]^2*x[1]^2 + 1]
        return polys

    def random(self, num_vars=4, degree=4, terms=6):
        ring = PolynomialRing(self.field, self.var_name, num_vars, order=self.mon_order)
        polys = [ring.random_element(degree, terms) for _ in range(num_vars)]
        return polys
