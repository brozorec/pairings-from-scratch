from sage.all import *

q = 13
a = 8
b = 8
F13 = GF(q)
# Elliptic Curve defined by y^2 = x^3 + 8*x + 8 over Finite Field of size 13
E_13_1 = EllipticCurve(F13, [a, b])
print(E_13_1)

# Curve order: 20 (factors in 5 * 2 * 2), embedding degree k wrt 5 == 4
order = E_13_1.order()
r = 5
k = 4

t = polygen(ZZ, 't')
F13t = F13['t']

# P_MOD = F13t(t**4 - 4*(t**2) + 5)
P_MOD = F13t(t**4 + 2)
F13_4 = GF(q**k, name='t', modulus=P_MOD)
E_13_4 = EllipticCurve(F13_4, [a, b])
# gen = E_13_4.gens()[1]
print(E_13_1.order())
# print(E_13_4)
E_13_4_5 = E_13_4(0).division_points(r)

def trace_map(point, q, k):
    E = point.curve()
    [x, y, z] = point
    for i in range(1, k):
        power = q ** i
        point += E(x ** power, y ** power, z ** power)

    return point

def find_g1(full_torsion, q, k):
    cardinality = full_torsion.cardinality()
    L_G1 = []
    for i in range(cardinality):
        point = full_torsion[i]
        if trace_map(point, q, k) == point * k:
            L_G1.append(point)

    return L_G1
        
def find_g2(full_torsion, q, k):
    cardinality = full_torsion.cardinality()
    L_G2 = []
    for i in range(cardinality):
        point = full_torsion[i]
        if trace_map(point, q, k) == point.curve()(0, 1, 0):
            L_G2.append(point)

    return L_G2

def linefunc(P, Q, T):
    x1, y1, _ = P
    x2, y2, _ = Q
    xt, yt, _ = T

    if x1 != x2:
        slope = (y2 - y1) / (x2 - x1)
        # print("slope add: ", slope)
        return slope * (xt - x1) - (yt - y1)
    elif y1 == y2:
        slope = (3 * x1**2 + a) / (2 * y1)
        # print("slope double: ", slope)
        return slope * (xt - x1) - (yt - y1)
    else:
        # print("slope none")
        # print(x1, x2, y1, y2)
        # print(xt - x1)
        # print(xt - x1 - yt + y1)
        return xt - x1

def pairing(P, Q):
    R = P
    f = F13_4(1)
    loop_bin = bin(r)[:1:-1] # bin(r) = '0b101' => cut and reverse =>'101'
    for i in range(1, len(loop_bin)):
        print(loop_bin[i])
        f_prim = linefunc(R, R, Q)
        # print("order of f_prim: ", f_prim.multiplicative_order())
        print("f_prim: ", f_prim)
        f = f * f * f_prim
        print("f_acc: ", f)
        R = 2 * R
        print("A: ", R)
        if loop_bin[i] == '1':
            f_prim = linefunc(R, P, Q)
            print("f_prim: ", f_prim)
            f = f * f_prim
            print("f_acc: ", f)
            R = R + P
            print("A: ", R)

    assert R == P * r
    
    # P1 = E_13_4(P[0] ** q, P[1] ** q)
    # nP2 = E_13_4(P1[0] ** q, -P1[1] ** q)
    # f = f * linefunc(R, P1, Q)
    # R = R + P1
    # f = f * linefunc(R, nP2, Q)

    final_exp = (q**k - 1) // r
    print("final_exp: ", final_exp)
    print("final f: ", f)
    return f ** final_exp
        

r_E_13_4_5 = Set(E_13_4_5) # 25 points
# print(len(r_E_13_4_5))

G1 = find_g1(r_E_13_4_5, q, k)
G2 = find_g2(r_E_13_4_5, q, k)
# print(G1)
# print(G2)

# P = G1[2]
P = E_13_4(8, 8)
# P = E_13_4(7, 2)
print("P: ", P)
# print("P trace_map: ", trace_map(P, q, k))
# print("2*P: ", 2 * P)
# Q = G2[3]
Q = E_13_4(4*t**2 + 7, 5*t**3 + 10*t)
# Q = E_13_4(9*t**2 + 7,12*t**3 + 2*t)
print("Q: ", Q)

print("------")
p = pairing(P, Q)
# print("order of pairing value: ", p ** p.multiplicative_order())
print(p)
print(P.tate_pairing(Q, Integer(5), Integer(4)))
print(P.weil_pairing(Q, 5))













