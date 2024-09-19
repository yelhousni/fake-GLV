## utils
def quo_rem_K(a, b, K):
    """
     Quotient and remainder in Eisenstein ring of integers K
    :param a: element in K
    :param b: element in K
    :return q, r: where a = q*b + r
    """
    num = a * b.conjugate()
    denum = b.norm()
    q1, r1 = Integer(num[0]).quo_rem(Integer(denum))
    q2, r2 = Integer(num[1]).quo_rem(Integer(denum))
    return K(q1+j*q2), K(r1+j*r2)/b.conjugate()

def halfgcd_K(a, b, K):
    """
     Rational reconstruction a.k.a half-GCD
    :param a: element in K
    :param b: element in K
    :return u, v : where u = b*v mod a
    """
    a_run, b_run = a, b
    u, v = 1, 0
    u_, v_ = v, u
    nbits = (Integer(a.norm()).nbits()/2.).ceil()
    while Integer(a_run.norm()).nbits() > nbits and Integer(b_run.norm()).nbits() > nbits:
        q, rem = quo_rem_K(a_run, b_run, K)
        (a_run, u, v, b_run, u_, v_) = (b_run, u_, v_, rem, u - q*u_, v - q*v_)
    assert a_run == a*u + b*v
    assert b_run == a*u_ + b*v_
    return a_run, v

def ee_basis(n, lam):
    """
     Implementation of extended euclidean as a generator of a sequence of r,s,t
    :param n: integer
    :param lam: integer
    :return r,s,t: where s*n+t*lam = r
    """
    r0, r1, r2 = n, lam, None
    s0, s1, s2 = 1, 0, None
    t0, t1, t2 = 0, 1, None
    seq = []
    while True:
        r2, q1 = r0 % r1, r0 // r1
        s2 = s0 - q1 * s1
        t2 = t0 - q1 * t1
        yield r1, s1, t1
        s1, s0 = s2, s1
        t1, t0 = t2, t1
        r1, r0 = r2, r1


def find_basis(k, n, lam):
    """
    :param k: integer
    :param n: modulus
    :param lam: fast scalar
    :return v1,v2: basis satisfying f(v1)=f(v2)=0
    where f(i,j) = i+lam*j (mod n)
    """
    r0, t0, r1, t1 = n, 0, lam, 1
    ee_gen = ee_basis(n, lam)
    for r, s, t in ee_gen:
        if r < sqrt(n):
            break
        r0, s0, t0 = r, s, t
    r1, s1, t1 = next(ee_gen)
    if r0**2 + t0**2 > r1**2 + t1**2:
        return vector((r, -t)), vector((r1, -t1))
    return vector((r, -t)), vector((r0, -t0))


def find_closest(k, vs: list) -> vector:
    """
    :param k: integer
    :param vs: list of vectors
    :return: a closest vector to k in the span of vs
    """
    m = [[QQ(i) for i in v] for v in vs]
    m = matrix(m)
    m = m.transpose()
    sol = m.solve_right(vector([QQ(k)] + [0] * (len(vs) - 1)))
    bs = list(map(lambda x: x.round(), sol))
    return sum(map(lambda x: x[0] * x[1], zip(bs, vs)))


def glv_decompose_simple(k, lam, n) -> vector:
    """
    Uses euclidean extended alg for glv decomposition k=k1+k2*lam (mod n)
    :param k: integer
    :param lam: fast scalar
    :param n: modulus
    :return: vector (k1,k2)
    """
    v1, v2 = find_basis(k, n, lam)
    # assert (v1[0]+v1[1]*lam)%n==0
    # assert (v2[0]+v2[1]*lam)%n==0
    v = find_closest(k, [v1, v2])
    # assert (v[0]+v[1]*lam)%n==0
    short = vector([k, 0]) - v
    assert (short[0] + short[1] * lam) % n == k % n
    return short

def phi0(E, P):
    """
    j=0 endomorphism
    :param E: Elliptic curve
    :param P: point on E
    :return: point lam*P
    """
    return E(P[0]*w, P[1])

## test
secp256k1 = {
    "name": "secp256k1", # Ethereum ecdsa curve
    "curve": EllipticCurve(
        GF(115792089237316195423570985008687907853269984665640564039457584007908834671663),
        [0, 7],
    ),
    "r": 115792089237316195423570985008687907852837564279074904382605163141518161494337,
    "r0": 64502973549206556628585045361533709077,
    "r1": -303414439467246543595250775667605759171,
    "w": 55594575648329892869085402983802832744385952214688224221778511981742606582254,
    "l": 37718080363155996902926221483475020450927657555482586988616620542887997980018,
    "ring": EisensteinIntegers(),
}
bn254 = {
    "name": "bn254", # Ethereum ec-precompiles curve
    "curve": EllipticCurve(
        GF(21888242871839275222246405745257275088696311157297823662689037894645226208583),
        [0, 3],
    ),
    "r": 21888242871839275222246405745257275088548364400416034343698204186575808495617,
    "r0": 9931322734385697763,
    "r1": -147946756881789319000765030803803410728,
    "w": 2203960485148121921418603742825762020974279258880205651966,
    "l": 4407920970296243842393367215006156084916469457145843978461,
    "ring": EisensteinIntegers(),
}
bls12_381 = {
    "name": "bls12-381", # Ethereum bls-sig curve
    "curve": EllipticCurve(
        GF(4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787),
        [0, 4],
    ),
    "r": 52435875175126190479447740508185965837690552500527637822603658699938581184513,
    "r0": 228988810152649578064853576960394133503,
    "r1": -1,
    "w": 4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436,
    "l": 228988810152649578064853576960394133503,
    "ring": EisensteinIntegers(),
}
bls12_377 = {
    "name": "bls12-377", # Linea aggregation curve
    "curve": EllipticCurve(
        GF(258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177),
        [0, 1],
    ),
    "r": 8444461749428370424248824938781546531375899335154063827935233455917409239041,
    "r0": 91893752504881257701523279626832445440,
    "r1": -1,
    "w": 80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945,
    "l": 91893752504881257701523279626832445440,
    "ring": EisensteinIntegers(),
}
bw6_761 = {
    "name": "bw6-761", # Linea aggregation curve
    "curve": EllipticCurve(
        GF(6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068299),
        [0, -1],
    ),
    "r": 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177,
    "r0": 293634935485640680722085584138834120315328839056164388863,
    "r1": -293634935485640680722085584138834120324914961969255022593,
    "w": 1968985824090209297278610739700577151397666382303825728450741611566800370218827257750865013421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868639090213645505120388400344940985710520836292650,
    "l": 80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945,
    "ring": EisensteinIntegers(),
}
curves = [secp256k1, bn254, bls12_381, bls12_377, bw6_761]

for curve in curves:
    print("\n*** Test on", curve["name"])
    E = curve["curve"]
    K.<j> = curve["ring"]
    r = curve["r"]
    l = curve["l"]
    w = curve["w"]
    r0, r1 = curve["r0"], curve["r1"]
    R=K(Integer(r0)+j*Integer(r1))
    assert((r0+r1*l)%r == 0)

    for i in range(10):
        s = Integer(GF(r).random_element())
        s0, s1 = glv_decompose_simple(s, l, r)
        assert((s0 + s1*l)%r == s)
        S = K(Integer(s0)+j*Integer(s1))
        U, V = halfgcd_K(R, S, K)
        assert((S*V-U).mod(R)==0)
        v0, v1 = Integer(V[0]), Integer(V[1])
        u0, u1 = Integer(U[0]), Integer(U[1])
        fourth_nbits = (Integer(R.norm())/4.).ceil().nbits()
        assert(u0.nbits() < fourth_nbits and u1.nbits() < fourth_nbits and v0.nbits() < fourth_nbits and v0.nbits() < fourth_nbits)
        P = E.random_point() * (E.order() // r)
        assert(r*P == E(0))
        Q = s*P
        assert(s*P - Q == E(0))
        assert((s0+l*s1)*P - Q == E(0))
        assert((v0+l*v1)*(s0+l*s1)*P - (v0+l*v1)*Q == E(0))
        assert((u0+l*u1)*P - (v0+l*v1)*Q == E(0))
        assert(u0*P + u1*phi0(E, P) - v0*Q - v1*phi0(E, Q) == E(0))
        print("Checking [{}]P == Q is equivalent to checking:".format(s))
        print("[{}]P + [{}]phi(P) - [{}]Q - [{}]phi(Q) == {}".format(u0, u1, v0, v1, E(0)))
