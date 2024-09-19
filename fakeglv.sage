## utils
def halfgcd(a, b):
    half_deg = (Integer(a).nbits()/2.).ceil()
    r = (a)
    r_ = (b)
    u = (1)
    v = (0)
    u_ = (0)
    v_ = (1)
    while Integer(r).nbits() > half_deg:
        q, rem = r.quo_rem(r_)
        (r, u, v, r_, u_, v_) = (r_, u_, v_, r-q*r_, u - q*u_, v - q*v_)
    assert r == a*u + b*v
    assert r_ == a*u_ + b*v_
    return r, u

## test
p256 = {
    "name": "P-256", # NIST curve
    "curve": EllipticCurve(
        GF(0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff),
        [0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc, 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b],
    ),
    "r": 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551,
}
frp256v1 = {
    "name": "FRP256v1", # ANSSI curve
    "curve": EllipticCurve(
        GF(0xf1fd178c0b3ad58f10126de8ce42435b3961adbcabc8ca6de8fcf353d86e9c03),
        [0xf1fd178c0b3ad58f10126de8ce42435b3961adbcabc8ca6de8fcf353d86e9c00, 0xee353fca5428a9300d4aba754a44c00fdfec0c9ae4b1a1803075ed967b7bb73f],
    ),
    "r": 0xf1fd178c0b3ad58f10126de8ce42435b53dc67e140d2bf941ffdd459c6d655e1,
}
curves = [p256, frp256v1]

for curve in curves:
    print("\n*** Test on", curve["name"])
    E = curve["curve"]
    r = curve["r"]

    for i in range(10):
        Fr = GF(r)
        s = Integer(Fr.random_element())
        u, v = halfgcd(s, r)
        assert(Fr(v*s) == Fr(u))
        nbits = (r.nbits()/2.).ceil()+1
        assert(v.nbits() <= nbits and u.nbits() <= nbits)
        P = E.random_point() * (E.order() // r)
        assert(r*P == E(0))
        Q = s*P
        assert(s*P - Q == E(0))
        assert(v*s*P - v*Q == E(0))
        assert(u*P - v*Q == E(0))
        print("Checking [{}]P == Q is equivalent to checking:".format(s))
        print("[{}]P - [{}]Q == {}".format(u, v, E(0)))
