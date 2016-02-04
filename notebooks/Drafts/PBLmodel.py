import math
import cmath


def lambdaCalc(epsilon, stratPar):
    """ function calculate dimensionless function lambda
        epsilon - main fitting parameter of the model
        stratPar - stratification parameter, mu
    """

    if stratPar > 0:
        lamb = (math.sqrt(1 + 40 * epsilon ** 2 * stratPar) - 1) / (
        10 * epsilon * stratPar)
    else:
        tmpLamb = 2 * epsilon
        mistake = 1
        while mistake > 0.01:
            lamb = 2 * epsilon * (
            1 - 16 * epsilon * tmpLamb * stratPar) ** 0.25
            mistake = abs(tmpLamb - lamb) / tmpLamb
            tmpLamb = lamb
    return lamb


def profileFunc_momentum_Calc(lambInv, epsilon, stratPar):
    """calculate profile dimensionless function for momentum
        epsilon - main fitting parameter of the model
        stratPar - stratification parameter, mu
        lambInv = 1 / lamb
    """
    if stratPar > 0:
        return - 5 * epsilon * stratPar / lambInv
    else:
        X = (1 - 16 * epsilon * stratPar / lambInv) ** 0.25
        return 2 * math.log((1 + X) / 2) + \
               math.log((1 + X ** 2) / 2) - 2 * math.atan(X) + math.pi / 2


def profileFunc_heat_Calc(lambInv, epsilon, stratPar):
    """calculate profile dimensionless function for temperature
        epsilon - main fitting parameter of the model
        stratPar - stratification parameter, mu
        lambInv = 1 / lamb
    """
    if stratPar > 0:
        return - 5 * epsilon * stratPar / lambInv
    else:
        X = (1 - 16 * epsilon * stratPar / lambInv) ** 0.25
        return 2 * math.log((1 + X ** 2) / 2)


def GeostrophicWindCalc(windCpx, tempSurf, tempSurf2, gradient, lat, M=2.05):
    """
    calculates geostrophic wind speed and direction
    input:
        windCpx - wind at 10 meters (complex),
        tempAtm - temperature at top of PBL,
        tempSurf - surface temperature,
        lat - latitude
    output:
        geostrWindCpx - geostrophic wind, complex number
    """

    GRAVITYACC = 9.8
    VONKARMANCONST = 0.41
    BETA = GRAVITYACC / 300.
    OMEGA = 7.3e-5
    epsilon = 0.15
    # M = 2.05
    CHARNOCKCONST = 0.018
    corPar = 2. * math.sin(math.radians(lat)) * OMEGA
    roughLength_h = 1e-5
    mistake = 1
    frVelCpx = math.sqrt(0.5e-3) * windCpx
    tmpfrVelCpx = frVelCpx
    while mistake > 0.001:
        roughLength_m = CHARNOCKCONST * abs(frVelCpx) ** 2. / GRAVITYACC
        frVelCpx = VONKARMANCONST * windCpx / (math.log(10. / roughLength_m))
        mistake = abs(tmpfrVelCpx - frVelCpx) / abs(tmpfrVelCpx)
        tmpfrVelCpx = frVelCpx
    it = 0
    stratPar = VONKARMANCONST ** 2. * BETA * (
    tempSurf2 + gradient * 100 - tempSurf) / (corPar * 1.5 * abs(windCpx))
    if abs(stratPar) > 10e-15:
        tmpD = 100
        mistake = 1
        while mistake > 0.01:
            it += 1
            lambInv = 1. / lambdaCalc(epsilon, stratPar)

            C = -2. * lambInv * (M - epsilon) +\
                profileFunc_heat_Calc(lambInv, epsilon, stratPar) - \
                math.log(VONKARMANCONST * epsilon / lambInv)
            EkmanDepth = tmpD
            tempScale = VONKARMANCONST * \
                        (tempSurf2 - tempSurf + gradient * M * EkmanDepth) /\
                        (math.log(abs(frVelCpx) /
                                  (corPar * roughLength_h)) - C)
            stratPar = VONKARMANCONST ** 2 * BETA * tempScale / \
                       (corPar * abs(frVelCpx))
            EkmanDepth = VONKARMANCONST * abs(frVelCpx) / (corPar * lambInv)

            mistake = abs(tmpD - EkmanDepth) / abs(tmpD)
            tmpD = 0.3 * EkmanDepth + 0.7 * tmpD

            if it > 100:
                break

    else:
        stratPar = 0
        lambInv = 1. / lambdaCalc(epsilon, stratPar)

    A = lambInv * cmath.tanh((1 + 1j) * (M - epsilon))
    B = - A + profileFunc_momentum_Calc(lambInv, epsilon, stratPar) - math.log(
        VONKARMANCONST * epsilon / lambInv)
    geostrWindCpx = frVelCpx / VONKARMANCONST * \
                    (math.log(
                        abs(frVelCpx) / (corPar * roughLength_m)) - B - 1j * A)

    return geostrWindCpx, stratPar


