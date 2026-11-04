#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
shs_mie.py  --  Diffusion de seconde harmonique (SHS) d'une sphere : theorie de Mie non lineaire
=====================================================================================================

Modele
------
* Sphere isotrope (rayon a, indice complexe m1 a w et m2 a 2w) dans un milieu d'indice n_m(w), n_m(2w).
* Onde plane incidente a la pulsation w, polarisee suivant x, se propageant suivant +z (e^{-i w t}).
* Champ fondamental : theorie de Mie lineaire (Bohren & Huffman, chap. 4, notations a_n, b_n, c_n, d_n,
  pi_n, tau_n, M_{o1n}, N_{e1n}).
* Source non lineaire : nappe de polarisation de SURFACE a 2w en r = a (sphere centrosymetrique -> pas de
  terme dipolaire electrique de volume ; les termes de volume "gradient/magnetique" sont negliges) :

      p_r     = chi_rrr * E_r^2 + chi_rtt * (E_theta^2 + E_phi^2)
      p_theta = 2 * chi_trt * E_r * E_theta
      p_phi   = 2 * chi_trt * E_r * E_phi            (p = P_surf / eps0)

  (chi_rrr = chi_perp-perp-perp, chi_rtt = chi_perp-par-par, chi_trt = chi_par-perp-par)

* Polarisation incidente variable : voir shs_vs_gamma (courbes I_V(gamma), I_H(gamma)).

Limites du modele / de l'interface en ligne de commande
-------------------------------------------------------
* Surface isotrope ACHIRALE : trois constantes seulement (chi_rrr, chi_rtt, chi_trt) ; les autres
  composantes du tenseur sont ignorees (un avertissement est affiche).
* L'interface en ligne de commande / GUI utilise le meme indice a w et a 2w (pas de dispersion) et
  ignore la partie imaginaire de l'indice du milieu ; shs_sphere() accepte des indices distincts.
* Les intensites sont r^2|E_2w|^2/|E_w|^4 (prefacteur radiatif inclus, voir shs_units.py).

* Probleme a 2w : on developpe la nappe sur les harmoniques spheriques vectorielles (Y_nm, B_nm, C_nm),
  puis on resout, pour chaque (n,m), deux systemes 2x2 (TM et TE) donnant les coefficients rayonnes.
  Conditions de saut a r = a (convention e^{-i w t}) :

      E_t(out) - E_t(in)             = - grad_t(p_r) / eps_sheet
      r_hat x [H(out) - H(in)]       = - i w P_t           ( -> -i k0 p_t avec H~ = Z0 H )

Regles de selection : le champ fondamental n'a que des composantes m = +-1  =>  la source SH n'a que
m = 0, +-2.

Unites : longueurs en nm (ou toute unite coherente), E0 = 1, chi sans dimension "effective".
Les intensites sont donc en unites arbitraires (|chi E0^2|^2).
"""

import math
import numpy as np

import shs_units
import shs_angles
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn, spherical_yn, lpmv, gammaln


# ----------------------------------------------------------------------------------------------
#  Fonctions de Bessel spheriques utiles
# ----------------------------------------------------------------------------------------------
def _jh(n, z):
    """
    Retourne j_n(z), h_n^(1)(z), D_j(z), D_h(z) avec D_z = [z f_n(z)]' / z = (f_n + z f_n') / z
    (z peut etre complexe).
    """
    j = spherical_jn(n, z)
    jp = spherical_jn(n, z, derivative=True)
    y = spherical_yn(n, z)
    yp = spherical_yn(n, z, derivative=True)
    h = j + 1j * y
    hp = jp + 1j * yp
    return j, h, (j + z * jp) / z, (h + z * hp) / z


def size_cutoff(x):
    """Critere de Wiscombe pour le nombre d'harmoniques."""
    x = abs(x)
    return int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))


# ----------------------------------------------------------------------------------------------
#  Theorie de Mie lineaire (Bohren & Huffman, mu = mu1 = 1)
# ----------------------------------------------------------------------------------------------
def mie_coefficients(n, x, m):
    """a_n, b_n (champ diffuse) et c_n, d_n (champ interne). x = k a (milieu), m = m_sphere / n_milieu."""
    jx, hx, Djx, Dhx = _jh(n, x)
    jm, _, Djm, _ = _jh(n, m * x)
    psi_x, dpsi_x = x * jx, x * Djx          # Riccati-Bessel psi_n et sa derivee
    xi_x, dxi_x = x * hx, x * Dhx            # xi_n = x h_n^(1) et sa derivee
    psi_m, dpsi_m = m * x * jm, m * x * Djm

    den_a = m * psi_m * dxi_x - xi_x * dpsi_m
    den_b = psi_m * dxi_x - m * xi_x * dpsi_m
    W = psi_x * dxi_x - xi_x * dpsi_x       # Wronskien

    an = (m * psi_m * dpsi_x - psi_x * dpsi_m) / den_a
    bn = (psi_m * dpsi_x - m * psi_x * dpsi_m) / den_b
    cn = m * W / den_b
    dn = m * W / den_a
    return an, bn, cn, dn


def pi_tau(nmax, mu):
    """pi_n(mu), tau_n(mu) par recurrence (BH eq. 4.47). Tableaux de forme (nmax+1, len(mu))."""
    mu = np.asarray(mu, dtype=float)
    pi = np.zeros((nmax + 1, mu.size))
    tau = np.zeros((nmax + 1, mu.size))
    if nmax >= 1:
        pi[1] = 1.0
        tau[1] = mu
    for n in range(2, nmax + 1):
        pi[n] = (2 * n - 1) / (n - 1) * mu * pi[n - 1] - n / (n - 1) * pi[n - 2]
        tau[n] = n * mu * pi[n] - (n + 1) * pi[n - 1]
    return pi, tau


def fundamental_fields(theta, phi, x, m, nmax, E0=1.0):
    """
    Champ fondamental TOTAL en r = a, juste a l'exterieur et juste a l'interieur de la sphere.
    theta : (nt,), phi : (np,).  Retour : dict 'out', 'in' -> (E_r, E_theta, E_phi), formes (nt, np).
    """
    theta = np.asarray(theta, float)
    phi = np.asarray(phi, float)
    st = np.sin(theta)[:, None]
    cp = np.cos(phi)[None, :]
    sp = np.sin(phi)[None, :]
    pi, tau = pi_tau(nmax, np.cos(theta))

    shape = (theta.size, phi.size)
    out = [np.zeros(shape, complex) for _ in range(3)]
    ins = [np.zeros(shape, complex) for _ in range(3)]

    for n in range(1, nmax + 1):
        En = E0 * 1j ** n * (2 * n + 1) / (n * (n + 1))
        an, bn, cn, dn = mie_coefficients(n, x, m)
        pin, taun = pi[n][:, None], tau[n][:, None]

        j, h, Dj, Dh = _jh(n, x)
        jm, _, Djm, _ = _jh(n, m * x)

        # (zM, zN, DN, rho) : combinaisons lineaires de M_{o1n} et N_{e1n} (BH eq. 4.37-4.38)
        sets = (
            (out, j - bn * h, -1j * j + 1j * an * h, -1j * Dj + 1j * an * Dh, x),
            (ins, cn * jm, -1j * dn * jm, -1j * dn * Djm, m * x),
        )
        for target, zM, zN, DN, rho in sets:
            target[0] += En * cp * n * (n + 1) * st * pin * zN / rho
            target[1] += En * (cp * pin * zM + cp * taun * DN)
            target[2] += En * (-sp * taun * zM - sp * pin * DN)

    return {"out": tuple(out), "in": tuple(ins)}


# ----------------------------------------------------------------------------------------------
#  Harmoniques spheriques vectorielles : Y_nm, B_nm = r grad Y, C_nm = grad Y x r_hat ... (composantes)
# ----------------------------------------------------------------------------------------------
def vsh_angular(n, m, theta, phi):
    """
    Retourne Y, (B_theta, B_phi), (C_theta, C_phi) avec
        B = theta_hat dY/dtheta + phi_hat (1/sin) dY/dphi
        C = theta_hat (1/sin) dY/dphi - phi_hat dY/dtheta        (= -r_hat x B)
    theta : (nt,1) ou (nt,), phi : (1,np) ou (np,).
    """
    theta = np.asarray(theta, float).reshape(-1, 1)
    phi = np.asarray(phi, float).reshape(1, -1)
    x = np.cos(theta)
    sint = np.sin(theta)
    sint = np.where(np.abs(sint) < 1e-12, 1e-12, sint)

    ma = abs(m)
    P = lpmv(ma, n, x)
    Pp = lpmv(ma + 1, n, x) if ma + 1 <= n else np.zeros_like(x)
    if ma >= 1:
        Pm = lpmv(ma - 1, n, x)
    else:                                     # P_n^{-1} = -P_n^1 / (n (n+1))
        Pm = -lpmv(1, n, x) / (n * (n + 1.0))
    dP = 0.5 * (Pp - (n + ma) * (n - ma + 1.0) * Pm)   # dP_n^m(cos th)/d th (phase de Condon-Shortley)

    norm = np.sqrt((2 * n + 1) / (4 * np.pi) * np.exp(gammaln(n - ma + 1.0) - gammaln(n + ma + 1.0)))
    s = (-1.0) ** ma if m < 0 else 1.0       # Y_{n,-m} = (-1)^m conj(Y_{n,m})
    eim = np.exp(1j * m * phi)

    Y = s * norm * P * eim
    dY = s * norm * dP * eim
    mY = 1j * m * Y / sint                   # (1/sin th) dY/dphi
    return Y, (dY, mY), (mY, -dY)


# ----------------------------------------------------------------------------------------------
#  Resolution SH
# ----------------------------------------------------------------------------------------------
def shs_sphere(radius, wavelength, m_sphere_w, m_sphere_2w,
               n_med_w=1.33, n_med_2w=1.33,
               chi_rrr=1.0, chi_rtt=1.0, chi_trt=1.0,
               eps_sheet=1.0, side="outside",
               theta_det=None, nmax_f=None, nmax_sh=None, nphi=16,
               extra_f=0, extra_sh=0, m_values=(-2, 0, 2)):
    """
    Calcule la diffusion de seconde harmonique d'une sphere.

    Parametres
    ----------
    radius, wavelength : rayon de la sphere et longueur d'onde DANS LE VIDE du fondamental (meme unite)
    m_sphere_w, m_sphere_2w : indices (complexes) de la sphere a w et a 2w
    n_med_w, n_med_2w : indices du milieu environnant a w et 2w
    chi_rrr, chi_rtt, chi_trt : composantes du tenseur chi^(2) de surface (isotrope)
    eps_sheet : permittivite lineaire de la nappe (1 par defaut ; ambiguite connue du modele)
    side : 'outside' | 'inside' | 'mean' : cote de l'interface ou l'on evalue le champ fondamental
    theta_det : angles de diffusion (rad) dans ]0, pi[ ; defaut 0.5 deg ... 179.5 deg
    nmax_f, nmax_sh : nombre de multipoles (n) du fondamental et de la SH. Par defaut : critere de Wiscombe.
    extra_f, extra_sh : multipoles SUPPLEMENTAIRES ajoutes au critere de Wiscombe (ou a nmax_f / nmax_sh
                        s'ils sont fournis).
    m_values : ordres azimutaux m inclus dans le developpement SH. Par regle de selection, seuls
               m = 0, +-2 sont non nuls pour une excitation lineaire ; on peut elargir (ex. range(-6, 7))
               pour verifier numeriquement que les autres ordres sont nuls.

    Retour : dict (coefficients rayonnes, diagrammes angulaires dans les plans phi=0 et phi=90 deg,
                   puissance SH totale (u.a.))
    """
    a = float(radius)
    k1 = 2 * np.pi * n_med_w / wavelength           # nombre d'onde fondamental dans le milieu
    x1 = k1 * a
    m1 = m_sphere_w / n_med_w

    lam2 = wavelength / 2.0
    k0 = 2 * np.pi / lam2                           # nombre d'onde SH dans le vide
    x2 = 2 * np.pi * n_med_2w * a / lam2
    m2 = m_sphere_2w / n_med_2w
    k2 = 2 * np.pi * n_med_2w / lam2

    if nmax_f is None:
        nmax_f = size_cutoff(max(abs(x1), abs(m1 * x1)))
    if nmax_sh is None:
        nmax_sh = size_cutoff(max(abs(x2), abs(m2 * x2)))
    nmax_f += int(extra_f)
    nmax_sh += int(extra_sh)
    m_values = tuple(int(mm) for mm in m_values)
    nphi = max(int(nphi), 4 * max(abs(mm) for mm in m_values) + 4)   # evite le repliement en phi

    # --- Grille de quadrature (Gauss-Legendre en cos(theta), uniforme en phi) ------------------
    nth = 2 * nmax_f + nmax_sh + 8
    mu_q, w_q = leggauss(nth)
    th_q = np.arccos(mu_q)
    ph_q = np.arange(nphi) * 2 * np.pi / nphi
    dphi = 2 * np.pi / nphi

    # --- Champ fondamental et nappe de polarisation SH ---------------------------------------
    F = fundamental_fields(th_q, ph_q, x1, m1, nmax_f)
    if side == "outside":
        Er, Et, Ep = F["out"]
    elif side == "inside":
        Er, Et, Ep = F["in"]
    elif side == "mean":
        Er, Et, Ep = [0.5 * (u + v) for u, v in zip(F["out"], F["in"])]
    else:
        raise ValueError("side doit etre 'outside', 'inside' ou 'mean'")

    p_r = chi_rrr * Er ** 2 + chi_rtt * (Et ** 2 + Ep ** 2)
    p_t = 2 * chi_trt * Er * Et
    p_p = 2 * chi_trt * Er * Ep

    W = w_q[:, None] * dphi                          # poids d'integration sur dOmega

    def proj(f, g):                                  # int f conj(g) dOmega
        return np.sum(W * f * np.conj(g))

    src_scale = max(np.max(np.abs(p_r)), np.max(np.abs(p_t)), np.max(np.abs(p_p)))

    # --- Resolution (n, m) -----------------------------------------------------------------
    sol = {}
    for n in range(1, nmax_sh + 1):
        j_i, _, Dj_i, _ = _jh(n, m2 * x2)
        _, h_o, _, Dh_o = _jh(n, x2)
        for m in m_values:
            if abs(m) > n:
                continue
            Y, (Bt, Bp), (Ct, Cp) = vsh_angular(n, m, th_q, ph_q)
            norm2 = n * (n + 1.0)
            pr_nm = proj(p_r, Y)
            alpha = (proj(p_t, Bt) + proj(p_p, Bp)) / norm2
            beta = (proj(p_t, Ct) + proj(p_p, Cp)) / norm2
            if max(abs(pr_nm), abs(alpha), abs(beta)) < 1e-10 * src_scale:   # composante nulle (regle de selection)
                continue

            # TM (coefficient de N) : inconnues (a_out, a_in)
            # H-jump terms use k0*n_med^2 (= k2*n_med) so that the TE/α channels
            # match the Born Green function in a medium of index n_med (factor n^4 on |b|^2).
            # The p_r jump (first TM row) is unchanged and already matches RGD for χ_rrr, χ_rtt.
            A_tm = np.array([[Dh_o, -Dj_i], [n_med_2w * h_o, -m_sphere_2w * j_i]], complex)
            r_tm = np.array([-pr_nm / (eps_sheet * a), k0 * (n_med_2w ** 2) * alpha], complex)
            a_o, a_i = np.linalg.solve(A_tm, r_tm)

            # TE (coefficient de M) : inconnues (b_out, b_in)
            A_te = np.array([[h_o, -j_i], [n_med_2w * Dh_o, -m_sphere_2w * Dj_i]], complex)
            r_te = np.array([0.0, -k0 * (n_med_2w ** 2) * beta], complex)
            b_o, b_i = np.linalg.solve(A_te, r_te)

            sol[(n, m)] = dict(a_out=a_o, a_in=a_i, b_out=b_o, b_in=b_i,
                               pr=pr_nm, alpha=alpha, beta=beta)

    # --- Puissance SH totale (u.a.) : int r^2 |E|^2 dOmega -----------------------------------
    power = sum(n * (n + 1.0) * (abs(c["a_out"]) ** 2 + abs(c["b_out"]) ** 2)
                for (n, m), c in sol.items()) / k2 ** 2

    # contribution de chaque multipole n (somme sur m) a la puissance SH
    power_n = {}
    for (n, m), c in sol.items():
        power_n[n] = power_n.get(n, 0.0) + n * (n + 1.0) * (abs(c["a_out"]) ** 2 + abs(c["b_out"]) ** 2) / k2 ** 2

    out = dict(sol=sol, power=power, power_n=power_n, params=dict(x1=x1, x2=x2, m1=m1, m2=m2, nmax_f=nmax_f, nmax_sh=nmax_sh,
                                                 k0=k0, k2=k2, eps_sheet=eps_sheet, radius=a),
               n_med_2w=n_med_2w, grid=(th_q, ph_q, W))

    # --- Diagramme de diffusion en champ lointain --------------------------------------------
    if theta_det is None:
        theta_det = np.deg2rad(np.linspace(0.5, 179.5, 180))
    theta_det = np.asarray(theta_det, float)
    out["theta_det"] = theta_det
    out["far_field"] = far_field(sol, theta_det, np.array([0.0, np.pi / 2]), k2)
    return out


def far_field(sol, theta, phi, k2):
    """
    Intensite SH diffusee I = r^2 |E|^2 (u.a.), composantes theta (dans le plan) et phi (perpendiculaire).
    Retourne F_theta, F_phi (amplitudes, formes (nt, np)) et les intensites correspondantes.
    """
    theta = np.asarray(theta, float).reshape(-1, 1)
    phi = np.asarray(phi, float).reshape(1, -1)
    Ft = np.zeros((theta.size, phi.size), complex)
    Fp = np.zeros_like(Ft)
    for (n, m), c in sol.items():
        _, (Bt, Bp), (Ct, Cp) = vsh_angular(n, m, theta, phi)
        ca = c["a_out"] * (-1j) ** n          # N_nm -> (-i)^n   e^{i rho}/rho  B_nm
        cb = c["b_out"] * (-1j) ** (n + 1)    # M_nm -> (-i)^{n+1} e^{i rho}/rho C_nm
        Ft += ca * Bt + cb * Ct
        Fp += ca * Bp + cb * Cp
    return dict(F_theta=Ft, F_phi=Fp,
                I_theta=np.abs(Ft) ** 2 / k2 ** 2, I_phi=np.abs(Fp) ** 2 / k2 ** 2)


def shs_polarisation(res):
    """
    Intensites SHS resolues en polarisation (indice 1 : polarisation incidente, 2 : detectee).
    Incident x-polarise.  Plan phi=0 : incident "H" (dans le plan de diffusion) ; plan phi=90 deg :
    incident "V" (perpendiculaire au plan).  Detection V = composante phi, H = composante theta.
    Retour : dict 'HH','HV','VH','VV' -> I(theta_det)
    """
    ff = res["far_field"]
    return {
        "HH": ff["I_theta"][:, 0], "HV": ff["I_phi"][:, 0],
        "VH": ff["I_theta"][:, 1], "VV": ff["I_phi"][:, 1],
    }


def shs_polarisation_signed(res, angles_deg=None):
    """
    Diagramme de diffusion SHS resolu en polarisation pour un angle de diffusion SIGNE dans [-180, 180] deg.

    Convention : angle = 0 -> diffusion vers l'avant (direction +z de propagation du faisceau),
                 |angle| = 180 deg -> retrodiffusion (backscattering).
    Un angle negatif correspond a l'autre cote du faisceau dans le meme plan de diffusion
    (azimut phi + 180 deg). Incident H : plan phi = 0 ; incident V : plan phi = 90 deg.

    Retour : dict 'angle' (deg), 'HH', 'HV', 'VH', 'VV' (1er indice : incident, 2e : detecte).
    """
    if angles_deg is None:
        angles_deg = np.linspace(-179.5, 179.5, 360)       # evite les poles theta = 0 et 180 exacts
    ang = np.asarray(angles_deg, float)
    th = np.deg2rad(np.abs(ang))
    ff = far_field(res["sol"], th, np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2]), res["params"]["k2"])
    neg = ang < 0
    iH = np.where(neg, 2, 0)          # colonne phi = 0 ou phi = pi
    iV = np.where(neg, 3, 1)          # colonne phi = pi/2 ou 3 pi/2
    r = np.arange(ang.size)
    return {"angle": ang,
            "HH": ff["I_theta"][r, iH], "HV": ff["I_phi"][r, iH],
            "VH": ff["I_theta"][r, iV], "VV": ff["I_phi"][r, iV]}


def shs_vs_gamma(res, theta_scat_deg=90.0, gamma_deg=None):
    """
    SHS resolue en polarisation en fonction de l'angle de polarisation incidente gamma, a angle de
    diffusion fixe.

    Convention : gamma = 0 -> polarisation incidente V (perpendiculaire au plan de diffusion),
                 gamma = 90 deg -> polarisation incidente H (dans le plan de diffusion).
    Detection : V = composante phi, H = composante theta.

    La sphere etant invariante par rotation autour de z, tourner la polarisation incidente d'un
    angle gamma revient a tourner le plan de diffusion : on evalue le champ lointain (calcule pour
    une polarisation incidente suivant x) a l'azimut phi_d = 90 deg - gamma.

    Retour : dict 'gamma' (deg), 'I_V', 'I_H' (u.a.), pour chaque angle de diffusion demande
             (tableaux de forme (n_gamma, n_theta)).
    """
    if gamma_deg is None:
        gamma_deg = np.linspace(0.0, 360.0, 361)
    gamma_deg = np.asarray(gamma_deg, float)
    theta_scat = np.deg2rad(np.atleast_1d(np.asarray(theta_scat_deg, float)))
    phi_d = np.pi / 2 - np.deg2rad(gamma_deg)
    ff = far_field(res["sol"], theta_scat, phi_d, res["params"]["k2"])   # formes (n_theta, n_gamma)
    return dict(gamma=gamma_deg, theta_scat=np.rad2deg(theta_scat),
                I_V=ff["I_phi"].T, I_H=ff["I_theta"].T)


def convergence_test(extra_sh_list=(-4, -2, 0, 2, 4, 8), extra_f=0, verbose=True, **kwargs):
    """
    Test de convergence en nombre de multipoles SH : calcule la puissance SH totale pour differents
    nmax_sh = (critere de Wiscombe) + extra_sh. kwargs = arguments de shs_sphere (radius, wavelength, ...).
    Retourne une liste de (nmax_sh, puissance, variation relative par rapport au point precedent).
    """
    rows, prev = [], None
    for e in extra_sh_list:
        r = shs_sphere(extra_f=extra_f, extra_sh=e, **kwargs)
        P = r["power"]
        d = np.nan if prev is None else abs(P - prev) / abs(P)
        rows.append((r["params"]["nmax_sh"], P, d))
        prev = P
        if verbose:
            print("nmax_sh = %3d   P = %.8e   variation relative = %.2e" % rows[-1])
    return rows


# ----------------------------------------------------------------------------------------------
#  Exemple
# ----------------------------------------------------------------------------------------------
if __name__ == "__main__" and len(__import__("sys").argv) <= 1:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Particule de polystyrene dans l'eau, excitation a 800 nm (valeurs d'indice typiques)
    res = shs_sphere(radius=50.0, wavelength=800.0,
                     m_sphere_w=1.50, m_sphere_2w=1.50,
                     n_med_w=1.33, n_med_2w=1.33,
                     chi_rrr=1.0, chi_rtt=0.0, chi_trt=0.0, side="outside")

    print("x(w) = %.3f   x(2w) = %.3f   n_max(w) = %d   n_max(2w) = %d"
          % (res["params"]["x1"].real, res["params"]["x2"].real,
             res["params"]["nmax_f"], res["params"]["nmax_sh"]))
    print("Puissance SH totale (u.a.) = %.6e" % res["power"])

    print("Contribution des multipoles n a la puissance SH :")
    for n_, p_ in sorted(res["power_n"].items()):
        print("   n = %2d : %.3e  (%.1f %%)" % (n_, p_, 100 * p_ / res["power"]))
    print("Test de convergence en n_max (SH) :")
    convergence_test(radius=100.0, wavelength=800.0, m_sphere_w=1.57, m_sphere_2w=1.60,
                     n_med_w=1.329, n_med_2w=1.339, chi_rrr=1.0, chi_rtt=0.5, chi_trt=0.5)

    I = shs_polarisation(res)
    th = np.rad2deg(res["theta_det"])

    # --- Figure 1 : diagrammes angulaires, angle de diffusion de -180 a +180 deg ---------------
    S = shs_polarisation_signed(res)
    ang = S["angle"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for key in ("HH", "HV", "VH", "VV"):
        ax[0].plot(ang, S[key], label=key)
    ax[0].set_xlim(-180, 180)
    ax[0].set_xticks(np.arange(-180, 181, 60))
    ax[0].axvline(0, color="grey", lw=0.5)
    ax[0].set_xlabel("angle de diffusion (deg)   [0 = avant, +-180 = retrodiffusion]")
    ax[0].set_ylabel("I_SHS (u.a.)")
    ax[0].legend()
    ax[0].set_title("SHS resolue en polarisation")
    ax[1].remove()
    axp = fig.add_subplot(1, 2, 2, projection="polar")
    axp.plot(np.deg2rad(ang), S["VV"] + S["VH"], label="inc. V")
    axp.plot(np.deg2rad(ang), S["HH"] + S["HV"], label="inc. H")
    axp.set_title("Diagramme polaire (0 = avant, 180 = arriere)")
    axp.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig("shs_mie_exemple.png", dpi=130)

    # --- Figure 2 : SHS en fonction de l'angle de polarisation incidente gamma ---------------
    angles = [90.0, 60.0, 120.0]                       # angles de diffusion (deg)
    g = shs_vs_gamma(res, theta_scat_deg=angles)
    fig2, axs = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
    for i, ang in enumerate(angles):
        axs[0].plot(g["gamma"], g["I_V"][:, i], label="theta = %g deg" % ang)
        axs[1].plot(g["gamma"], g["I_H"][:, i], label="theta = %g deg" % ang)
    axs[0].set_title("Detection verticale (V)")
    axs[1].set_title("Detection horizontale (H)")
    for a_ in axs:
        a_.set_xlabel("angle de polarisation incidente gamma (deg)  [0 = V, 90 = H]")
        a_.set_ylabel("I_SHS (u.a.)")
        a_.set_xticks(np.arange(0, 361, 45))
        a_.legend()
    fig2.tight_layout()
    fig2.savefig("shs_mie_gamma.png", dpi=130)

    # Export des courbes (theta = 90 deg) : gamma, I_V, I_H
    np.savetxt("shs_mie_gamma_90deg.csv",
               np.column_stack([g["gamma"], g["I_V"][:, 0], g["I_H"][:, 0]]),
               delimiter=",", header="gamma_deg,I_V,I_H", comments="")


# ----------------------------------------------------------------------------------------------
#  CLI wrapper for PySHS GUI (writes out_plot compatible with graph_polar / graph_angle)
# ----------------------------------------------------------------------------------------------
def _parse_beta_file(path):
    """Parse PySHS beta tensor file → dict label→value (labels like zzz, zxx, xzx)."""
    beta = {}
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    beta[parts[0].lower()] = float(parts[1])
                except ValueError:
                    continue
    return beta


def _chi_from_beta(beta):
    """
    Map surface hyperpolarizability (local frame: z = surface normal) to Mie chi.

    Isotropic achiral surface (standard convention, consistent with NLRGD):
      χ_rrr = β_zzz
      χ_rtt = β_zxx = β_zyy          (each equals χ_⊥∥∥)
      χ_trt = β_xzx = β_xxz = β_yzy = β_yyz   (each equals χ_∥⊥∥)

    Extraction uses the mean of the symmetry-related components so that
    p_r = χ_rtt (E_t²+E_p²) and p_t = 2 χ_trt E_r E_t match the tensor contraction
    P_i = Σ_jk β_ijk E_j E_k used in the RGD/NLRGD code.
    """
    chi_rrr = beta.get("zzz", 0.0)
    zxx, zyy = beta.get("zxx", 0.0), beta.get("zyy", 0.0)
    chi_rtt = 0.5 * (zxx + zyy) if (zxx or zyy) else 0.0
    comps = [beta.get(k, 0.0) for k in ("xzx", "xxz", "yzy", "yyz")]
    chi_trt = (sum(comps) / 4.0) if any(comps) else 0.0
    return chi_rrr, chi_rtt, chi_trt




_USED_COMPONENTS = ("zzz", "zxx", "zyy", "xzx", "xxz", "yzy", "yyz")


def _check_beta_usage(beta):
    """Messages telling which parts of the tensor the Mie model cannot use.

    The nonlinear Mie model has only three independent surface constants (isotropic achiral
    surface): chi_rrr = zzz, chi_rtt = (zxx + zyy)/2, chi_trt = mean(xzx, xxz, yzy, yyz).
    Everything else (chiral components such as xyz, in-plane anisotropy, ...) is ignored, and
    zxx != zyy or unequal chi_trt components are averaged.
    """
    msgs = []
    scale = max([abs(v) for v in beta.values()] + [0.0])
    if scale == 0.0:
        return msgs
    tol = 1e-9 * scale
    ignored = sorted(k for k, v in beta.items() if k not in _USED_COMPONENTS and abs(v) > tol)
    if ignored:
        msgs.append("WARNING: the Mie model only uses zzz, zxx, zyy, xzx, xxz, yzy, yyz (isotropic achiral "
                    "surface); these non-zero components are IGNORED: %s" % " ".join(ignored))
    zxx, zyy = beta.get("zxx", 0.0), beta.get("zyy", 0.0)
    if abs(zxx - zyy) > tol:
        msgs.append("WARNING: zxx != zyy (%.6g, %.6g): their mean is used for chi_rtt." % (zxx, zyy))
    comps = [beta.get(k, 0.0) for k in ("xzx", "xxz", "yzy", "yyz")]
    if max(comps) - min(comps) > tol:
        msgs.append("WARNING: xzx, xxz, yzy, yyz differ: their mean is used for chi_trt.")
    return msgs


def run_mie_cli(keyword, beta_path, out_path, answers, work_dir="."):
    import os
    """
    Run Mie SHS calculation and write out_plot + output.txt for the GUI.

    answers (list of strings, same order as GUI stdin for Sphere):
      wavelength, n_p_real, n_p_imag, n_m_real, n_m_imag, radius,
      then depending on keyword:
        polarplot_single: theta
        polarplot_integrate: theta_start, theta_end, n_points
        angle_scattering: n_points
    """
    import os
    wavelength = float(answers[0])
    np_real, np_imag = float(answers[1]), float(answers[2])
    nm_real, nm_imag = float(answers[3]), float(answers[4])
    radius = float(answers[5])

    m_sphere_w = complex(np_real, np_imag)
    m_sphere_2w = complex(np_real, np_imag)  # same index at 2ω if not specified separately
    n_med_w = complex(nm_real, nm_imag)
    n_med_2w = complex(nm_real, nm_imag)

    beta = _parse_beta_file(beta_path)
    for _m in _check_beta_usage(beta):
        print(_m)
    chi_rrr, chi_rtt, chi_trt = _chi_from_beta(beta)

    print("PySHS – Sphere SHS (nonlinear Mie theory)")
    print("wavelength = %.2f nm, radius = %.2f nm" % (wavelength, radius))
    print("n_p = %.5f + i*%.5f" % (np_real, np_imag))
    print("n_m = %.5f + i*%.5f" % (nm_real, nm_imag))
    print("chi_rrr = %.6g  chi_rtt = %.6g  chi_trt = %.6g" % (chi_rrr, chi_rtt, chi_trt))
    print("keyword = %s" % keyword)
    # Units: the nonlinear Mie intensity is already r^2|E_2w|^2 (radiative prefactor
    # included, "rad").  PYSHS_UNITS=raw removes the prefactor (see shs_units.py).
    unit = shs_units.settings(4.0 * math.pi * complex(nm_real, nm_imag).real / wavelength, "rad")
    i_scale = unit.i_scale if unit.mode == "rad" else unit.i_scale / (unit.pref * unit.pref)
    for _l in unit.lines():
        print(_l)
    print("Computing fundamental + SH multipoles…")

    res = shs_sphere(
        radius=radius, wavelength=wavelength,
        m_sphere_w=m_sphere_w, m_sphere_2w=m_sphere_2w,
        n_med_w=n_med_w.real, n_med_2w=n_med_2w.real,  # medium real for Mie (as in theory)
        chi_rrr=chi_rrr, chi_rtt=chi_rtt, chi_trt=chi_trt,
        side="outside",
    )
    print("x(w) = %.3f  x(2w) = %.3f  nmax_f = %d  nmax_sh = %d"
          % (res["params"]["x1"].real, res["params"]["x2"].real,
             res["params"]["nmax_f"], res["params"]["nmax_sh"]))
    print("Total SH power (a.u.) = %.6e" % res["power"])

    out_plot_path = os.path.join(work_dir, "out_plot")

    with open(out_path, "w", encoding="utf-8") as fout:
        fout.write("********************************************************\n")
        fout.write("****     Py_SHS: Sphere SHS – nonlinear Mie theory  ****\n")
        fout.write("********************************************************\n")
        fout.write("Mode                : Nonlinear Mie\n")
        fout.write("Wavelength          : %.2f nm\n" % wavelength)
        fout.write("n_particle (n_p)    : %.5f + i*%.5f\n" % (np_real, np_imag))
        fout.write("n_medium   (n_m)    : %.5f + i*%.5f\n" % (nm_real, nm_imag))
        fout.write("Radius              : %.2f nm\n" % radius)
        fout.write("chi_rrr=%.6g  chi_rtt=%.6g  chi_trt=%.6g\n" % (chi_rrr, chi_rtt, chi_trt))
        fout.write("x(w)=%.4f  x(2w)=%.4f  power=%.6e\n" % (
            res["params"]["x1"].real, res["params"]["x2"].real, res["power"]))
        fout.write("********************************************************\n\n")

        if keyword in ("polarplot_single", "polarplot_integrate"):
            if keyword == "polarplot_single":
                theta_list = [float(answers[6])]
            else:
                t0 = float(answers[6])
                t1 = float(answers[7])
                npt = max(int(answers[8]), 2)
                theta_list = list(np.linspace(t0, t1, npt))

            # Average I_V(gamma), I_H(gamma) over requested angles
            g_all = None
            IV_sum = None
            IH_sum = None
            for th in theta_list:
                print("  Mie polar plot at Theta = %.1f deg" % th)
                g = shs_vs_gamma(res, theta_scat_deg=th, gamma_deg=np.linspace(0.0, 360.0, 101))
                if g_all is None:
                    g_all = g["gamma"]
                    IV_sum = np.array(g["I_V"][:, 0], dtype=float)
                    IH_sum = np.array(g["I_H"][:, 0], dtype=float)
                else:
                    IV_sum += np.array(g["I_V"][:, 0], dtype=float)
                    IH_sum += np.array(g["I_H"][:, 0], dtype=float)
            n = len(theta_list)
            IV_sum /= n
            IH_sum /= n

            fout.write("Mie polar plot (averaged over %d angle(s))\n" % n)
            with open(out_plot_path, "w", encoding="utf-8") as fp:
                fp.write("gamma IV(gamma) IH(gamma)\n")
                for i in range(len(g_all)):
                    # out_plot uses gamma in radians for polarplot (GUI converts if <= 2pi)
                    gam_rad = g_all[i] * np.pi / 180.0
                    fp.write("%.8e %.8e %.8e\n" % (gam_rad, IV_sum[i]*i_scale, IH_sum[i]*i_scale))
                    fout.write("%.8e %.8e %.8e\n" % (gam_rad, IV_sum[i]*i_scale, IH_sum[i]*i_scale))

        elif keyword == "angle_scattering":
            npt = max(int(answers[6]), 2)
            # signed angles from -180 to +180, matching C code style (-pi … +pi)
            angles = np.degrees(shs_angles.scattering_angles(npt))
            print("  Mie angle_scattering: %d points" % len(angles))
            S = shs_polarisation_signed(res, angles_deg=angles)
            # PySHS convention: I_(detected)(gamma_incident)
            #   IH(0°)  = IVH  (incident V, detect H)
            #   IH(90°) = IHH  (incident H, detect H)
            #   IV(0°)  = IVV  (incident V, detect V)
            #   IV(90°) = IHV  (incident H, detect V)
            # shs_mie keys are Pin-Pout: HH, HV, VH, VV
            fout.write(" grandtheta  IPSS-IH(0)  IPPP-IH(90)  ISSS-IV(0)  ISPP-IV(90)\n")
            with open(out_plot_path, "w", encoding="utf-8") as fp:
                fp.write(" grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n")
                for i, ang in enumerate(S["angle"]):
                    ih0 = float(S["VH"][i])   # IVH → IH(0°)
                    ih90 = float(S["HH"][i])  # IHH → IH(90°)
                    iv0 = float(S["VV"][i])   # IVV → IV(0°)
                    iv90 = float(S["HV"][i])  # IHV → IV(90°)
                    # Only positive half in C-style listing is optional; we keep full signed range
                    fp.write(" %.6f    %.8e   %.8e   %.8e   %.8e\n" % (ang, ih0*i_scale, ih90*i_scale, iv0*i_scale, iv90*i_scale))
                    fout.write(" %.6f    %.8e   %.8e   %.8e   %.8e\n" % (ang, ih0*i_scale, ih90*i_scale, iv0*i_scale, iv90*i_scale))
        else:
            raise ValueError("Unknown keyword: %s" % keyword)

    print("calculation finished, out_plot written")
    return res


if __name__ == "__main__" and len(__import__("sys").argv) > 1:
    # Usage: python3 shs_mie.py <keyword> <beta_file> <output_file>
    # stdin: wavelength n_p_r n_p_i n_m_r n_m_i radius [angles…]
    import sys
    import os
    kw = sys.argv[1]
    beta_f = sys.argv[2]
    out_f = sys.argv[3]
    ans = [line.strip() for line in sys.stdin if line.strip()]
    work = os.path.dirname(os.path.abspath(out_f)) or "."
    run_mie_cli(kw, beta_f, out_f, ans, work_dir=work)
