/*********************************************************************************
******************PROGRAM sphere_SHS_angle_scattering*****************************
**********************************************************************************
   Py_SHS - An open source software about Second Harmonic Scattering 
   developed at Institut Charles Gerhardt Montpellier - ENSCM
   Dr Pierre-Marie GASSIN
   Dr Gaelle GASSIN
    (june 2021)  

   Modified version with optional NLWKB approximation
   (Nonlinear Wentzel-Kramers-Brillouin)
*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

// CONVENTION ZYZ pour matrice Euler!!!

double Txx(double phi, double theta, double psi) {
    return -sin(phi)*sin(psi) + cos(theta)*cos(phi)*cos(psi);
}
double Txy(double phi, double theta, double psi) {
    return -sin(phi)*cos(psi) - cos(theta)*cos(phi)*sin(psi);
}
double Txz(double phi, double theta, double psi) {
    return sin(theta)*cos(phi);
}
double Tyx(double phi, double theta, double psi) {
    return cos(phi)*sin(psi) + cos(theta)*sin(phi)*cos(psi);
}
double Tyy(double phi, double theta, double psi) {
    return cos(phi)*cos(psi) - cos(theta)*sin(psi)*sin(phi);
}
double Tyz(double phi, double theta, double psi) {
    return sin(theta)*sin(phi);
}
double Tzx(double phi, double theta, double psi) {
    return -sin(theta)*cos(psi);
}
double Tzy(double phi, double theta, double psi) {
    return sin(theta)*sin(psi);
}
double Tzz(double phi, double theta, double psi) {
    return cos(theta);
}

double T(double phi, double theta, double psi, int i, int j) {
    if (i==0 && j==0) return Txx(phi,theta,psi);
    if (i==0 && j==1) return Txy(phi,theta,psi);
    if (i==0 && j==2) return Txz(phi,theta,psi);
    if (i==1 && j==0) return Tyx(phi,theta,psi);
    if (i==1 && j==1) return Tyy(phi,theta,psi);
    if (i==1 && j==2) return Tyz(phi,theta,psi);
    if (i==2 && j==0) return Tzx(phi,theta,psi);
    if (i==2 && j==1) return Tzy(phi,theta,psi);
    if (i==2 && j==2) return Tzz(phi,theta,psi);
    return 0.0;
}

double kimesoijk(double phi, double theta, double psi, int i, int j, int k, double beta[3][3][3]) {
    double uuu = 0.0;
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[0][0][0];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[0][0][1];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[0][0][2];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[0][1][0];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[0][1][1];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[0][1][2];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[0][2][0];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[0][2][1];
    uuu += T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[0][2][2];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[1][0][0];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[1][0][1];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[1][0][2];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[1][1][0];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[1][1][1];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[1][1][2];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[1][2][0];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[1][2][1];
    uuu += T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[1][2][2];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[2][0][0];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[2][0][1];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[2][0][2];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[2][1][0];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[2][1][1];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[2][1][2];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[2][2][0];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[2][2][1];
    uuu += T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[2][2][2];
    return uuu;
}

// ========== Fonctions grandgamma ==========

double complex grandgammaxxx(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,0,0,0,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammazxx(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,2,0,0,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammayxx(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,1,0,0,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammazyy(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,2,1,1,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammayyy(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,1,1,1,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammaxyy(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,0,1,1,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammaxxy(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,0,0,1,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammazxy(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,2,0,1,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

double complex grandgammayxy(double R, double phic, double thetac, double grandtheta,
                             double beta[3][3][3], double complex petitk, double complex grandk) {
    return R*R * kimesoijk(phic,thetac,0,1,0,1,beta) * sin(thetac)
           * cexp(((-2*petitk + grandk*cos(grandtheta))*R*cos(thetac)
                  + (grandk*sin(grandtheta)*R*sin(thetac)*sin(phic))) * I);
}

// ======================== MAIN ========================

int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    FILE *fichier_beta = NULL;
    FILE *fichier_out = NULL;
    FILE *fichier_grace1 = NULL;
    char chaine[121];

    printf("********************************************************\n");
    printf("****           Py_SHS: Sphere SHS module            ****\n");
    printf("****         with optional NLWKB approximation      ****\n");
    printf("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n");
    printf("****                  ICGM - ENSCM                  ****\n");
    printf("********************************************************\n\n");

    double pi = 3.141592653589793;
    double beta[3][3][3] = {0};
    double longueuronde, R;
    int use_nlwkb = 0;
    int pas = 10;
    int pas2 = 10;

    // ---------- Choix du mode ----------
    printf("Do you want to use NLWKB approximation? (1 = yes, 0 = no [classic RGD]): ");
    scanf("%d", &use_nlwkb);
    printf("\n");

    double complex n, n_p, n_m;
    double realpart, imagpart;

    printf("Enter the wavelength in nm: ");
    scanf("%lf", &longueuronde);

    if (use_nlwkb) {
        printf("\n=== NLWKB mode activated ===\n");
        printf("Enter REAL part of particle refractive index n_p: ");
        scanf("%lf", &realpart);
        printf("Enter IMAGINARY part of particle refractive index n_p: ");
        scanf("%lf", &imagpart);
        n_p = realpart + imagpart * I;

        printf("Enter REAL part of medium refractive index n_m: ");
        scanf("%lf", &realpart);
        printf("Enter IMAGINARY part of medium refractive index n_m: ");
        scanf("%lf", &imagpart);
        n_m = realpart + imagpart * I;

        printf("\nn_p (particle) = %.5f + i*%.5f\n", creal(n_p), cimag(n_p));
        printf("n_m (medium)   = %.5f + i*%.5f\n", creal(n_m), cimag(n_m));
    } else {
        printf("\n=== Classic RGD mode (single index) ===\n");
        printf("Enter REAL part of the refractive index n: ");
        scanf("%lf", &realpart);
        printf("Enter IMAGINARY part of the refractive index n: ");
        scanf("%lf", &imagpart);
        n = realpart + imagpart * I;
        n_p = n;
        n_m = n;
        printf("\nn = %.5f + i*%.5f\n", creal(n), cimag(n));
    }

    printf("Enter the radius of the sphere in nm: ");
    scanf("%lf", &R);
    printf("\n");

    // ---------- Calcul des vecteurs d'onde ----------
    double complex petitk, grandk;

    if (use_nlwkb) {
        petitk = 2.0 * pi * n_p / longueuronde;   // k_ω inside particle
        grandk = 4.0 * pi * n_m / longueuronde;   // k_2ω in medium
        printf("NLWKB vectors:\n");
    } else {
        petitk = 2.0 * pi * n / longueuronde;
        grandk = 4.0 * pi * n / longueuronde;
        printf("Classic RGD vectors:\n");
    }

    printf("petitk (k_omega) = %.6f + i*%.6f\n", creal(petitk), cimag(petitk));
    printf("grandk (k_2omega)= %.6f + i*%.6f\n\n", creal(grandk), cimag(grandk));

    // ---------- Vérification des arguments ----------
    if (argc < 4) {
        printf("Usage: %s <keyword> <ki2_inputfile> <outputfile>\n", argv[0]);
        printf("keyword = polarplot_single | polarplot_integrate | angle_scattering\n");
        return 1;
    }

    // ---------- Lecture du tenseur β ----------
    fichier_beta = fopen(argv[2], "rt");
    if (fichier_beta == NULL) {
        printf("Problem to open the input ki(2) file!\n");
        return 1;
    }

    fgets(chaine, 121, fichier_beta); // header
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                fgets(chaine, 121, fichier_beta);
                sscanf(chaine, "%*s%lf", &beta[i][j][k]);
            }
        }
    }
    fclose(fichier_beta);

    fichier_out = fopen(argv[3], "w");
    if (fichier_out == NULL) {
        printf("Problem to open the output file!\n");
        return 1;
    }

    // ---------- En-tête du fichier de sortie ----------
    fprintf(fichier_out, "********************************************************\n");
    fprintf(fichier_out, "****           Py_SHS: Sphere SHS module            ****\n");
    fprintf(fichier_out, "****         with optional NLWKB approximation      ****\n");
    fprintf(fichier_out, "********************************************************\n");
    fprintf(fichier_out, "Mode                : %s\n", use_nlwkb ? "NLWKB" : "Classic RGD");
    fprintf(fichier_out, "Wavelength          : %.2f nm\n", longueuronde);
    if (use_nlwkb) {
        fprintf(fichier_out, "n_particle (n_p)    : %.5f + i*%.5f\n", creal(n_p), cimag(n_p));
        fprintf(fichier_out, "n_medium   (n_m)    : %.5f + i*%.5f\n", creal(n_m), cimag(n_m));
    } else {
        fprintf(fichier_out, "Refractive index n  : %.5f + i*%.5f\n", creal(n), cimag(n));
    }
    fprintf(fichier_out, "Radius              : %.2f nm\n", R);
    fprintf(fichier_out, "petitk              : %.6f + i*%.6f\n", creal(petitk), cimag(petitk));
    fprintf(fichier_out, "grandk              : %.6f + i*%.6f\n", creal(grandk), cimag(grandk));
    fprintf(fichier_out, "********************************************************\n\n");

    // Affichage du tenseur
    fprintf(fichier_out, "*********          SURFACE KI(2) TENSOR      ***********\n\n");
    fprintf(fichier_out, "    xxx=%lf   xxy=%lf   xxz= %lf\n", beta[0][0][0], beta[0][0][1], beta[0][0][2]);
    fprintf(fichier_out, "    xyx=%lf   xyy=%lf   xyz= %lf\n", beta[0][1][0], beta[0][1][1], beta[0][1][2]);
    fprintf(fichier_out, "    xzx=%lf   xzy=%lf   xzz= %lf\n", beta[0][2][0], beta[0][2][1], beta[0][2][2]);
    fprintf(fichier_out, "\n");
    fprintf(fichier_out, "    yxx=%lf   yxy=%lf   yxz= %lf\n", beta[1][0][0], beta[1][0][1], beta[1][0][2]);
    fprintf(fichier_out, "    yyx=%lf   yyy=%lf   yyz= %lf\n", beta[1][1][0], beta[1][1][1], beta[1][1][2]);
    fprintf(fichier_out, "    yzx=%lf   yzy=%lf   yzz= %lf\n", beta[1][2][0], beta[1][2][1], beta[1][2][2]);
    fprintf(fichier_out, "\n");
    fprintf(fichier_out, "    zxx=%lf   zxy=%lf   zxz= %lf\n", beta[2][0][0], beta[2][0][1], beta[2][0][2]);
    fprintf(fichier_out, "    zyx=%lf   zyy=%lf   zyz= %lf\n", beta[2][1][0], beta[2][1][1], beta[2][1][2]);
    fprintf(fichier_out, "    zzx=%lf   zzy=%lf   zzz= %lf\n", beta[2][2][0], beta[2][2][1], beta[2][2][2]);
    fprintf(fichier_out, "\n*******************************************************\n\n");

    char str1[30] = "polarplot_single";
    char str2[30] = "polarplot_integrate";
    char str3[30] = "angle_scattering";

    double hx, hy, xi, yj, compteur;
    int nx, ny, i, j, iii, ic, mmm, mmmf, point;
    double complex grandgammazxx_moy, grandgammayyy_moy, grandgammazyy_moy;
    double complex grandgammaxxx_moy, grandgammayxx_moy, grandgammayxy_moy;
    double complex grandgammaxyy_moy, grandgammaxxy_moy, grandgammazxy_moy;
    double av, bv, cv, ah, bh, ch;
    double I4v, I2v, I2h, I4h;
    double Iv, Ih, gammaaa, nomc, noms;
    double grandtheta, grandthetadegree, grandtheta_1, grandtheta_2;
    double grandthetadegree_1, grandthetadegree_2;
    double IH_0, IH_90, IV_0, IV_90;

    // =====================================================================
    // MODE 1 : polarplot_single
    // =====================================================================
    if (strcmp(str1, argv[1]) == 0) {

        printf("Enter the scattered angle in degree (0 = transmission): ");
        scanf("%lf", &grandthetadegree);
        grandtheta = grandthetadegree * 2.0 * pi / 360.0;
        printf("\n********** Begin calculation for sphere polar plot at %.1f ° ********\n", grandthetadegree);

        av = bv = cv = ah = bh = ch = 0.0;

        hx = (2.0 * pi) / (2.0 * pas);
        hy = pi / pas;
        nx = 2 * pas;
        ny = pas;

        grandgammazxx_moy = grandgammayyy_moy = grandgammazyy_moy = 0.0 + 0.0*I;
        grandgammaxxx_moy = grandgammayxx_moy = grandgammayxy_moy = 0.0 + 0.0*I;
        grandgammaxyy_moy = grandgammaxxy_moy = grandgammazxy_moy = 0.0 + 0.0*I;

        for (i = 0; i < nx; i++) {
            compteur = i * 100.0 / nx;
            printf("calculation in progress : %.1f %%\n", compteur);
            for (j = 0; j < ny; j++) {
                xi = hx/2.0 + i * hx;
                yj = hy/2.0 + j * hy;

                grandgammazxx_moy += hx*hy * grandgammazxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammaxxx_moy += hx*hy * grandgammaxxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammayxx_moy += hx*hy * grandgammayxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammaxyy_moy += hx*hy * grandgammaxyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammazyy_moy += hx*hy * grandgammazyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammayyy_moy += hx*hy * grandgammayyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammaxxy_moy += hx*hy * grandgammaxxy(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammayxy_moy += hx*hy * grandgammayxy(R, xi, yj, grandtheta, beta, petitk, grandk);
                grandgammazxy_moy += hx*hy * grandgammazxy(R, xi, yj, grandtheta, beta, petitk, grandk);
            }
        }

        av = creal(grandgammaxxx_moy * conj(grandgammaxxx_moy));
        cv = creal(grandgammaxyy_moy * conj(grandgammaxyy_moy));
        bv = creal(4.0*grandgammaxxy_moy*conj(grandgammaxxy_moy)
                 + grandgammaxxx_moy*conj(grandgammaxyy_moy)
                 + grandgammaxyy_moy*conj(grandgammaxxx_moy));

        ah = creal(grandgammayxx_moy*conj(grandgammayxx_moy)*cos(grandtheta)*cos(grandtheta)
                 + sin(grandtheta)*sin(grandtheta)*grandgammazxx_moy*conj(grandgammazxx_moy)
                 - cos(grandtheta)*sin(grandtheta)*(grandgammayxx_moy*conj(grandgammazxx_moy)
                                                   + grandgammazxx_moy*conj(grandgammayxx_moy)));

        ch = creal(grandgammayyy_moy*conj(grandgammayyy_moy)*cos(grandtheta)*cos(grandtheta)
                 + sin(grandtheta)*sin(grandtheta)*grandgammazyy_moy*conj(grandgammazyy_moy)
                 - cos(grandtheta)*sin(grandtheta)*(grandgammayyy_moy*conj(grandgammazyy_moy)
                                                   + grandgammazyy_moy*conj(grandgammayyy_moy)));

        bh = creal(cos(grandtheta)*cos(grandtheta)*(4.0*grandgammayxy_moy*conj(grandgammayxy_moy)
                                                  + grandgammayxx_moy*conj(grandgammayyy_moy)
                                                  + grandgammayyy_moy*conj(grandgammayxx_moy))
                 + sin(grandtheta)*sin(grandtheta)*(4.0*grandgammazxy_moy*conj(grandgammazxy_moy)
                                                  + grandgammazxx_moy*conj(grandgammazyy_moy)
                                                  + grandgammazyy_moy*conj(grandgammazxx_moy))
                 - cos(grandtheta)*sin(grandtheta)*(4.0*grandgammayxy_moy*conj(grandgammazxy_moy)
                                                  + 4.0*grandgammazxy_moy*conj(grandgammayxy_moy)
                                                  + grandgammayxx_moy*conj(grandgammazyy_moy)
                                                  + grandgammazyy_moy*conj(grandgammayxx_moy)
                                                  + grandgammazxx_moy*conj(grandgammayyy_moy)
                                                  + grandgammayyy_moy*conj(grandgammazxx_moy)));

        printf("\n********** calculation complete ************\n");
        printf("av = %lf\nbv = %lf\ncv = %lf\nah = %lf\nbh = %lf\nch = %lf\n", av, bv, cv, ah, bh, ch);

        I4v = (av - bv + cv) / (3.0*av + bv + 3.0*cv);
        I2v = 4.0 * (av - cv) / (3.0*av + bv + 3.0*cv);
        I4h = (ah - bh + ch) / (3.0*ah + bh + 3.0*ch);
        I2h = 4.0 * (ah - ch) / (3.0*ah + bh + 3.0*ch);

        printf("I2v = %lf\nI4v = %lf\nI2h = %lf\nI4h = %lf\n", I2v, I4v, I2h, I4h);

        fichier_grace1 = fopen("out_plot", "w");
        fprintf(fichier_out, "*********************************************\n");
        fprintf(fichier_out, "* Sphere POLAR PLOT - CONFIGURATION AT %.1f °\n", grandthetadegree);
        fprintf(fichier_out, "*********************************************\n");
        fprintf(fichier_out, "I2v=%lf  I4v=%lf\nI2h=%lf  I4h=%lf\n", I2v, I4v, I2h, I4h);
        fprintf(fichier_out, "av=%lf   bv=%lf   cv=%lf\nah=%lf   bh=%lf   ch=%lf\n", av, bv, cv, ah, bh, ch);

        fprintf(fichier_grace1, "gamma IV(gamma) IH(gamma)\n");
        for (iii = 0; iii < 100; iii++) {
            gammaaa = 2.0 * pi * iii / 100.0;
            nomc = cos(gammaaa);
            noms = sin(gammaaa);
            Iv = av*nomc*nomc*nomc*nomc + bv*nomc*nomc*noms*noms + cv*noms*noms*noms*noms;
            Ih = ah*nomc*nomc*nomc*nomc + bh*nomc*nomc*noms*noms + ch*noms*noms*noms*noms;
            fprintf(fichier_grace1, "%lf %lf %lf\n", gammaaa, Iv, Ih);
        }
        fclose(fichier_grace1);
        printf("-----------------------------------------------------------------------------\n");
        printf("calculation finished, open the out_plot file with gnuplot\n");
        printf("-----------------------------------------------------------------------------\n");
    }

    // =====================================================================
    // MODE 2 : polarplot_integrate
    // =====================================================================
    else if (strcmp(str2, argv[1]) == 0) {

        printf("Enter the beginning scattered angle in degree (0 = transmission): ");
        scanf("%lf", &grandthetadegree_1);
        printf("Enter the final scattered angle in degree (0 = transmission): ");
        scanf("%lf", &grandthetadegree_2);
        printf("Number of integrale points: ");
        scanf("%d", &point);

        grandtheta_1 = grandthetadegree_1 * 2.0 * pi / 360.0;
        grandtheta_2 = grandthetadegree_2 * 2.0 * pi / 360.0;

        av = bv = cv = ah = bh = ch = 0.0;

        hx = (2.0 * pi) / (2.0 * pas);
        hy = pi / pas;
        nx = 2 * pas;
        ny = pas;

        for (ic = 0; ic < point; ic++) {
            grandtheta = grandtheta_1 + ic * (grandtheta_2 - grandtheta_1) / (point - 1.0);
            grandthetadegree = grandtheta * 360.0 / (2.0 * pi);
            printf("********** Begin calculation for sphere polar plot at %.1f ° ********\n", grandthetadegree);

            grandgammazxx_moy = grandgammayyy_moy = grandgammazyy_moy = 0.0 + 0.0*I;
            grandgammaxxx_moy = grandgammayxx_moy = grandgammayxy_moy = 0.0 + 0.0*I;
            grandgammaxyy_moy = grandgammaxxy_moy = grandgammazxy_moy = 0.0 + 0.0*I;

            for (i = 0; i < nx; i++) {
                for (j = 0; j < ny; j++) {
                    xi = hx/2.0 + i * hx;
                    yj = hy/2.0 + j * hy;

                    grandgammazxx_moy += hx*hy * grandgammazxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammaxxx_moy += hx*hy * grandgammaxxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammayxx_moy += hx*hy * grandgammayxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammaxyy_moy += hx*hy * grandgammaxyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammazyy_moy += hx*hy * grandgammazyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammayyy_moy += hx*hy * grandgammayyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammaxxy_moy += hx*hy * grandgammaxxy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammayxy_moy += hx*hy * grandgammayxy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammazxy_moy += hx*hy * grandgammazxy(R, xi, yj, grandtheta, beta, petitk, grandk);
                }
            }

            av += creal(grandgammaxxx_moy * conj(grandgammaxxx_moy));
            cv += creal(grandgammaxyy_moy * conj(grandgammaxyy_moy));
            bv += creal(4.0*grandgammaxxy_moy*conj(grandgammaxxy_moy)
                      + grandgammaxxx_moy*conj(grandgammaxyy_moy)
                      + grandgammaxyy_moy*conj(grandgammaxxx_moy));

            ah += creal(grandgammayxx_moy*conj(grandgammayxx_moy)*cos(grandtheta)*cos(grandtheta)
                      + sin(grandtheta)*sin(grandtheta)*grandgammazxx_moy*conj(grandgammazxx_moy)
                      - cos(grandtheta)*sin(grandtheta)*(grandgammayxx_moy*conj(grandgammazxx_moy)
                                                        + grandgammazxx_moy*conj(grandgammayxx_moy)));

            ch += creal(grandgammayyy_moy*conj(grandgammayyy_moy)*cos(grandtheta)*cos(grandtheta)
                      + sin(grandtheta)*sin(grandtheta)*grandgammazyy_moy*conj(grandgammazyy_moy)
                      - cos(grandtheta)*sin(grandtheta)*(grandgammayyy_moy*conj(grandgammazyy_moy)
                                                        + grandgammazyy_moy*conj(grandgammayyy_moy)));

            bh += creal(cos(grandtheta)*cos(grandtheta)*(4.0*grandgammayxy_moy*conj(grandgammayxy_moy)
                                                       + grandgammayxx_moy*conj(grandgammayyy_moy)
                                                       + grandgammayyy_moy*conj(grandgammayxx_moy))
                      + sin(grandtheta)*sin(grandtheta)*(4.0*grandgammazxy_moy*conj(grandgammazxy_moy)
                                                       + grandgammazxx_moy*conj(grandgammazyy_moy)
                                                       + grandgammazyy_moy*conj(grandgammazxx_moy))
                      - cos(grandtheta)*sin(grandtheta)*(4.0*grandgammayxy_moy*conj(grandgammazxy_moy)
                                                       + 4.0*grandgammazxy_moy*conj(grandgammayxy_moy)
                                                       + grandgammayxx_moy*conj(grandgammazyy_moy)
                                                       + grandgammazyy_moy*conj(grandgammayxx_moy)
                                                       + grandgammazxx_moy*conj(grandgammayyy_moy)
                                                       + grandgammayyy_moy*conj(grandgammazxx_moy)));
        }

        av /= point; cv /= point; bv /= point;
        ah /= point; ch /= point; bh /= point;

        printf("\n********** calculation complete ************\n");
        printf("av = %lf\nbv = %lf\ncv = %lf\nah = %lf\nbh = %lf\nch = %lf\n", av, bv, cv, ah, bh, ch);

        I4v = (av - bv + cv) / (3.0*av + bv + 3.0*cv);
        I2v = 4.0 * (av - cv) / (3.0*av + bv + 3.0*cv);
        I4h = (ah - bh + ch) / (3.0*ah + bh + 3.0*ch);
        I2h = 4.0 * (ah - ch) / (3.0*ah + bh + 3.0*ch);

        printf("I2v = %lf\nI4v = %lf\nI2h = %lf\nI4h = %lf\n", I2v, I4v, I2h, I4h);

        fichier_grace1 = fopen("out_plot", "w");
        fprintf(fichier_out, "***********************************************************************\n");
        fprintf(fichier_out, "* Sphere POLAR PLOT - Average between %.1f ° and %.1f °\n", grandthetadegree_1, grandthetadegree_2);
        fprintf(fichier_out, "***********************************************************************\n");
        fprintf(fichier_out, "I2v=%lf  I4v=%lf\nI2h=%lf  I4h=%lf\n", I2v, I4v, I2h, I4h);
        fprintf(fichier_out, "av=%lf   bv=%lf   cv=%lf\nah=%lf   bh=%lf   ch=%lf\n", av, bv, cv, ah, bh, ch);

        fprintf(fichier_grace1, "gamma IV(gamma) IH(gamma)\n");
        for (iii = 0; iii < 100; iii++) {
            gammaaa = 2.0 * pi * iii / 100.0;
            nomc = cos(gammaaa);
            noms = sin(gammaaa);
            Iv = av*nomc*nomc*nomc*nomc + bv*nomc*nomc*noms*noms + cv*noms*noms*noms*noms;
            Ih = ah*nomc*nomc*nomc*nomc + bh*nomc*nomc*noms*noms + ch*noms*noms*noms*noms;
            fprintf(fichier_grace1, "%lf %lf %lf\n", gammaaa, Iv, Ih);
        }
        fclose(fichier_grace1);
        printf("-----------------------------------------------------------------------------\n");
        printf("calculation finished, open the out_plot file with gnuplot\n");
        printf("-----------------------------------------------------------------------------\n");
    }

    // =====================================================================
    // MODE 3 : angle_scattering
    // =====================================================================
    else if (strcmp(str3, argv[1]) == 0) {

        fprintf(fichier_out, "***************************************************************\n");
        fprintf(fichier_out, "*****            Sphere SHS angular distribution          *****\n");
        fprintf(fichier_out, " ****************OUTPUT****************************************\n");
        fprintf(fichier_out, " grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n");

        fichier_grace1 = fopen("out_plot", "w");
        fprintf(fichier_grace1, " grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n");

        printf("Enter the number of points between 0 and Pi for the angular SHS analysis: ");
        scanf("%d", &pas2);
        mmmf = 2 * pas2 + 1;

        printf("********** Begin Sphere SHS Angular Distribution calculation **************\n");

        for (mmm = 0; mmm < mmmf; mmm++) {
            grandtheta = -pi + mmm * (pi / pas2);
            grandthetadegree = grandtheta * 360.0 / (2.0 * pi);
            printf("angle value: %.1f °\n", grandthetadegree);

            hx = (2.0 * pi) / (2.0 * pas);
            hy = pi / pas;
            nx = 2 * pas;
            ny = pas;

            grandgammazxx_moy = grandgammayyy_moy = grandgammazyy_moy = 0.0 + 0.0*I;
            grandgammaxxx_moy = grandgammayxx_moy = 0.0 + 0.0*I;
            grandgammaxyy_moy = 0.0 + 0.0*I;

            for (i = 0; i < nx; i++) {
                for (j = 0; j < ny; j++) {
                    xi = hx/2.0 + i * hx;
                    yj = hy/2.0 + j * hy;

                    grandgammazxx_moy += hx*hy * grandgammazxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammaxxx_moy += hx*hy * grandgammaxxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammayxx_moy += hx*hy * grandgammayxx(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammaxyy_moy += hx*hy * grandgammaxyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammazyy_moy += hx*hy * grandgammazyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                    grandgammayyy_moy += hx*hy * grandgammayyy(R, xi, yj, grandtheta, beta, petitk, grandk);
                }
            }

            IV_0  = creal(grandgammaxxx_moy * conj(grandgammaxxx_moy));
            IV_90 = creal(grandgammaxyy_moy * conj(grandgammaxyy_moy));

            IH_0 = creal(grandgammayxx_moy*conj(grandgammayxx_moy)*cos(grandtheta)*cos(grandtheta)
                       + sin(grandtheta)*sin(grandtheta)*grandgammazxx_moy*conj(grandgammazxx_moy)
                       - cos(grandtheta)*sin(grandtheta)*(grandgammayxx_moy*conj(grandgammazxx_moy)
                                                         + grandgammazxx_moy*conj(grandgammayxx_moy)));

            IH_90 = creal(grandgammayyy_moy*conj(grandgammayyy_moy)*cos(grandtheta)*cos(grandtheta)
                        + sin(grandtheta)*sin(grandtheta)*grandgammazyy_moy*conj(grandgammazyy_moy)
                        - cos(grandtheta)*sin(grandtheta)*(grandgammayyy_moy*conj(grandgammazyy_moy)
                                                          + grandgammazyy_moy*conj(grandgammayyy_moy)));

            printf("IH_0=%lf    IH_90=%lf   IV_0=%lf    IV_90=%lf\n", IH_0, IH_90, IV_0, IV_90);
            fprintf(fichier_out, " %lf    %lf   %lf   %lf   %lf\n", grandthetadegree, IH_0, IH_90, IV_0, IV_90);
            fprintf(fichier_grace1, " %lf    %lf   %lf   %lf   %lf\n", grandthetadegree, IH_0, IH_90, IV_0, IV_90);
        }

        printf("-----------------------------------------------------------------------------\n");
        printf("calculation finished, open the out_plot file with gnuplot\n");
        printf("-----------------------------------------------------------------------------\n");
        fclose(fichier_grace1);
    }

    else {
        printf("Unknown keyword. Use: polarplot_single | polarplot_integrate | angle_scattering\n");
    }

    fclose(fichier_out);
    return 0;
}