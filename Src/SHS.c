/*********************************************************************************
******************PROGRAM Py_SHS***************************************************
**********************************************************************************
   Py_SHS - An open source software about Second Harmonic Scattering 
   developed at Institut Charles Gerhardt Montpellier - ENSCM
   Dr Pierre-Marie GASSIN
   Dr Gassin GASSIN
    (june 2021)  

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

*******************************************************************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<string.h> 

// CONVENTION ZYZ pour matrice Euler!!!

double Txx(double phi,double theta,double psi){
double uu; 
uu=-sin(phi)*sin(psi)+cos(theta)*cos(phi)*cos(psi);
return uu;
} 

double Txy(double phi,double theta,double psi){
double uu; 
uu=-sin(phi)*cos(psi)-cos(theta)*cos(phi)*sin(psi);
return uu;
} 

double Txz(double phi,double theta,double psi){
double uu; 
uu=sin(theta)*cos(phi);
return uu;
} 

double Tyx(double phi,double theta,double psi){
double uu; 
uu=cos(phi)*sin(psi)+cos(theta)*sin(phi)*cos(psi);
return uu;
} 

double Tyy(double phi,double theta,double psi){
double uu; 
uu=cos(phi)*cos(psi)-cos(theta)*sin(psi)*sin(phi);
return uu;
} 

double Tyz(double phi,double theta,double psi){
double uu; 
uu=sin(theta)*sin(phi);
return uu;
} 

double Tzx(double phi,double theta,double psi){
double uu; 
uu=-sin(theta)*cos(psi);
return uu;
} 

double Tzy(double phi,double theta,double psi){
double uu; 
uu=sin(theta)*sin(psi);
return uu;
} 

double Tzz(double phi,double theta,double psi){
double uu; 
uu=cos(theta);
return uu;
} 

double T(double phi,double theta,double psi,int i,int j){
double uu;
if (i==0 && j==0){
uu=Txx(phi,theta,psi);
}
if (i==0 && j==1){
uu=Txy(phi,theta,psi);
}
if (i==0 && j==2){
uu=Txz(phi,theta,psi);
}
if (i==1 && j==0){
uu=Tyx(phi,theta,psi);
}
if (i==1 && j==1){
uu=Tyy(phi,theta,psi);
}
if (i==1 && j==2){
uu=Tyz(phi,theta,psi);
}
if (i==2 && j==0){
uu=Tzx(phi,theta,psi);
}
if (i==2 && j==1){
uu=Tzy(phi,theta,psi);
}
if (i==2 && j==2){
uu=Tzz(phi,theta,psi);
}
return uu;
}

double betalabijk(double phi,double theta,double psi,int i,int j,int k,double beta[3][3][3]){
double uuu;
uuu=(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[0][0][0]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[0][0][1]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[0][0][2]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[0][1][0]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[0][1][1]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[0][1][2]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[0][2][0]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[0][2][1]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[0][2][2]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[1][0][0]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[1][0][1]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[1][0][2]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[1][1][0]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[1][1][1]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[1][1][2]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[1][2][0]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[1][2][1]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[1][2][2]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[2][0][0]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[2][0][1]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[2][0][2]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[2][1][0]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[2][1][1]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[2][1][2]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[2][2][0]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[2][2][1]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[2][2][2]);
    return uuu;
}


double xyzlabijk(double phi, double theta,double psi, int i, double x, double y,double z){
// c'est i qui dit si c'est x (i=0) y(i=1) ou z (i=2)...
double uu;
uu=T(phi,theta,psi,i,0)*x+T(phi,theta,psi,i,1)*y +T(phi,theta,psi,i,2)*z;
return uu;
}

double complex betalabzyy_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk, double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,2,1,1,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabxxx_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double  complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,0,0,0,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabzxx_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,2,0,0,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabxyy_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,0,1,1,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabxxy_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,0,0,1,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabzxy_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,2,0,1,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabyyy_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,1,1,1,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}

double complex betalabyxx_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,1,0,0,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}


double complex betalabyxy_retard(double phi,double theta,double psi,double beta[3][3][3],double x,double y,double z,double complex petitk,double complex grandk,double grandtheta){
double complex uuuvv;
    uuuvv=betalabijk(phi,theta,psi,1,0,1,beta)*(cexp((xyzlabijk(phi,theta,psi,2,x,y,z)*2*petitk-xyzlabijk(phi,theta,psi,1,x,y,z)*grandk*sin(grandtheta)-xyzlabijk(phi,theta,psi,2,x,y,z)*grandk*cos(grandtheta))*I));
    return uuuvv;}



int main(int argc, char *argv[])
{
/* Sortie ligne par ligne (pour suivi en direct depuis la GUI) */
setvbuf(stdout, NULL, _IOLBF, 0);
setvbuf(stderr, NULL, _IONBF, 0);

FILE* fichier_in=NULL;
FILE* fichier_beta=NULL;
FILE* fichier_out=NULL;
FILE* fichier_grace1=NULL;
FILE* fichier_grace2=NULL;
char chaine[121];
char nome[20]; 
int nombredipole=0;
int* noomdip=NULL;
printf("********************************************************\n");
printf("****           Py_SHS: An open source software      ****\n");
printf("****         about Second Harmonic Scattering       ****\n");
printf("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n");          
printf("****                  ICGM - ENSCM                  ****\n"); 
printf("****                   june 2021                    ****\n");
printf("********************************************************\n");
printf("                                                        \n");
printf("           Don't forget to cite this work:              \n");
printf("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n");
printf("                                                        \n");
printf("********************************************************\n");
printf("********************************************************\n");
printf("********       Enter the input parameters     **********\n");
printf("********************************************************\n");
printf("Number of dipoles ? ");
scanf("%d", &nombredipole);
printf("The number of dipoles is: %d \n",nombredipole);
if (nombredipole > 0) 
    {
        noomdip = malloc(nombredipole * sizeof(int)); // On alloue de la mémoire pour le tableau
        if (noomdip == NULL) // On vérifie si l'allocation a marché ou non
        {
         fprintf(stderr, "Échec allocation : %zu octets demandés\n", nombredipole * sizeof(int));
            exit(0); // On arrête tout
        }


int i,j,k,kkk,ii,iii,ic;
int nom;
double complex grandk,petitk;
double pi=3.1415926535897;
double pisur2;
double beta[3][3][3],betameso[3][3][3];
double nn[nombredipole],phi[nombredipole],theta[nombredipole],psi[nombredipole],x[nombredipole],y[nombredipole],z[nombredipole];
double av=0,bv=0,cv=0,ah=0,bh=0,ch=0,bv1=0,bv2=0,bh1=0,bh2=0;
double betamesoxxx[nombredipole],betamesoxxy[nombredipole],betamesoxxz[nombredipole],betamesoxyx[nombredipole],betamesoxyy[nombredipole],betamesoxyz[nombredipole],betamesoxzx[nombredipole],
betamesoxzy[nombredipole],betamesoxzz[nombredipole];
double betamesoyxx[nombredipole],betamesoyxy[nombredipole],betamesoyxz[nombredipole],betamesoyyx[nombredipole],betamesoyyy[nombredipole],betamesoyyz[nombredipole],betamesoyzx[nombredipole],
betamesoyzy[nombredipole],betamesoyzz[nombredipole];
double betamesozxx[nombredipole],betamesozxy[nombredipole],betamesozxz[nombredipole],betamesozyx[nombredipole],betamesozyy[nombredipole],betamesozyz[nombredipole],betamesozzx[nombredipole],
betamesozzy[nombredipole],betamesozzz[nombredipole];
short nb_lignes_lues, nb_val_lues;
int pas,nx,nz,ny;
double hx,hy,hz,xi,yj,zk;
double compteur;
double complex betalabozyyint,betalaboxxxint,betalabozxxint,betalaboxyyint,betalaboxxyint,betalabozxyint,betalaboyyyint,betalaboyxyint,betalaboyxxint;
double I4v,I2v,I2h,I4h;
double Iv,Ih,gammaaa,nomc,noms,deltax,deltay,deltaz;
char str1[20]="polarplot_single";
char str2[20]="polarplot_integrate";
char str3[20]="angle_scattering";
double grandtheta,grandthetadegree,grandtheta_1,grandtheta_2,grandthetadegree_1,grandthetadegree_2,IH_0,IH_90,IV_0,IV_90;
int pas2,mmm,mmmf,mmmfbis,point;
double beta_y_moy_0,beta_z_moy_0,beta_y_moy_90,beta_z_moy_90,beta_xxx_moy,beta_xyy_moy;
double complex n=1.33+0*I;
double realpartn,imagpartn;
double longueuronde=800;
pisur2=pi/2;
pas2=10;
pas=10;
nx=2*pas;
ny=pas;
nz=2*pas;
hx = (2*pi)/(nx);
hy = (pi)/(ny);
hz = (2*pi)/(nz);

printf("Enter the wavelength in nm? ");
scanf("%lf", &longueuronde);
printf("The wavelength is: %lf \n",longueuronde);

printf("Enter the real part of the refractive index? ");
scanf("%lf", &realpartn);

printf("Enter the imaginary part of the refractive index? ");
scanf("%lf", &imagpartn);

n=realpartn+imagpartn*I;
printf("The refractive index is: %lf + I*%lf \n",creal(n),cimag(n));

if( argc == 5 ) {
      printf("The input beta file is %s, the input orientation file is %s and the output file is %s \n", argv[2], argv[3], argv[4]);
   }
   else if( argc > 5 ) {
      printf("Too many arguments supplied.\n");
   }
   else {
      printf("Computation failed: \n");
      printf("4 arguments expected : the keyword of the calculation, the beta and orientation inputfiles and the outputfile. \n");
      printf("For more details see the documentation \n");
   }

fichier_beta=fopen(argv[2],"rt");
if (fichier_beta==NULL)
    {
        puts("Problem to open the input beta file \n!");
        exit(0) ;
    }

fichier_in = fopen(argv[3], "rt");
 if (fichier_in==NULL)
    {
        puts("Problem to open the input orientation file \n!");
        exit(0) ;
    }
fichier_out = fopen(argv[4], "w");
   if (fichier_out==NULL){
	printf("problem to open the output file to write results!\n");
}

grandk=4*pi*n/longueuronde;
petitk=2*pi*n/longueuronde;

fgets (chaine, 121, fichier_beta);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][0][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][0][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][0][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][1][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][1][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][1][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][2][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][2][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][2][2]);

fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][0][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][0][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][0][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][1][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][1][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][1][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][2][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][2][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][2][2]);

fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][0][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][0][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][0][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][1][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][1][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][1][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][2][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][2][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][2][2]);

fclose(fichier_beta);


fichier_grace2 = fopen("out_arrows_gnuplot", "w");
   if (fichier_grace2==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire fleche gnuplot!");
}

for (k=0;k<nombredipole;k++){	
	fgets (chaine, 121, fichier_in);
	nb_lignes_lues++ ;
	nb_val_lues = sscanf (chaine, "%lf%lf%lf%lf%lf%lf", &phi[k],&theta[k],&psi[k],&x[k],&y[k],&z[k]);
	deltax=cos(phi[k])*sin(theta[k]);
	deltay=sin(phi[k])*sin(theta[k]);
	deltaz=cos(theta[k]);
	fprintf(fichier_grace2,"%lf %lf %lf %lf %lf %lf \n",x[k],y[k],z[k],deltax,deltay,deltaz);

	betamesoxxx[k]=betalabijk(phi[k],theta[k],psi[k],0,0,0,beta);
	betamesoxxy[k]=betalabijk(phi[k],theta[k],psi[k],0,0,1,beta);
	betamesoxxz[k]=betalabijk(phi[k],theta[k],psi[k],0,0,2,beta);
	betamesoxyx[k]=betalabijk(phi[k],theta[k],psi[k],0,1,0,beta);
	betamesoxyy[k]=betalabijk(phi[k],theta[k],psi[k],0,1,1,beta);
	betamesoxyz[k]=betalabijk(phi[k],theta[k],psi[k],0,1,2,beta);
	betamesoxzx[k]=betalabijk(phi[k],theta[k],psi[k],0,2,0,beta);
	betamesoxzy[k]=betalabijk(phi[k],theta[k],psi[k],0,2,1,beta);
	betamesoxzz[k]=betalabijk(phi[k],theta[k],psi[k],0,2,2,beta);

	betamesoyxx[k]=betalabijk(phi[k],theta[k],psi[k],1,0,0,beta);
	betamesoyxy[k]=betalabijk(phi[k],theta[k],psi[k],1,0,1,beta);
	betamesoyxz[k]=betalabijk(phi[k],theta[k],psi[k],1,0,2,beta);
	betamesoyyx[k]=betalabijk(phi[k],theta[k],psi[k],1,1,0,beta);
	betamesoyyy[k]=betalabijk(phi[k],theta[k],psi[k],1,1,1,beta);
	betamesoyyz[k]=betalabijk(phi[k],theta[k],psi[k],1,1,2,beta);
	betamesoyzx[k]=betalabijk(phi[k],theta[k],psi[k],1,2,0,beta);
	betamesoyzy[k]=betalabijk(phi[k],theta[k],psi[k],1,2,1,beta);
	betamesoyzz[k]=betalabijk(phi[k],theta[k],psi[k],1,2,2,beta);

	betamesozxx[k]=betalabijk(phi[k],theta[k],psi[k],2,0,0,beta);
	betamesozxy[k]=betalabijk(phi[k],theta[k],psi[k],2,0,1,beta);
	betamesozxz[k]=betalabijk(phi[k],theta[k],psi[k],2,0,2,beta);
	betamesozyx[k]=betalabijk(phi[k],theta[k],psi[k],2,1,0,beta);
	betamesozyy[k]=betalabijk(phi[k],theta[k],psi[k],2,1,1,beta);
	betamesozyz[k]=betalabijk(phi[k],theta[k],psi[k],2,1,2,beta);
	betamesozzx[k]=betalabijk(phi[k],theta[k],psi[k],2,2,0,beta);
	betamesozzy[k]=betalabijk(phi[k],theta[k],psi[k],2,2,1,beta);
	betamesozzz[k]=betalabijk(phi[k],theta[k],psi[k],2,2,2,beta);
	}
fclose(fichier_in);
fclose(fichier_grace2);




fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"****           Py_SHS: An open source software      ****\n");
fprintf(fichier_out,"****         about Second Harmonic Scattering       ****\n");
fprintf(fichier_out,"**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n");          
fprintf(fichier_out,"****                  ICGM - ENSCM                  ****\n"); 
fprintf(fichier_out,"****                   june 2021                    ****\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"                                                        \n");
fprintf(fichier_out,"           Don't forget to cite this work:              \n");
fprintf(fichier_out,"    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n");
fprintf(fichier_out,"                                                        \n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"******          SIMULATION RESULTS               *******\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"******     SURFACE HYPERPOLARISABILITY TENSOR     ******\n");
fprintf(fichier_out,"\n");
fprintf(fichier_out,"    xxx=%lf   xxy=%lf   xxz= %lf            \n",beta[0][0][0],beta[0][0][1],beta[0][0][2]);
fprintf(fichier_out,"    xyx=%lf   xyy=%lf   xyz= %lf  \n",beta[0][1][0],beta[0][1][1],beta[0][1][2]);
fprintf(fichier_out,"    xzx=%lf   xzy=%lf   xzz= %lf  \n",beta[0][2][0],beta[0][2][1],beta[0][2][2]);
fprintf(fichier_out,"\n");
fprintf(fichier_out,"    yxx=%lf   yxy=%lf   yxz= %lf  \n",beta[1][0][0],beta[1][0][1],beta[1][0][2]);
fprintf(fichier_out,"    yyx=%lf   yyy=%lf   yyz= %lf  \n",beta[1][1][0],beta[1][1][1],beta[1][1][2]);
fprintf(fichier_out,"    yzx=%lf   yzy=%lf   yzz= %lf  \n",beta[1][2][0],beta[1][2][1],beta[1][2][2]);
fprintf(fichier_out,"\n");
fprintf(fichier_out,"    zxx=%lf   zxy=%lf   zxz= %lf  \n",beta[2][0][0],beta[2][0][1],beta[2][0][2]);
fprintf(fichier_out,"    zyx=%lf   zyy=%lf   zyz= %lf  \n",beta[2][1][0],beta[2][1][1],beta[2][1][2]);
fprintf(fichier_out,"    zzx=%lf   zzy=%lf   zzz= %lf  \n",beta[2][2][0],beta[2][2][1],beta[2][2][2]);
fprintf(fichier_out,"\n");
fprintf(fichier_out,"*******************************************************\n");
if (strcmp(str1,argv[1]) == 0){
printf("Enter the scattered angle in degree (0°=transmission) ? ");
scanf("%lf", &grandthetadegree);
printf("      ******************************     \n");
grandtheta=grandthetadegree*2*pi/360;
printf("********** Begin calculation polar plot at  %lf ° ********\n", grandthetadegree);
//printf("calculation in progress:");
for (i=0;i<nx;i++){
	compteur=i*100/nx;
	printf("calculation in progress : %lf pourcent \n",compteur);
	for (j=0;j<ny;j++){
		for (kkk=0;kkk<nz;kkk++){
            		xi = 0 + hx/2 + i*hx;
            		yj = 0 + hy/2 + j*hy;
            		zk = 0 + hz/2 + kkk*hz;
            		betalabozyyint=0+0*I;
            		betalaboxxxint=0+0*I;
           		    betalabozxxint=0+0*I;
            		betalaboxyyint=0+0*I;
            		betalaboxxyint=0+0*I;
            		betalabozxyint=0+0*I;
            		betalaboyxyint=0+0*I;
            		betalaboyxxint=0+0*I;
            		betalaboyyyint=0+0*I;
            		
			for (ii=0;ii<nombredipole;ii++){
				betameso[0][0][0]=betamesoxxx[ii];
				betameso[0][0][1]=betamesoxxy[ii];	
				betameso[0][0][2]=betamesoxxz[ii];
				betameso[0][1][0]=betamesoxyx[ii];
				betameso[0][1][1]=betamesoxyy[ii];	
				betameso[0][1][2]=betamesoxyz[ii];
				betameso[0][2][0]=betamesoxzx[ii];
				betameso[0][2][1]=betamesoxzy[ii];	
				betameso[0][2][2]=betamesoxzz[ii];
				
				betameso[1][0][0]=betamesoyxx[ii];
				betameso[1][0][1]=betamesoyxy[ii];	
				betameso[1][0][2]=betamesoyxz[ii];
				betameso[1][1][0]=betamesoyyx[ii];
				betameso[1][1][1]=betamesoyyy[ii];	
				betameso[1][1][2]=betamesoyyz[ii];
				betameso[1][2][0]=betamesoyzx[ii];
				betameso[1][2][1]=betamesoyzy[ii];	
				betameso[1][2][2]=betamesoyzz[ii];

				betameso[2][0][0]=betamesozxx[ii];
				betameso[2][0][1]=betamesozxy[ii];	
				betameso[2][0][2]=betamesozxz[ii];
				betameso[2][1][0]=betamesozyx[ii];
				betameso[2][1][1]=betamesozyy[ii];	
				betameso[2][1][2]=betamesozyz[ii];
				betameso[2][2][0]=betamesozzx[ii];
				betameso[2][2][1]=betamesozzy[ii];	
				betameso[2][2][2]=betamesozzz[ii];
				//printf("%lf \n",betameso[0][0][0]);
				
				
				
				betalabozyyint=betalabozyyint+betalabzyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalaboxxxint=betalaboxxxint+betalabxxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalabozxxint=betalabozxxint+betalabzxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalaboxyyint=betalaboxyyint+betalabxyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
               	betalaboxxyint=betalaboxxyint+betalabxxy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalabozxyint=betalabozxyint+betalabzxy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				betalaboyxyint=betalaboyxyint+betalabyxy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				betalaboyxxint=betalaboyxxint+betalabyxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				betalaboyyyint=betalaboyyyint+betalabyyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				//printf("%lf+%lfi\n",creal(betalabozyyint),cimag(betalabozyyint));
				}
					av=av+hx*hy*hz*(betalaboxxxint*conj(betalaboxxxint))*sin(yj);
            		cv=cv+hx*hy*hz*(betalaboxyyint*conj(betalaboxyyint))*sin(yj);
            		bv=bv+hx*hy*hz*sin(yj)*(4*(betalaboxxyint*conj(betalaboxxyint))+betalaboxxxint*(conj(betalaboxyyint))+(betalaboxyyint*(conj(betalaboxxxint))));
            		bv1=bv1+hx*hy*hz*sin(yj)*(2*(betalaboxxxint*(conj(betalaboxxyint)))+2*(betalaboxxyint*(conj(betalaboxxxint))));
            		bv2=bv2+hx*hy*hz*sin(yj)*(2*(betalaboxyyint*(conj(betalaboxxyint)))+2*(betalaboxxyint*(conj(betalaboxyyint))));
            				
            		ah=ah+hx*hy*hz*(betalabozxxint*conj(betalabozxxint)*sin(grandtheta)*sin(grandtheta)+betalaboyxxint*conj(betalaboyxxint)*cos(grandtheta)*cos(grandtheta)-sin(grandtheta)*cos(grandtheta)*(betalaboyxxint*conj(betalabozxxint)+betalabozxxint*conj(betalaboyxxint)))*sin(yj);			
            		ch=ch+hx*hy*hz*(betalaboyyyint*conj(betalaboyyyint)*cos(grandtheta)*cos(grandtheta)+betalabozyyint*conj(betalabozyyint)*sin(grandtheta)*sin(grandtheta)-sin(grandtheta)*cos(grandtheta)*(betalaboyyyint*conj(betalabozyyint)+betalabozyyint*conj(betalaboyyyint)))*sin(yj);
            		bh=bh+hx*hy*hz*sin(yj)*(cos(grandtheta)*cos(grandtheta)*(4*betalaboyxyint*conj(betalaboyxyint)+betalaboyxxint*conj(betalaboyyyint)+betalaboyyyint*conj(betalaboyxxint))+sin(grandtheta)*sin(grandtheta)*(4*betalabozxyint*conj(betalabozxyint)+betalabozxxint*conj(betalabozyyint)+betalabozyyint*conj(betalabozxxint))-sin(grandtheta)*cos(grandtheta)*(4*betalaboyxyint*conj(betalabozxyint)+betalaboyxxint*conj(betalabozyyint)+betalabozyyint*conj(betalaboyxxint)+4*betalabozxyint*conj(betalaboyxyint)+betalabozxxint*conj(betalaboyyyint)+betalaboyyyint*conj(betalabozxxint)));
            		bh1=bh1+hx*hy*hz*sin(yj)*(cos(grandtheta)*cos(grandtheta)*(2*betalaboyxxint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyxxint))+sin(grandtheta)*sin(grandtheta)*(2*betalabozxxint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozxxint))-sin(grandtheta)*cos(grandtheta)*(2*betalabozxxint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozxxint)+2*betalaboyxxint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyxxint)));
            		bh2=bh2+hx*hy*hz*sin(yj)*(cos(grandtheta)*cos(grandtheta)*(2*betalaboyyyint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyyyint))+sin(grandtheta)*sin(grandtheta)*(2*betalabozyyint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozyyint))-sin(grandtheta)*cos(grandtheta)*(2*betalabozyyint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozyyint)+2*betalaboyyyint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyyyint)));
            		//ATTENTION A FINIR D4ECRIRE LES EXPRESSIONS DE BH1 ET BH2!
            		//bh1=bh1+hx*hy*hz*sin(yj)*(2*(betalabozxxint*(conj(betalabozxyint)))+2*(betalabozxyint*(conj(betalabozxxint))));
           		    //bh2=bh2+hx*hy*hz*sin(yj)*(2*(betalabozyyint*(conj(betalabozxyint)))+2*(betalabozxyint*(conj(betalabozyyint))));
            		//ch=ch+hx*hy*hz*(betalabozyyint*conj(betalabozyyint))*sin(yj);
            		//printf("%lf\n",av);
            					}
		}
	}
printf("\n"); 
printf("**********calculation complete ************\n");
av=av/(8*pi*pi);
bv=bv/(8*pi*pi);
cv=cv/(8*pi*pi);
ah=ah/(8*pi*pi);
bh=bh/(8*pi*pi);
ch=ch/(8*pi*pi);           
printf("av = %lf\n",av);
printf("bv = %lf\n",bv);
printf("bv1 = %lf\n",bv1);
printf("bv2 = %lf\n",bv2);
printf("cv = %lf\n",cv);
printf("ah = %lf\n",ah);
printf("bh = %lf\n",bh);
printf("bh1 = %lf\n",bh1);
printf("bh2 = %lf\n",bh2);
printf("ch = %lf\n",ch);


printf("********************************************\n");   
I4v=(av-bv+cv)/((3*av+bv+3*cv));
I2v=(4)*(av-cv)/(3*av+bv+3*cv);
I4h=(ah-bh+ch)/((3*ah+bh+3*ch));
I2h=(4)*(ah-ch)/(3*ah+bh+3*ch);

printf("I2v = %lf \n",I2v);
printf("I4v = %lf \n",I4v);
printf("I2h = %lf \n",I2h);
printf("I4h = %lf \n",I4h);

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}

fprintf(fichier_out,"                                      \n");
fprintf(fichier_out," POLAR PLOT - CONFIGURATION AT %lf °  \n", grandthetadegree);
fprintf(fichier_out,"                                      \n");
fprintf(fichier_out,"nombre de dipole = %d  \n",nombredipole);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"I2v=%lf  I4v=%lf \n",I2v,I4v);
fprintf(fichier_out,"I2h=%lf  I4h=%lf \n",I2h,I4h);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"av=%lf   bv=%lf   bv1=%lf  bv2=%lf  cv= %lf \n",av,bv,bv1,bv2,cv);
fprintf(fichier_out,"ah=%lf   bh=%lf   bh1=%lf  bh2=%lf  ch= %lf \n",ah,bh,bh1,bh2,ch);

fclose(fichier_out);

fprintf(fichier_grace1,"gamma IV(gamma) IH(gamma)\n");
for (iii=0;iii<100;iii++){
gammaaa=2*pi*iii/100;
nomc=cos(gammaaa);
noms=sin(gammaaa);
Iv=av*(nomc*nomc*nomc*nomc)+bv*(nomc*nomc*noms*noms)+cv*(noms*noms*noms*noms)+bv1*nomc*nomc*nomc*noms+bv2*nomc*noms*noms*noms;
Ih=ah*(nomc*nomc*nomc*nomc)+bh*(nomc*nomc*noms*noms)+ch*(noms*noms*noms*noms)+bh1*nomc*nomc*nomc*noms+bh2*nomc*noms*noms*noms;
fprintf(fichier_grace1,"%lf %lf %lf\n",gammaaa,Iv,Ih);
}


fclose(fichier_grace1);
}

else if (strcmp(str2,argv[1]) == 0){
printf("Enter the beginning scattered angle in degree (0°=transmission) ? ");
scanf("%lf", &grandthetadegree_1);
printf("      ******************************     \n");
printf("Enter the final scattered angle in degree (0°=transmission) ? ");
scanf("%lf", &grandthetadegree_2);
printf("      ******************************     \n");
grandtheta_1=grandthetadegree_1*2*pi/360;
grandtheta_2=grandthetadegree_2*2*pi/360;
printf(" Number of integrale points ? ");
scanf("%d", &point);
printf("      ******************************     \n");


for (ic=0;ic<point;ic++){
grandtheta=grandtheta_1+ic*(grandtheta_2-grandtheta_1)/(point-1);
grandthetadegree=grandtheta*360/(2*pi);



for (i=0;i<nx;i++){
	compteur=(float)(ic*100/point)+((i*100/nx)/point);
	printf("calculation in progress : %lf pourcent \n",compteur);
	for (j=0;j<ny;j++){
		for (kkk=0;kkk<nz;kkk++){
            		xi = 0 + hx/2 + i*hx;
            		yj = 0 + hy/2 + j*hy;
            		zk = 0 + hz/2 + kkk*hz;
            		betalabozyyint=0+0*I;
            		betalaboxxxint=0+0*I;
           		    betalabozxxint=0+0*I;
            		betalaboxyyint=0+0*I;
            		betalaboxxyint=0+0*I;
            		betalabozxyint=0+0*I;
            		betalaboyxyint=0+0*I;
            		betalaboyxxint=0+0*I;
            		betalaboyyyint=0+0*I;
            		
			for (ii=0;ii<nombredipole;ii++){
				betameso[0][0][0]=betamesoxxx[ii];
				betameso[0][0][1]=betamesoxxy[ii];	
				betameso[0][0][2]=betamesoxxz[ii];
				betameso[0][1][0]=betamesoxyx[ii];
				betameso[0][1][1]=betamesoxyy[ii];	
				betameso[0][1][2]=betamesoxyz[ii];
				betameso[0][2][0]=betamesoxzx[ii];
				betameso[0][2][1]=betamesoxzy[ii];	
				betameso[0][2][2]=betamesoxzz[ii];
				
				betameso[1][0][0]=betamesoyxx[ii];
				betameso[1][0][1]=betamesoyxy[ii];	
				betameso[1][0][2]=betamesoyxz[ii];
				betameso[1][1][0]=betamesoyyx[ii];
				betameso[1][1][1]=betamesoyyy[ii];	
				betameso[1][1][2]=betamesoyyz[ii];
				betameso[1][2][0]=betamesoyzx[ii];
				betameso[1][2][1]=betamesoyzy[ii];	
				betameso[1][2][2]=betamesoyzz[ii];

				betameso[2][0][0]=betamesozxx[ii];
				betameso[2][0][1]=betamesozxy[ii];	
				betameso[2][0][2]=betamesozxz[ii];
				betameso[2][1][0]=betamesozyx[ii];
				betameso[2][1][1]=betamesozyy[ii];	
				betameso[2][1][2]=betamesozyz[ii];
				betameso[2][2][0]=betamesozzx[ii];
				betameso[2][2][1]=betamesozzy[ii];	
				betameso[2][2][2]=betamesozzz[ii];

				betalabozyyint=betalabozyyint+betalabzyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalaboxxxint=betalaboxxxint+betalabxxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalabozxxint=betalabozxxint+betalabzxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalaboxyyint=betalaboxyyint+betalabxyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
               	betalaboxxyint=betalaboxxyint+betalabxxy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
                betalabozxyint=betalabozxyint+betalabzxy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				betalaboyxyint=betalaboyxyint+betalabyxy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				betalaboyxxint=betalaboyxxint+betalabyxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				betalaboyyyint=betalaboyyyint+betalabyyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
				}
					av=av+hx*hy*hz*(betalaboxxxint*conj(betalaboxxxint))*sin(yj);
            		cv=cv+hx*hy*hz*(betalaboxyyint*conj(betalaboxyyint))*sin(yj);
            		bv=bv+hx*hy*hz*sin(yj)*(4*(betalaboxxyint*conj(betalaboxxyint))+betalaboxxxint*(conj(betalaboxyyint))+(betalaboxyyint*(conj(betalaboxxxint))));
            		bv1=bv1+hx*hy*hz*sin(yj)*(2*(betalaboxxxint*(conj(betalaboxxyint)))+2*(betalaboxxyint*(conj(betalaboxxxint))));
            		bv2=bv2+hx*hy*hz*sin(yj)*(2*(betalaboxyyint*(conj(betalaboxxyint)))+2*(betalaboxxyint*(conj(betalaboxyyint))));
            				
            		ah=ah+hx*hy*hz*(betalabozxxint*conj(betalabozxxint)*sin(grandtheta)*sin(grandtheta)+betalaboyxxint*conj(betalaboyxxint)*cos(grandtheta)*cos(grandtheta)-sin(grandtheta)*cos(grandtheta)*(betalaboyxxint*conj(betalabozxxint)+betalabozxxint*conj(betalaboyxxint)))*sin(yj);			
            		ch=ch+hx*hy*hz*(betalaboyyyint*conj(betalaboyyyint)*cos(grandtheta)*cos(grandtheta)+betalabozyyint*conj(betalabozyyint)*sin(grandtheta)*sin(grandtheta)-sin(grandtheta)*cos(grandtheta)*(betalaboyyyint*conj(betalabozyyint)+betalabozyyint*conj(betalaboyyyint)))*sin(yj);
            		bh=bh+hx*hy*hz*sin(yj)*(cos(grandtheta)*cos(grandtheta)*(4*betalaboyxyint*conj(betalaboyxyint)+betalaboyxxint*conj(betalaboyyyint)+betalaboyyyint*conj(betalaboyxxint))+sin(grandtheta)*sin(grandtheta)*(4*betalabozxyint*conj(betalabozxyint)+betalabozxxint*conj(betalabozyyint)+betalabozyyint*conj(betalabozxxint))-sin(grandtheta)*cos(grandtheta)*(4*betalaboyxyint*conj(betalabozxyint)+betalaboyxxint*conj(betalabozyyint)+betalabozyyint*conj(betalaboyxxint)+4*betalabozxyint*conj(betalaboyxyint)+betalabozxxint*conj(betalaboyyyint)+betalaboyyyint*conj(betalabozxxint)));
            		
            		bh1=bh1+hx*hy*hz*sin(yj)*(cos(grandtheta)*cos(grandtheta)*(2*betalaboyxxint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyxxint))+sin(grandtheta)*sin(grandtheta)*(2*betalabozxxint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozxxint))-sin(grandtheta)*cos(grandtheta)*(2*betalabozxxint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozxxint)+2*betalaboyxxint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyxxint)));
            		bh2=bh2+hx*hy*hz*sin(yj)*(cos(grandtheta)*cos(grandtheta)*(2*betalaboyyyint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyyyint))+sin(grandtheta)*sin(grandtheta)*(2*betalabozyyint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozyyint))-sin(grandtheta)*cos(grandtheta)*(2*betalabozyyint*conj(betalabozxyint)+2*betalabozxyint*conj(betalabozyyint)+2*betalaboyyyint*conj(betalaboyxyint)+2*betalaboyxyint*conj(betalaboyyyint)));
            		
            		//ATTENTION A FINIR D4ECRIRE LES EXPRESSIONS DE BH1 ET BH2!
            		//bh1=bh1+hx*hy*hz*sin(yj)*(2*(betalabozxxint*(conj(betalabozxyint)))+2*(betalabozxyint*(conj(betalabozxxint))));
           		    //bh2=bh2+hx*hy*hz*sin(yj)*(2*(betalabozyyint*(conj(betalabozxyint)))+2*(betalabozxyint*(conj(betalabozyyint))));
            		//ch=ch+hx*hy*hz*(betalabozyyint*conj(betalabozyyint))*sin(yj);
            					}
		}
	}
}	

av=av/point;
cv=cv/point;
bv=bv/point;
bv1=bv1/point;
bv2=bv2/point;
ah=ah/point;
ch=ch/point;
bh=bh/point;
bh1=bh1/point;
bh2=bh2/point;
	
printf("\n"); 
printf("**********calculation complete ************\n");   
av=av/(8*pi*pi);
bv=bv/(8*pi*pi);
cv=cv/(8*pi*pi);
ah=ah/(8*pi*pi);
bh=bh/(8*pi*pi);
ch=ch/(8*pi*pi);        
printf("av = %lf\n",av);
printf("bv = %lf\n",bv);
printf("bv1 = %lf\n",bv1);
printf("bv2 = %lf\n",bv2);
printf("cv = %lf\n",cv);
printf("ah = %lf\n",ah);
printf("bh = %lf\n",bh);
printf("bh1 = %lf\n",bh1);
printf("bh2 = %lf\n",bh2);
printf("ch = %lf\n",ch);


printf("********************************************\n");   
I4v=(av-bv+cv)/((3*av+bv+3*cv));
I2v=(4)*(av-cv)/(3*av+bv+3*cv);
I4h=(ah-bh+ch)/((3*ah+bh+3*ch));
I2h=(4)*(ah-ch)/(3*ah+bh+3*ch);

printf("I2v = %lf \n",I2v);
printf("I4v = %lf \n",I4v);
printf("I2h = %lf \n",I2h);
printf("I4h = %lf \n",I4h);

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}

fprintf(fichier_out,"                                      \n");
fprintf(fichier_out," POLAR PLOT - CONFIGURATION AT %lf °  \n", grandthetadegree);
fprintf(fichier_out,"                                      \n");
fprintf(fichier_out,"nombre de dipole = %d  \n",nombredipole);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"I2v=%lf  I4v=%lf \n",I2v,I4v);
fprintf(fichier_out,"I2h=%lf  I4h=%lf \n",I2h,I4h);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"av=%lf   bv=%lf   bv1=%lf  bv2=%lf  cv= %lf \n",av,bv,bv1,bv2,cv);
fprintf(fichier_out,"ah=%lf   bh=%lf   bh1=%lf  bh2=%lf  ch= %lf \n",ah,bh,bh1,bh2,ch);

fclose(fichier_out);

fprintf(fichier_grace1,"gamma IV(gamma) IH(gamma)\n");
for (iii=0;iii<100;iii++){
gammaaa=2*pi*iii/100;
nomc=cos(gammaaa);
noms=sin(gammaaa);
Iv=av*(nomc*nomc*nomc*nomc)+bv*(nomc*nomc*noms*noms)+cv*(noms*noms*noms*noms)+bv1*nomc*nomc*nomc*noms+bv2*nomc*noms*noms*noms;
Ih=ah*(nomc*nomc*nomc*nomc)+bh*(nomc*nomc*noms*noms)+ch*(noms*noms*noms*noms)+bh1*nomc*nomc*nomc*noms+bh2*nomc*noms*noms*noms;
fprintf(fichier_grace1,"%lf %lf %lf\n",gammaaa,Iv,Ih);
}


fclose(fichier_grace1);
}





else if (strcmp(str3,argv[1]) == 0){


fprintf(fichier_out,"                                      \n");
fprintf(fichier_out,"* SHS angular distribution            \n");
fprintf(fichier_out,"                                      \n");
fprintf(fichier_out,"nombre de dipole = %d  \n",nombredipole);
fprintf(fichier_out," ****************OUTPUT****************************************\n");
fprintf(fichier_out,  " grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n"); 

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}
fprintf(fichier_grace1," grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n");
printf("Enter the number of points between 0 and 180° for the angular SHS analysis ?"); 
scanf("%d", &pas2);
mmmf=2*(pas2)+1;
mmmfbis=pas2+1;

printf("********** Begin SHS Angular Distribution calculation**************\n");
for (mmm=0;mmm<mmmfbis;mmm++){
	grandtheta=-pi+mmm*(pi/pas2);
        grandthetadegree=(grandtheta*360)/(2*pi);
        printf("angle value: %lf \n",grandthetadegree);
        beta_y_moy_0=0;  //ah A CHANGER!!
        beta_z_moy_0=0;  //ch
        beta_y_moy_90=0;
        beta_z_moy_90=0;
        beta_xxx_moy=0;
        beta_xyy_moy=0;
	for (i=0;i<nx;i++){
		//compteur=i*100/nx;
	        //printf("calculation in progress: %lf \n",compteur);
		for (j=0;j<ny;j++){
			for (kkk=0;kkk<nz;kkk++){
                xi = 0 + hx/2 + i*hx;
                yj = 0 + hy/2 + j*hy;
                zk = 0 + hz/2 + kkk*hz;
                    
                
				betalabozyyint=0+0*I;
            	betalabozxxint=0+0*I;
           		betalaboyxxint=0+0*I;
            	betalaboyyyint=0+0*I;
            	betalaboxyyint=0+0*I;
            	betalaboxxxint=0+0*I;
            	
		for (ii=0;ii<nombredipole;ii++){
			betameso[0][0][0]=betamesoxxx[ii];
			betameso[0][0][1]=betamesoxxy[ii];	
			betameso[0][0][2]=betamesoxxz[ii];
			betameso[0][1][0]=betamesoxyx[ii];
			betameso[0][1][1]=betamesoxyy[ii];	
			betameso[0][1][2]=betamesoxyz[ii];
			betameso[0][2][0]=betamesoxzx[ii];
			betameso[0][2][1]=betamesoxzy[ii];	
			betameso[0][2][2]=betamesoxzz[ii];
				
			betameso[1][0][0]=betamesoyxx[ii];
			betameso[1][0][1]=betamesoyxy[ii];	
			betameso[1][0][2]=betamesoyxz[ii];
			betameso[1][1][0]=betamesoyyx[ii];
			betameso[1][1][1]=betamesoyyy[ii];	
			betameso[1][1][2]=betamesoyyz[ii];
			betameso[1][2][0]=betamesoyzx[ii];
			betameso[1][2][1]=betamesoyzy[ii];	
			betameso[1][2][2]=betamesoyzz[ii];

			betameso[2][0][0]=betamesozxx[ii];
			betameso[2][0][1]=betamesozxy[ii];	
			betameso[2][0][2]=betamesozxz[ii];
			betameso[2][1][0]=betamesozyx[ii];
			betameso[2][1][1]=betamesozyy[ii];	
			betameso[2][1][2]=betamesozyz[ii];
			betameso[2][2][0]=betamesozzx[ii];
			betameso[2][2][1]=betamesozzy[ii];	
			betameso[2][2][2]=betamesozzz[ii];

			betalaboyyyint=betalaboyyyint+betalabyyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
            betalaboxxxint=betalaboxxxint+betalabxxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
            betalabozxxint=betalabozxxint+betalabzxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
            betalabozyyint=betalabozyyint+betalabzyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
            betalaboxyyint=betalaboxyyint+betalabxyy_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
			betalaboyxxint=betalaboyxxint+betalabyxx_retard(xi, yj, zk, betameso, x[ii],y[ii],z[ii],petitk, grandk,grandtheta);
			
				}
				beta_xxx_moy=beta_xxx_moy+hx*hy*hz*(betalaboxxxint*conj(betalaboxxxint))*sin(yj);
                beta_xyy_moy=beta_xyy_moy+hx*hy*hz*(betalaboxyyint*conj(betalaboxyyint))*sin(yj);
                beta_y_moy_0=beta_y_moy_0+hx*hy*hz*(betalaboyxxint*cos(grandtheta)*cos(grandtheta)*conj(betalaboyxxint)+sin(grandtheta)*sin(grandtheta)*betalabozxxint*conj(betalabozxxint)-cos(grandtheta)*sin(grandtheta)*(betalabozxxint*conj(betalaboyxxint)+betalaboyxxint*conj(betalabozxxint)))*sin(yj);
                beta_z_moy_0=beta_z_moy_0+hx*hy*hz*(betalaboyyyint*conj(betalaboyyyint)*cos(grandtheta)*cos(grandtheta)+betalabozyyint*conj(betalabozyyint)*sin(grandtheta)*sin(grandtheta)-sin(grandtheta)*cos(grandtheta)*(betalaboyyyint*conj(betalabozyyint)+betalabozyyint*conj(betalaboyyyint)))*sin(yj);//hx*hy*hz*(betalaboyyyint*sin(grandtheta)*sin(grandtheta)*conj(betalaboyyyint)+cos(grandtheta)*cos(grandtheta)*betalabozyyint*conj(betalabozyyint)-cos(grandtheta)*sin(grandtheta)*(betalabozyyint*conj(betalaboyyyint)+betalaboyyyint*conj(betalabozyyint)))*sin(yj);
                
			}
		}
	}
    
    
    IH_0=beta_y_moy_0/(8*pi*pi);
    IH_90=beta_z_moy_0/(8*pi*pi);
    IV_0=beta_xxx_moy/(8*pi*pi);
    IV_90=beta_xyy_moy/(8*pi*pi);
    printf("IH_0=%lf    IH_90=%lf   IV_0=%lf    IV_90=%lf \n",IH_0,IH_90,IV_0,IV_90);
    fprintf(fichier_out," %lf    %lf   %lf   %lf   %lf \n",grandthetadegree,IH_0,IH_90,IV_0,IV_90);
    fprintf(fichier_grace1," %lf    %lf   %lf   %lf   %lf \n",grandthetadegree,IH_0,IH_90,IV_0,IV_90);
}




printf("calculation finished \n");
fclose(fichier_out);
fclose(fichier_grace1);


}

else {
printf("--------------------- ----------------------------------------------------\n"); 
printf("             ERROR in the name of the keyword                             \n"); 
printf(" USE ONLY polarplot_single OR polarplot_integrate OR angle_scattering  !! \n"); 
printf(" ! no calculation has been done ! Retry with the correct keyword          \n");
printf("--------------------------------------------------------------------------\n"); 
}


free(noomdip);
    }


}
