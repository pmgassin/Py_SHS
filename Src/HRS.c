/*********************************************************************************
******************PROGRAM Py_SHS***************************************************
**********************************************************************************
   Py_SHS - An open source software about Second Harmonic Scattering 
   developed at Institut Charles Gerhardt Montpellier - ENSCM
   Dr Pierre-Marie GASSIN
   Dr Gassin GASSIN
    (september 2023)  

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

double betalabzyy_sin(double phi,double theta,double psi,double beta[3][3][3]){
double uuuvv;
    uuuvv=betalabijk(phi,theta,psi,2,1,1,beta)*betalabijk(phi,theta,psi,2,1,1,beta)*sin(theta);
    return uuuvv;}

double betalabxxx_sin(double phi,double theta,double psi,double beta[3][3][3]){
double uuuvv;
    uuuvv=betalabijk(phi,theta,psi,0,0,0,beta)*betalabijk(phi,theta,psi,0,0,0,beta)*sin(theta);
    return uuuvv;}

double betalabzxx_sin(double phi,double theta,double psi,double beta[3][3][3]){
double uuuvv;
    uuuvv=betalabijk(phi,theta,psi,2,0,0,beta)*betalabijk(phi,theta,psi,2,0,0,beta)*sin(theta);
    return uuuvv;}

double betalabxyy_sin(double phi,double theta,double psi,double beta[3][3][3]){
double uuuvv;
    uuuvv=betalabijk(phi,theta,psi,0,1,1,beta)*betalabijk(phi,theta,psi,0,1,1,beta)*sin(theta);
return uuuvv;}

double betalabyxx_sin(double phi,double theta,double psi,double beta[3][3][3]){
double uuuvv;
    uuuvv=betalabijk(phi,theta,psi,1,0,0,beta)*betalabijk(phi,theta,psi,1,0,0,beta)*sin(theta);
    return uuuvv;}

double betalabyyy_sin(double phi,double theta,double psi,double beta[3][3][3]){
double uuuvv;
    uuuvv=betalabijk(phi,theta,psi,1,1,1,beta)*betalabijk(phi,theta,psi,1,1,1,beta)*sin(theta);
return uuuvv;}

double beta_bv(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=(4*(betalabijk(phi,theta,psi,0,0,1,beta)*betalabijk(phi,theta,psi,0,0,1,beta))+2*betalabijk(phi,theta,psi,0,0,0,beta)*betalabijk(phi,theta,psi,0,1,1,beta))*sin(theta);
    return uuu;}

double beta_bh_90(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=(4*(betalabijk(phi,theta,psi,2,0,1,beta)*betalabijk(phi,theta,psi,2,0,1,beta))+2*betalabijk(phi,theta,psi,2,0,0,beta)*betalabijk(phi,theta,psi,2,1,1,beta))*sin(theta);
    return uuu;}

double beta_bh_0(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=(4*(betalabijk(phi,theta,psi,1,0,1,beta)*betalabijk(phi,theta,psi,1,0,1,beta))+2*betalabijk(phi,theta,psi,1,0,0,beta)*betalabijk(phi,theta,psi,1,1,1,beta))*sin(theta);
    return uuu;}


double betalaboyxxint(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=betalabijk(phi,theta,psi,1,0,0,beta);
    return uuu;}
    
double betalaboyyyint(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=betalabijk(phi,theta,psi,1,1,1,beta);
    return uuu;}

double betalabozxxint(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=betalabijk(phi,theta,psi,2,0,0,beta);
    return uuu;} 
      
double betalabozyyint(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=betalabijk(phi,theta,psi,2,1,1,beta);
    return uuu;} 

double betalaboxxxint(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=betalabijk(phi,theta,psi,0,0,0,beta);
    return uuu;} 
      
double betalaboxyyint(double phi,double theta,double psi,double beta[3][3][3]){
double   uuu;
uuu=betalabijk(phi,theta,psi,0,1,1,beta);
    return uuu;} 


int main(int argc, char *argv[])
{
/* Sortie ligne par ligne (pour suivi en direct depuis la GUI) */
setvbuf(stdout, NULL, _IOLBF, 0);
setvbuf(stderr, NULL, _IONBF, 0);

FILE* fichier_beta=NULL;
FILE* fichier_out=NULL;
FILE* fichier_grace1=NULL;
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
int i,j,k,kkk,ii,iii;
int nom;
double pi=3.1415926535897;
double pisur2;
double beta[3][3][3];
double phi,theta,psi;
double av=0,bv=0,cv=0,ah=0,bh=0,ch=0,bv1=0,bv2=0,bh1=0,bh2=0;
short nb_lignes_lues, nb_val_lues;
int pas,nx,nz,ny;
double hx,hy,hz,xi,yj,zk;
double compteur;
double I4v,I2v,I2h,I4h;
double Iv,Ih,gammaaa,nomc,noms,deltax,deltay,deltaz;
char str1[20]="polarplot_single";
char str2[20]="polarplot_integrate";
char str3[20]="angle_scattering";
double grandtheta,grandthetadegree,grandtheta_1,grandtheta_2,grandthetadegree_1,grandthetadegree_2,IH_0,IH_90,IV_0,IV_90,Iz_0,Iz_90,Iy_90,Iy_0,IV_0b,IV_90b;
int pas2,mmm,mmmf,mmmfbis,point,ic;
double ui1,ui2,ui3,ui4;
double beta_xxx_moy,beta_xyy_moy,beta_y_moy_0,beta_y_moy_90,beta_z_moy_0,beta_z_moy_90=0;
pisur2=pi/2;
double Dep,betaHRS,betaHRS_esu;
//pas2=10;
pas=15;
pas2=30;
nx=2*pas;
ny=nx;
nz=nx;
hx = (2*pi)/(nx);
hy = (pi)/(ny);
hz = (2*pi)/(nz);

if( argc == 4 ) {
      printf("The input beta file is %s, and the output file is %s \n", argv[2], argv[3]);
   }
   else if( argc > 4 ) {
      printf("Too many arguments supplied.\n");
   }
   else {
      printf("Computation failed: \n");
      printf("3 arguments expected : the keyword of the calculation, the beta inputfiles and the outputfile. \n");
      printf("For more details see the documentation \n");
   }

fichier_beta=fopen(argv[2],"rt");
if (fichier_beta==NULL)
    {
        puts("Problem to open the input beta file \n!");
        exit(0) ;
    }

fichier_out = fopen(argv[3], "w");
   if (fichier_out==NULL){
	printf("problem to open the output file to write results!\n");
}

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


fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"****           Py_SHS: An open source software      ****\n");
fprintf(fichier_out,"****         about Second Harmonic Scattering       ****\n");
fprintf(fichier_out,"**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n");          
fprintf(fichier_out,"****                  ICGM - ENSCM                  ****\n"); 
fprintf(fichier_out,"****                   september 2023               ****\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"                                                        \n");
fprintf(fichier_out,"           Don't forget to cite this work:              \n");
fprintf(fichier_out,"    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n");
fprintf(fichier_out,"                                                        \n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"*********       SIMULATION RESULTS        **************\n");
fprintf(fichier_out,"********************************************************\n");

if (strcmp(str1,argv[1]) == 0){
printf("Enter the scattered angle in degree (0°=transmission) ? ");
scanf("%lf", &grandthetadegree);
printf("      ******************************     \n");
grandtheta=grandthetadegree*2*pi/360;
printf("********** Begin calculation polar plot at  %lf ° ********\n", grandthetadegree);
//printf("calculation in progress:");

for (i=0;i<nx;i++){
	compteur=i*100/nx;
	printf("calculation in progress: %lf\n",compteur);
	for (j=0;j<ny;j++){
		for (kkk=0;kkk<nz;kkk++){
            		xi = 0 + hx/2 + i*hx;
            		yj = 0 + hy/2 + j*hy;
            		zk = 0 + hz/2 + kkk*hz;
            		
					av=av+hx*hy*hz*betalabxxx_sin(xi, yj, zk,beta);
            		bv=bv+hx*hy*hz*beta_bv(xi, yj, zk,beta);
            		cv=cv+hx*hy*hz*betalabxyy_sin(xi, yj, zk,beta);
            		ah=ah+hx*hy*hz*betalabzxx_sin(xi, yj, zk,beta);
            		bh=bh+hx*hy*hz*beta_bh_90(xi, yj, zk,beta);
            		ch=ch+hx*hy*hz*betalabzyy_sin(xi, yj, zk,beta);
            		
			}
		}
	}
printf("\n"); 

av=av/(8*pi*pi);
bv=bv/(8*pi*pi);
cv=cv/(8*pi*pi);
ah=ah/(8*pi*pi);
bh=bh/(8*pi*pi);
ch=ch/(8*pi*pi);

printf("**********calculation complete ************\n");           
printf("av = %lf\n",av);
printf("bv = %lf\n",bv);
printf("cv = %lf\n",cv);
printf("ah = %lf\n",ah);
printf("bh = %lf\n",bh);
printf("ch = %lf\n",ch);


printf("********************************************\n");   
I4v=(av-bv+cv)/((3*av+bv+3*cv));
I2v=(4)*(av-cv)/(3*av+bv+3*cv);
I4h=(ah-bh+ch)/((3*ah+bh+3*ch));
I2h=(4)*(ah-ch)/(3*ah+bh+3*ch);
Dep=cv/av;
betaHRS=sqrt((av)+(ch));
betaHRS_esu=betaHRS*8.639*pow(10,-3);

printf("I2v = %lf \n",I2v);
printf("I4v = %lf \n",I4v);
printf("I2h = %lf \n",I2h);
printf("I4h = %lf \n",I4h);
printf("DEPOLARIZATION = %lf \n",Dep);
printf("beta_HRS = %lf au\n",betaHRS);
printf("beta_HRS = %lf *10^(-30)esu\n",betaHRS_esu);
//printf("av = %lf  au^2 \n",av);
//printf("cv = %lf  au^2 \n",cv);


fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}

fprintf(fichier_out,"                                      \n");
fprintf(fichier_out," POLAR PLOT - CONFIGURATION AT %lf °  \n", grandthetadegree);
fprintf(fichier_out,"                                    \n");
fprintf(fichier_out,"nombre de dipole = %d  \n",nombredipole);
fprintf(fichier_out,"*********          Hyperpolarizability Tensor      ***********\n");
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
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"I2v=%lf  I4v=%lf \n",I2v,I4v);
fprintf(fichier_out,"I2h=%lf  I4h=%lf \n",I2h,I4h);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"av=%lf   bv=%lf    cv= %lf \n",av,bv,cv);
fprintf(fichier_out,"ah=%lf   bh=%lf    ch= %lf \n",ah,bh,ch);

fclose(fichier_out);

fprintf(fichier_grace1,"gamma IV(gamma) IH(gamma)\n");
for (iii=0;iii<100;iii++){
gammaaa=2*pi*iii/100;
nomc=cos(gammaaa);
noms=sin(gammaaa);
Iv=av*(nomc*nomc*nomc*nomc)+bv*(nomc*nomc*noms*noms)+cv*(noms*noms*noms*noms);
Ih=ah*(nomc*nomc*nomc*nomc)+bh*(nomc*nomc*noms*noms)+ch*(noms*noms*noms*noms);
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
for (i=0;i<nx;i++){
//	compteur=i*100/nx;
//printf("..");
	for (j=0;j<ny;j++){
		for (kkk=0;kkk<nz;kkk++){
            		xi = 0 + hx/2 + i*hx;
            		yj = 0 + hy/2 + j*hy;
            		zk = 0 + hz/2 + kkk*hz;
            		
			av=av+hx*hy*hz*betalabxxx_sin(xi, yj, zk,beta);
            		bv=bv+hx*hy*hz*beta_bv(xi, yj, zk,beta);
            		cv=cv+hx*hy*hz*betalabxyy_sin(xi, yj, zk,beta);
            		ah=ah+hx*hy*hz*betalabyxx_sin(xi, yj, zk,beta);
            		bh=bh+hx*hy*hz*beta_bh_0(xi, yj, zk,beta);
            		ch=ch+hx*hy*hz*betalabyyy_sin(xi, yj, zk,beta);
            		
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
printf("cv = %lf\n",cv);
printf("ah = %lf\n",ah);
printf("bh = %lf\n",bh);
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
fprintf(fichier_out,"* POLAR PLOT - CONFIGURATION AT 0° *\n");
fprintf(fichier_out,"                                    \n");
fprintf(fichier_out,"nombre de dipole = %d  \n",nombredipole);
fprintf(fichier_out,"*********          Hyperpolarizability Tensor      ***********\n");
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
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"I2v=%lf  I4v=%lf \n",I2v,I4v);
fprintf(fichier_out,"I2h=%lf  I4h=%lf \n",I2h,I4h);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"av=%lf   bv=%lf    cv= %lf \n",av,bv,cv);
fprintf(fichier_out,"ah=%lf   bh=%lf    ch= %lf \n",ah,bh,ch);

fclose(fichier_out);

fprintf(fichier_grace1,"gamma IV(gamma) IH(gamma)\n");
for (iii=0;iii<100;iii++){
gammaaa=2*pi*iii/100;
nomc=cos(gammaaa);
noms=sin(gammaaa);
Iv=av*(nomc*nomc*nomc*nomc)+bv*(nomc*nomc*noms*noms)+cv*(noms*noms*noms*noms);
Ih=ah*(nomc*nomc*nomc*nomc)+bh*(nomc*nomc*noms*noms)+ch*(noms*noms*noms*noms);
fprintf(fichier_grace1,"%lf %lf %lf\n",gammaaa,Iv,Ih);
}


fclose(fichier_grace1);

}


else if (strcmp(str3,argv[1]) == 0){


fprintf(fichier_out,"                                      \n");
fprintf(fichier_out,"* SHS angular distribution            \n");
fprintf(fichier_out,"                                      \n");
fprintf(fichier_out,"*********          Hyperpolarizability Tensor      ***********\n");
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
fprintf(fichier_out,  " grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n"); 

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}
fprintf(fichier_grace1," grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n");
printf("Enter the number of points between 0 and Pi for the angular SHS analysis ?"); 
scanf("%d", &pas2);
mmmf=2*(pas2)+1;
mmmfbis=pas2+1;
printf("********** Begin SHS Angular Distribution calculation**************\n");
for (mmm=0;mmm<mmmfbis;mmm++){
	grandtheta=-pi+mmm*(pi/pas2);
        grandthetadegree=(grandtheta*360)/(2*pi);
        printf("angle value: %lf \n",grandthetadegree);
   IV_0b=0;
   IV_90b=0;
   beta_xxx_moy=0;
   beta_xyy_moy=0;
   beta_y_moy_0=0;
   beta_y_moy_90=0;
   beta_z_moy_0=0;
   beta_z_moy_90=0;
	for (i=0;i<nx;i++){
		//compteur=i*100/nx;
	        //printf("calculation in progress: %lf \n",compteur);
		for (j=0;j<ny;j++){
			for (kkk=0;kkk<nz;kkk++){
                xi = 0 + hx/2 + i*hx;
                yj = 0 + hy/2 + j*hy;
                zk = 0 + hz/2 + kkk*hz;
				beta_xxx_moy=beta_xxx_moy+hx*hy*hz*(betalaboxxxint(xi, yj, zk,beta)*(betalaboxxxint(xi, yj, zk,beta)))*sin(yj);
                beta_xyy_moy=beta_xyy_moy+hx*hy*hz*(betalaboxyyint(xi, yj, zk,beta)*(betalaboxyyint(xi, yj, zk,beta)))*sin(yj);
				
				ui1=betalaboyxxint(xi,yj,zk,beta)*betalaboyxxint(xi,yj,zk,beta)*cos(grandtheta)*cos(grandtheta)-2*betalabozxxint(xi,yj,zk,beta)*betalaboyxxint(xi,yj,zk,beta)*cos(grandtheta)*sin(grandtheta)+betalabozxxint(xi,yj,zk,beta)*betalabozxxint(xi,yj,zk,beta)*sin(grandtheta)*sin(grandtheta);
				ui2=betalaboyyyint(xi,yj,zk,beta)*betalaboyyyint(xi,yj,zk,beta)*cos(grandtheta)*cos(grandtheta)-2*betalabozyyint(xi,yj,zk,beta)*betalaboyyyint(xi,yj,zk,beta)*cos(grandtheta)*sin(grandtheta)+betalabozyyint(xi,yj,zk,beta)*betalabozyyint(xi,yj,zk,beta)*sin(grandtheta)*sin(grandtheta);
				
                beta_y_moy_0=beta_y_moy_0+hx*hy*hz*ui1*sin(yj);
                beta_z_moy_0=beta_z_moy_0+hx*hy*hz*ui2*sin(yj);
	
			}
		}
	}
	
	IH_0=beta_y_moy_0;
    IH_90=beta_z_moy_0;
    IV_0=beta_xxx_moy;
    IV_90=beta_xyy_moy;
    printf("IH_0=%lf    IH_90=%lf   IV_0=%lf    IV_90=%lf \n",IH_0,IH_90,IV_0,IV_90);
    fprintf(fichier_out," %lf    %lf   %lf   %lf   %lf \n",grandtheta,IH_0,IH_90,IV_0,IV_90);
    fprintf(fichier_grace1," %lf    %lf   %lf   %lf   %lf \n",grandthetadegree,IH_0,IH_90,IV_0,IV_90);
}
printf("calculation finished \n");
fclose(fichier_out);
fclose(fichier_grace1);


}
   /* IH_0=creal(csqrt(beta_y_moy_0*beta_y_moy_0+beta_z_moy_0*beta_z_moy_0));
    IH_90=creal(csqrt(beta_y_moy_90*beta_y_moy_90+beta_z_moy_90*beta_z_moy_90));
    IV_0=beta_xxx_moy;
    IV_90=beta_xyy_moy;
    printf("IH_0=%lf    IH_90=%lf   IV_0=%lf    IV_90=%lf \n",IH_0,IH_90,IV_0,IV_90);
    fprintf(fichier_out," %lf    %lf   %lf   %lf   %lf \n",grandtheta,IH_0,IH_90,IV_0,IV_90);
    fprintf(fichier_grace1," %lf    %lf   %lf   %lf   %lf \n",grandtheta,IH_0,IH_90,IV_0,IV_90);
*/

else {
printf("--------------------- ----------------------------------------------------\n"); 
printf("             ERROR in the name of the keyword                             \n"); 
printf(" USE ONLY polarplot_single OR polarplot_integrate OR angle_scattering  !! \n"); 
printf(" ! no calculation has been done ! Retry with the correct keyword          \n");
printf("--------------------------------------------------------------------------\n"); 
}






}
