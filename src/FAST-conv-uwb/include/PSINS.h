
/* PSINS(Precise Strapdown Inertial Navigation System) C++ algorithm hearder file PSINS.h

Copyright(c) 2015-2021, by YanGongmin, All rights reserved.
Northwestern Polytechnical University, Xi'an, P.R.China.
Date: 17/02/2015, 19/07/2017, 11/12/2018, 27/12/2019, 12/12/2020, 18/02/2021
*/

#ifndef _PSINS_H
#define _PSINS_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

#pragma pack(4)

/************** compiling control !!! ***************/
#define PSINS_MATRIX_MAX_DIM	34
#define PSINS_IO_FILE
#define PSINS_RMEMORY
//#define PSINS_AHRS_MEMS
//#define PSINS_psinsassert
//#define MAT_COUNT_STATISTIC
//#define PSINS_STACK
//#define PSINS_UART_PUSH_POP
//#define PSINS_VC_AFX_HEADER

// type re-define
#ifndef BOOL
typedef int		BOOL;
#endif

#ifndef BYTE
typedef unsigned char BYTE;
#endif

// constant define
#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif

#ifndef NULL
#define NULL	((void *)0)
#endif

#ifndef PI
#define PI		3.14159265358979
#endif
#define PI_2	(PI/2.0)
#define PI_4	(PI/4.0)
#define _2PI	(2.0*PI)

#define sqrt2	1.414213562373095	// sqrt(2) ...
#define sqrt3	1.732050807568877
#define sqrt5	2.236067977499790
#define sqrt6	2.449489742783178
#define sqrt7	2.645751311064591
#define sqrt8	2.828427124746190

#define DEG		(PI/180.0)		// arcdeg
#define MIN		(DEG/60.0)		// arcmin
#define SEC		(MIN/60.0)		// arcsec
#define HUR		3600.0			// hur
#define SHUR	60.0			// sqrt(hur)
#define DPS		(DEG/1.0)		// deg/s
#define DPH		(DEG/HUR)		// deg/h
#define DPSH	(DEG/SHUR)		// deg/sqrt(h)
#define G0		9.7803267714
#define MG		(G0/1.0e3)		// mg
#define UG		(G0/1.0e6)		// ug
#define UGPSHZ	(UG/1)			// ug/sqrt(Hz)
#define RE		6378137.0		//���򳤰���
#define PPM		1.0e-6

#ifndef EPS
#define EPS		(2.220446049e-16)
#endif
#ifndef INF
#define INF		(3.402823466e+30)
#endif
#define INFp5	(INF*0.5)
#define fEND	(10.0*INF)

#define FRQ1		1				// sampling frequency (FRQ** Hz)
#define FRQ50		50
#define FRQ100		100
#define FRQ125		125
#define FRQ200		200
#define FRQ400		400
#define FRQ500		500
#define FRQ800		800
#define FRQ1000		1000
#define TS1			(1.0/FRQ1000)	// sampling interval (TS** ms)
#define TS1p25		(1.0/FRQ800)
#define TS2			(1.0/FRQ500)
#define TS2p5		(1.0/FRQ400)
#define TS5			(1.0/FRQ200)
#define TS8			(1.0/FRQ125)
#define TS10		(1.0/FRQ100)
#define TS20		(1.0/FRQ50)
#define TS1000		(1.0/FRQ1)

// constant define for short in KF P/Q/R setting
#define fXYZU(X,Y,Z,U)	1.0*(X)*(U),1.0*(Y)*(U),1.0*(Z)*(U)
#define fXYZ(X,Y,Z)		fXYZU(X,Y,Z,1.0)
#define fXXZ(X,Z)		fXYZ(X,X,Z)
#define fXXX(X)			fXYZ(X,X,X)
#define fXX6(X)			fXXX(X),fXXX(X)
#define fXX9(X)			fXX6(X),fXXX(X)
#define fOOO			fXXX(0.0)
#define fOO6			fXX6(0.0)
#define fOO9			fXX9(0.0)
#define fIII			fXXX(1.0)
#define fII6			fXX6(1.0)
#define fII9			fXX9(1.0)
#define fINF3			fXXX(INF)
#define fINF6			fXX6(INF)
#define fINF9			fXX9(INF)
#define fPHI(EN,U)		fXYZU(EN,EN,U,MIN) //�Ƕ�ת���ȣ���λΪ��
#define fLLH(LL,H)		fXXZ((LL)/RE,(H))
#define fPOS(LLH)		fLLH(LLH,LLH)
#define fDEG3(X)		fXXX(X*DEG)
#define fMIN3(X)		fXXX(X*DEG/60.0)
#define fSEC3(X)		fXXX(X*SEC)
#define fDPS3(X)		fXXX(X*DPS)
#define fDPH3(X)		fXXX(X*DPH)
#define fDPSH3(X)		fXXX(X*DPSH)
#define fMG3(X)			fXXX(X*MG)
#define fUG3(X)			fXXX(X*UG)
#define fUGPSHZ3(X)		fXXX(X*UGPSHZ)
#define fdKG1(dkii)			(dkii)*PPM			// dkzz
#define fdKG3(dkii,dkij)	(dkij)*SEC,(dkij)*SEC,(dkii)*PPM		// dkxz,dkyz,dkzz
#define fdKG9(dkii,dkij)	(dkii)*PPM,(dkij)*SEC,(dkij)*SEC,(dkij)*SEC,(dkii)*PPM,(dkij)*SEC,(dkij)*SEC,(dkij)*SEC,(dkii)*PPM
#define fdKA6(dkii,dkij)	(dkii)*PPM,(dkij)*SEC,(dkij)*SEC,(dkii)*PPM,(dkij)*SEC,(dkii)*PPM

#define dbsize(datatype)  ((int)(sizeof(datatype)/sizeof(double)))

#ifdef PSINS_psinsassert
	BOOL	psinsassert(BOOL b);
#else
	#define psinsassert(b)  {};
#endif

#ifndef max
#define max(x,y)        ( (x)>=(y)?(x):(y) )
#endif
#ifndef min
#define min(x,y)        ( (x)<=(y)?(x):(y) )
#endif

#define CC180C360(yaw)  ( (yaw)>0.0 ? (_2PI-(yaw)) : -(yaw) )   // counter-clockwise +-180deg -> clockwise 0~360deg for yaw
#define C360CC180(yaw)  ( (yaw)>=PI ? (_2PI-(yaw)) : -(yaw) )   // clockwise 0~360deg -> counter-clockwise +-180deg for yaw
#define LLH(latitude,longitude,height)	CVect3((latitude)*DEG,(longitude)*DEG,height)
#define PRY(pitch,roll,yaw)				CVect3((pitch)*DEG,(roll)*DEG,(yaw)*PI/180)

#ifdef PSINS_STACK
extern int	psinsstack0, psinsstacksize;
#define stackstart()	{ int stack0=0; psinsstack0=(int)&stack0; }
#define stacksize()		{ int stack1=0; psinsstacksize=max(psinsstacksize,psinsstack0-(int)&stack1); }
#else
#define stackstart()	NULL
#define stacksize()		NULL
#endif

#define disp(i, FRQ, n)  if((i)%((n)*(FRQ))==0) printf("%d\n", (i)/(FRQ))


// class define
class CGLV;
class CVect3;		class CMat3;		class CQuat;		class CVect;		class CMat;
class CEarth;		class CIMU;			class CSINS;		class CAligni0;
class CMaxMin;		class CMaxMinn;		class CSINSGNSS;    class CSINSTDKF;	class CKalman;
		
// function define
double  r2dm(double r);
double	dm2r(double dm);
BOOL	logtrigger(int n, double f0=1.0);
BOOL	IsZero(double f, double eps=EPS);
int		sign(double val, double eps=EPS);
double	range(double val, double minVal, double maxVal);
double	atan2Ex(double y, double x);
double  diffYaw(double yaw, double yaw0);
double	MKQt(double sR, double tau);
double	randn(double mu, double sigma);
#define Swap(a, b, tpe)  { tpe tmp=a; a=b; b=tmp; };
//#define pow2(x)			((x)*(x))
#define asinEx(x)		asin(range(x, -1.0, 1.0))
#define acosEx(x)		acos(range(x, -1.0, 1.0))
#define hit0   (NULL)
#define hit1(t, t10, t11)   ( t10<t && t<=t11 )
#define hit2(t, t10, t11, t20, t21)   ( hit1(t,t10,t11) || hit1(t,t20,t21) )
#define hit3(t, t10, t11, t20, t21, t30, t31)   ( hit2(t,t10,t11,t20,t21) || hit1(t,t30,t31) )
#define hit4(t, t10, t11, t20, t21, t30, t31, t40, t41)   ( hit3(t, t10, t11, t20, t21, t30, t31) || hit1(t,t40,t41) )
CVect3	q2att(const CQuat &qnb);
CMat3	diag(const CVect3 &v);
CMat3	pro_mat(const CVect3 & v);
void	IMURFU(CVect3 *pwm, int nSamples, const char *str);
void	IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str);
vector<string> split(const string & str, const string & delim);

// Matrix Max Dimension define
#define MMD		PSINS_MATRIX_MAX_DIM
#define MMD2	(MMD*MMD)

// global variables
extern const CVect3	O31, One31, Ipos, posNWPU;
extern const CQuat	qI;
extern const CMat3	I33, O33, One33;
extern const CVect  On1, O1n, Onen1;
extern const CGLV	glv;
extern int			psinslasterror;

class CGLV
{
public:
	double X0, Y0;
	double Re, f, g0, wie, L0, UTM;										// the Earth's parameters �����򳤰��ᡢ���ʡ���������ת������
	double e, e2, ep, ep2, Rp;
	double Alpha, Beta, Gamma, Delta;
	double mg, ug, deg, min, sec, hur, ppm, ppmpsh;					// commonly used units
	double dps, dph, dpsh, dphpsh, dph2, dphpg, ugpsh, ugpsHz, ugpg2, mpsh, mpspsh, secpsh;
	//1.0 / 298.257223563
	CGLV(double Re=6378137.0, double f=(1.0/298.257223563), double wie0=7.2921151467e-5, double g0=9.780325,double L0 = 2.042035225);
	//CGLV(double Re = 6378137.0, double f = (1.0 / 298.257223563), double wie0 = 7.2921151467e-5, double g0 = 9.7803267714);
};
class CVect3 
{
public:
	double i, j, k;

	CVect3(void);
	CVect3(double xyz);
	CVect3(double xx, double yy, double zz);
	CVect3(const double *pdata);
	CVect3(const float *pdata);

	CVect3& operator=(double f);							// every element equal to a same double
	CVect3& operator=(const double *pf);					// vector equal to a array
	friend BOOL IsZero(const CVect3 &v, double eps);	// psinsassert if all elements are zeros
	friend BOOL IsZeroXY(const CVect3 &v, double eps);	// psinsassert if x&&y-elements are zeros
	friend BOOL IsNaN(const CVect3 &v);						// psinsassert if any element is NaN
	CVect3 operator+(const CVect3 &v) const;				// vector addition
	CVect3 operator-(const CVect3 &v) const;				// vector subtraction
	CVect3 operator*(const CVect3 &v) const;				// vector cross multiplication
	CVect3 operator*(const CMat3 &m) const;					// row-vector multiply matrix
	CVect3 operator*(double f) const;						// vector multiply scale
	CVect3 operator/(double f) const;						// vector divide scale
	CVect3 operator/(const CVect3 &v) const;				// vector divide vect3 element by element
	CVect3& operator+=(const CVect3 &v);					// vector addition
	CVect3& operator-=(const CVect3 &v);					// vector subtraction
	CVect3& operator*=(double f);							// vector multiply scale
	CVect3& operator/=(double f);							// vector divide scale
	CVect3& operator/=(const CVect3 &v);					// vector divide vect3 element by element
	double& operator()(int r);					// vector element
	friend CVect3 operator*(double f, const CVect3 &v);		// scale multiply vector
	friend CVect3 operator-(const CVect3 &v);				// minus
	friend CMat3 vxv(const CVect3 &v1, const CVect3 &v2);	// column-vector multiply row-vector, v1*v2'
	friend CVect3 abs(const CVect3 &v);						// abs for each element
	friend double norm(const CVect3 &v);					// vector norm
	friend double normInf(const CVect3 &v);					// vector inf-norm
	friend double normXY(const CVect3 &v);					// vector norm of X & Y components
	friend CVect3 sqrt(const CVect3 &v);					// sqrt
	friend CVect3 pow(const CVect3 &v, int k);			// power
	friend double pow_2(const CVect3 v);
	friend double dot(const CVect3 &v1, const CVect3 &v2);	// vector dot multiplication
	friend CVect3 dotmul(const CVect3 &v1, const CVect3 &v2);	// vector dot multiplication '.*'
	friend CMat3 a2mat(const CVect3 &att);					// Euler angles to DCM 
	friend CVect3 m2att(const CMat3 &Cnb);					// DCM to Euler angles 
	friend CQuat a2qua(double pitch, double roll, double yaw);	// Euler angles to quaternion
	friend CQuat a2qua(const CVect3 &att);					// Euler angles to quaternion
	friend CVect3 q2att(const CQuat &qnb);					// quaternion to Euler angles 
	friend CQuat rv2q(const CVect3 &rv);					// rotation vector to quaternion
	friend CVect3 q2rv(const CQuat &q);						// quaternion to rotation vector
	friend CMat3 askew(const CVect3 &v);					// askew matrix;
	friend double sinAng(const CVect3 &v1, const CVect3 &v2); // |sin(angle(v1,v2))|
	friend CMat3 pos2Cen(const CVect3 &pos);				// to geographical position matrix
	friend CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts, CEarth *pEth);  // position difference to velocity
	friend CVect3 MKQt(const CVect3 &sR, const CVect3 &tau);// first order Markov white-noise variance calculation
	friend CMat3 dv2att(const CVect3 &va1, const CVect3 &va2, const CVect3 &vb1, const CVect3 &vb2);  // attitude determination using double-vector
	friend CVect3 vn2att(const CVect3 &vn);  // trans ENU velocity to attitude (pitch & yaw)
	friend CVect3 Alignsb(const CVect3 wmm, const CVect3 vmm, double latitude);  // align in static-base
	friend double MagYaw(const CVect3 &mag, const CVect3 &att, double declination);
	friend CVect3 xyz2blh(const CVect3 &xyz);				// ECEF X/Y/Z to latitude/longitude/height
	friend CVect3 blh2xyz(const CVect3 &blh);				// latitude/longitude/height to ECEF X/Y/Z 
	friend CVect3 xyh2blh(const CVect3& xyh);
	friend CVect3 blh2xyh(const CVect3& blh);
	friend CVect3 Vxyz2enu(const CVect3 &Vxyz, const CVect3 &pos);  // ECEF Vx/Vy/Vz to Ve/Vn/Vu
	friend CVect3 randn(const CVect3 &mu, const CVect3 &sigma);
	friend CVect3 sort(const CVect3 &v);
};

class CQuat
{
public:
	double q0, q1, q2, q3;

	CQuat(void);
	CQuat(double qq0, double qq1=0.0, double qq2=0.0, double qq3=0.0);
	CQuat(const double *pdata);

	CQuat operator+(const CVect3 &phi) const;	// true quaternion add misalign angles
	CQuat operator-(const CVect3 &phi) const;	// calculated quaternion delete misalign angles
	CVect3 operator-(CQuat &quat) const;		// get misalign angles from calculated quaternion & true quaternion
	CQuat operator*(const CQuat &q) const;		// quaternion multiplication
	CVect3 operator*(const CVect3 &v) const;	// quaternion multiply vector
	CQuat& operator*=(const CQuat &q);			// quaternion multiplication
	CQuat& operator-=(const CVect3 &phi);		// calculated quaternion delete misalign angles
	void SetYaw(double yaw=0.0);				// set Euler angles to designated yaw
	friend void normlize(CQuat *q);				// quaternion norm
	friend CQuat operator~(const CQuat &q);		// quaternion conjugate
	friend CVect3 qq2phi(const CQuat &qcalcu, const CQuat &qreal);
	friend CQuat UpDown(const CQuat &q);		// Up-Down the quaternion represented attitide
};

class CMat3 
{
public:
	double e00, e01, e02, e10, e11, e12, e20, e21, e22;

	CMat3(void);
	CMat3(double xyz);
	CMat3(double xx, double yy, double zz);
	CMat3(double xx, double xy, double xz,
		  double yx, double yy, double yz,
		  double zx, double zy, double zz );
	CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2, BOOL isrow=1);  // M = [v0; v1; v2]

	CMat3 operator+(const CMat3 &m) const;					// matrix addition
	CMat3 operator-(const CMat3 &m) const;					// matrix subtraction
	CMat3 operator*(const CMat3 &m) const;					// matrix multiplication
	CMat3 operator*(double f) const;						// matrix multiply scale
	CVect3 operator*(const CVect3 &v) const;				// matrix multiply vector
	CMat3& operator+=(const CMat3 &m);						// matrix +=
	CMat3 operator+(const CVect3 &v) const;					// matrix addition
	CMat3& operator+=(const CVect3 &v);						// matrix + diag(vector)
	void SetRow(int i, CVect3 &v);							// set i-row from vector
	void SetClm(int i, CVect3 &v);							// set i-column from vector
	CVect3 GetRow(int i) const;								// get i-row from matrix
	CVect3 GetClm(int i) const;								// get i-column from matrix
	friend CMat3 rcijk(const CMat3 &m, int ijk);			// re-arrange row/clm indexed by ijk
	friend CMat3 operator-(const CMat3 &m);					// minus
	friend CMat3 operator~(const CMat3 &m);					// matrix transposition
	friend CMat3 operator*(double f, const CMat3 &m);		// scale multiply matrix
	friend void symmetry(CMat3 &m);							// matrix symmetrization
	friend CMat3 pow(const CMat3 &m, int k);				// k^th power
	friend double trace(const CMat3 &m);					// matrix trace
	friend double det(const CMat3 &m);						// matrix determinat
	friend CMat3 adj(const CMat3 &m);						// 3x3 adjoint matrix
	friend CMat3 inv(const CMat3 &m);						// 3x3 matrix inverse
	friend CVect3 diag(const CMat3 &m);						// the diagonal of a matrix
	friend CMat3 diag(const CVect3 &v);						// diagonal matrix
	friend CMat3 MMT(const CMat3 &m1, const CMat3 &m2);	// m=m1*m2^T
	friend double trMMT(const CMat3 &m1, const CMat3 &m2);	// trace(m1*m2^T)
	friend double norm(const CMat3 &m);						// matrix norm
	friend CQuat m2qua(const CMat3 &Cnb);					// DCM to quaternion
	friend CMat3 q2mat(const CQuat &qnb);					// attitude quaternion to DCM
	friend CMat3 sfoam(const CMat3 &B, int iter);		// Supper Fast Optimal Attitude Matrix(SFOAM)
	friend CMat3 randn(const CMat3 &mu, const CMat3 &sigma);
	friend CMat3 pro_mat(const CVect3 &v);
};

class CVect
{
public:
	int row, clm, rc;//�����С����Լ�Ԫ������
	double dd[MMD] = {0};//�洢����Ԫ�ص�����

	CVect(void);
	CVect(int row0, int clm0=1);
	CVect(int row0, double f);
	CVect(int row0, double f, double f1, ...);
	CVect(int row0, const double *pf);
	CVect(const CVect3 &v);
	CVect(const CVect3 &v1, const CVect3 v2);

	void Set(double f, ...);
	void Set2(double f, ...);
	CVect operator+(const CVect &v) const;		// vector addition
	CVect operator-(const CVect &v) const;		// vector subtraction
	CVect operator*(double f) const;			// vector multiply scale
	CVect& operator=(double f);					// every element equal to a same double
	CVect& operator=(const double *pf);			// vector equal to a array
	CVect& operator=(const CMat3 &m);			// vector equal to matrix 3x3
	CVect& operator+=(const CVect &v);			// vector addition
	CVect& operator-=(const CVect &v);			// vector subtraction
	CVect& operator*=(double f);				// vector multiply scale
	CVect operator*(const CMat &m) const;		// row-vector multiply matrix
	CMat operator*(const CVect &v) const;		// 1xn vector multiply nx1 vector, or nx1 vector multiply 1xn vector
	double& operator()(int r);					// vector element
	friend double dot(const CVect &v1, const CVect &v2);	// vector dot multiplication
	friend CVect dotmul(const CVect &v1, const CVect &v2);	// vector dot multiplication '.*'
	friend CVect operator~(const CVect &v);		// vector transposition
	friend CVect abs(const CVect &v);			// vector abs for each element
	friend double norm(const CVect &v);			// vector norm
	friend double normInf(const CVect &v);		// inf-norm
	friend CVect pow(const CVect &v, int k);	// vector element power
	friend CVect randn(const CVect &mu, const CVect &sigma);
	friend CVect sort(const CVect &v);
};

class CMat
{
public:
	int row, clm, rc;
	double dd[MMD2] = {0};

	CMat(void);
	CMat(int row0, int clm0);
	CMat(int row0, int clm0, double f);
	CMat(int row0, int clm0, double f, double f1, ...);
	CMat(int row0, int clm0, const double *pf);

	void Clear(void);
	void SetDiag(double f, ...);
	void SetDiag2(double f, ...);
	void IsSquareMat(void);
	void ChangeRow(int r1,int r2);
	void Elim(void);//����Ԫ��ȥ��
	double Arg(void);
	bool IsPtv(void);
	CMat Inv(void);
	CMat operator+(const CMat &m) const;				// matrix addition
	CMat operator-(const CMat &m) const;				// matrix subtraction
	CMat operator*(double f) const;						// matrix multiply scale
	CVect operator*(const CVect &v) const;				// matrix multiply vector
	CMat operator*(const CMat &m) const;				// matrix multiplication
	CMat& operator=(double f);							// every element equal to a same double
	CMat& operator+=(const CMat &m0);					// matrix addition
	CMat& operator+=(const CVect &v);					// matrix + diag(vector)
	CMat& operator-=(const CMat &m0);					// matrix subtraction
	CMat& operator*=(double f);							// matrix multiply scale
	CMat& operator++();									// 1.0 + diagonal
	double& operator()(int r, int c=-1);				// get element m(r,c)
	void ZeroRow(int i);								// set i-row to 0
	void ZeroClm(int j);								// set j-column to 0
	void SetRow(int i, double f, ...);					// set i-row from n-double
	void SetRow(int i, const CVect &v);					// set i-row from vector
	void SetClm(int j, double f, ...);					// set j-column from n-double
	void SetClm(int j, const CVect &v);					// set j-column from vector
	CVect GetRow(int i) const;							// get i-row from matrix
	CVect GetClm(int j) const;							// get j-column from matrix
	void SetRowVect3(int i, int j, const CVect3 &v);	// set i-row&j...(j+2)-column from CVect3
	void SetClmVect3(int i, int j, const CVect3 &v);	// set i...(i+2)-row&j-column from CVect3
	void SetDiagVect3(int i, int j, const CVect3 &v);	// m(i,j)=v.i, m(i+1,j+1)=v.j, m(i+2,j+2)=v.k;
	CVect3 GetDiagVect3(int i, int j=-1);				// return CVect3(m(i,j), m(i+1,j+1), m(i+2,j+2))
	void SetAskew(int i, int j, const CVect3 &v);		// set i...(i+2)-row&j...(j+2)-comumn from askew CVect3
	void SetMat3(int i, int j, const CMat3 &m);			// set i...(i+2)-row&j...(j+2)-comumn from CMat3
	CMat3 GetMat3(int i, int j=-1) const;				// get CMat3 from i...(i+2)-row&j...(j+2)-comumn
	void SubAddMat3(int i, int j, const CMat3 &m);		// add i...(i+2)-row&j...(j+2)-comumn with CMat3 m
	friend CMat operator~(const CMat &m);				// matrix transposition
	friend void symmetry(CMat &m);						// matrix symmetrization
	friend double normInf(CMat &m);						// inf-norm
	friend CMat dotmul(const CMat &m1, const CMat &m2);	// matrix dot multiplication '.*'
	friend CVect diag(const CMat &m);					// diagonal of a matrix
	friend CMat diag(const CVect &v);					// diagonal matrix
	friend void RowMul(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1
	friend void RowMulT(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1'
	friend void DVMDVafa(const CVect &V, CMat &M, double afa);	// M = diag(V)*M*diag(V)*afa
	
#ifdef MAT_COUNT_STATISTIC
	static int iCount, iMax;
	~CMat(void);
#endif
};

class CEarth
{
public:
	double a, b;
	double f, e, e2;//e�������һƫ����
	double wie;//������ת���ٶ�

	double sl, sl2, sl4; //γ�����ҡ�ƽ���Լ��Ĵη�
	double cl; //γ������
	double tl; //γ������

	double RM, RN;//RM������Ȧ���ʰ뾶��RN��î��Ȧ���ʰ뾶�� RM=a(1-e^2)(1-e^2*sin^2 B)^(-3/2); RN=a/((1-e^2*sin^2 B)^(1/2))
	double RMh, RNh;//���ʰ뾶�Ӵ�ظ�h
//	double clRNh, f_RMh, f_RNh, f_clRNh;
	CVect3 pos, vn, wnie, wnen, wnin, gn, gcc, *pgn;//pos��[γ�ȣ����ȣ���ظ�]
	CVect3 vn_m, wnie_m, wnen_m, gn_m;
	CVect3 vn_p, wnie_p, wnen_p, gn_p;
	CMat3 Mpv, Mpv_m, Mpv_p;
	CEarth(double a0=glv.Re, double f0=glv.f, double g0=glv.g0);
	void Update(const CVect3 &pos, const CVect3 &vn=O31, int isMemsgrade = 0);
	CVect3 vn2dpos(const CVect3 &vn, double ts=1.0) const;
	void Extrap();
};


class CIMU
{
	CMat3 *pgSens, gSens, *pgSens2, gSens2, *pgSensX, gSensX;  // gyro g-sensitivity
	CVect3 *pKa2, Ka2;  // acc quadratic nonlinearity
	char *prfu, rfu[3];
public:
	int nSamples;
	bool preFirst, onePlusPre;
	CVect3 phim, dvbm, wm_1, vm_1;//�������Ч��תʸ������������������ǰһ�ν�/�ٶ���������

	CIMU(void);
	void SetgSens(const CMat3 &gSens0, const CMat3 &gSens20=O33, const CMat3 &gSensX0=O33);
	void SetKa2(const CVect3 &Ka20);
	void SetRFU(const char *rfu0);
	void Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts=0.0);
	friend void IMURFU(CVect3 *pwm, int nSamples, const char *str);
	friend void IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str);
	friend void IMUStatic(CVect3 &wm, CVect3 &vm, CVect3 &att0, CVect3 &pos0, double ts); 
};

class CSINS// sizeof(CSINS)~=3k bytes
{
public:
	double ts, nts, tk, velMax, hgtMin, hgtMax;//���������n��������֮����ܼ����
	CEarth eth;
	CIMU imu;
	CQuat qnb;
	CMat3 Cnb, Cnb0, Cbn, Kg, Ka;
	CVect3 wib, fb, fn, an, web, wnb, att, vn, vb, pos, eb, db, Ka2, tauGyro, tauAcc, _betaGyro, _betaAcc;
	CMat3 Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp;	// for etm
	CVect3 lvr, vnL, posL; CMat3 CW, MpvCnb;		// for lever arm
	CQuat qnbE; CVect3 attE, vnE, posE;				// for extrapolation
	BOOL isOpenloop, isMemsgrade, isNocompasseffect, isOutlever;

	CSINS(const CVect3 &att0, const CVect3 &vn0=O31, const CVect3 &pos0=O31, double tk0=0.0);
	CSINS(const CQuat &qnb0=qI, const CVect3 &vn0=O31, const CVect3 &pos0=O31, double tk0=0.0);
	void Init(const CQuat &qnb0=qI, const CVect3 &vn0=O31, const CVect3 &pos0=O31, double tk0=0.0);    // initialization using quat attitude, velocity & position
	void SetTauGA(const CVect3 &tauG, const CVect3 &tauA);
	void Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts);		// SINS update using Gyro&Acc samples
	void lever(const CVect3 &dL=O31, CVect3 *ppos=NULL, CVect3* pvn=NULL);		// lever arm
	void etm(void);							// SINS error transform matrix coefficients
};

class CAligni0
{
public:
	int velAid, t0, t1, t2;
	CVect3 pos0, vel0, wmm, vmm, vib0, vi0, Pib01, Pib02, Pi01, Pi02, tmpPib0, tmpPi0;
	CQuat qib0b;
	CEarth eth;
	CIMU imu;
	double tk;
	CQuat qnb0, qnb, qnbsb;

	CAligni0(const CVect3 &pos0=O31, const CVect3 &vel0=O31, int velAid=0);
	void Init(const CVect3 &pos0=O31, const CVect3 &vel0=O31, int velAid=0);
	CQuat Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts, const CVect3 &vel=O31);
};

class CKalman
{
public:
	CSINS sins;
	double kftk, zfdafa;//kftk���˲���ʱ
	int nq, nr, measflag, measflaglog, measstop;//״̬ά��������ά���������ʶ�����������ʶ��ֹͣ�������
	CMat Ft,//ϵͳ���� 
		Pk,//Э������
		Hk,//ϵ������ 
		Kk,//�������˲��������
		Fading;
	CVect Xk, //״̬��������19ά����������ʦ�α�p206��0-2 phi ʧ׼�����; 3-5 dvn �ٶ����;
			  //6-8 dpos λ�����; 9-11 eb ����Ư��; 12-14 db �ӱ���ƫ; 15-17 lever �˱����; 18 dt
		Zk, //��������
		Qt, //��������
		Rt,//�������� 
		Vk,//��Ϣ����
		rts, RtTau, measlost, Xmax, //�������ڡ���������ʱ�䡢���ⶪʧʱ�䡢״̬����
		Pmax, Pmin, //������������
		Pset, Zfd, Zfd0,
		Rmax, Rmin, Rbeta, Rb,				// measurement noise R adaptive
		FBTau, FBMax, FBOne, FBOne1, FBXk, FBTotal;	// feedback control ����ʱ�䳣�����������ޡ�����������ǰ����ֵ���ܷ�����
	int Rmaxcount[MMD], Rmaxcount0[MMD];

	CKalman(void);
	CKalman(int nq0, int nr0, const CSINS &sins0);
	void Init(int nq0, int nr0, const CSINS &sins0);				// initialize Qk,Rk,P0...
	void SetRmaxcount(int cnt = 5);
	void SetFt(int nnq=19);			// process matrix setting
	void SetHk(int nnq=19);			// measurement matrix setting
	void SetMeas(const CVect3& pgps, const CVect3& vgps);				// set measurement
	void Update();
	void SetMeasMC(void);	// set measurement
	void Feedback(int nnq, double ts);	// state feedback
	void RtFading(int i, double fdts);			// Rt growing if no measurment
	void TimeUpdate(double kfts, int fback = 1);	// time update
	int MeasUpdate(double fading = 1.0);			// measurement update
	int RAdaptive(int i, double r, double Pr);	// Rt adaptive
	void RPkFading(int i);						// multiple fading
	void SetMeasFlag(int flag);					// measurement flag setting
	void XPConstrain(void);						// Xk & Pk constrain: -Xmax<Xk<Xmax, Pmin<diag(Pk)<Pmax
};

class MyCKF :public CKalman
{
public:
	CMat3 Cne;
	CVect3 Xn, Xr, Vk_1;//nϵλ��������eϵλ����������һʱ���ٶ�
	MyCKF(const CSINS& sins0);
	void Init();
	void CneUpdate();
	void XnUpdate();
	void XrUpdate();

};

//class GNSS_QC
//{
//public:
//	vector<int> sat;
//	vector<int> num, errsat;
//	vector<double> ek;
//	int errnum;
//	CVect X, Z, Vk; 
//	double R, P, Q;
//	CVect3 rk, Xr, Xs;
//	GNSS_QC();
//	bool Update1(MyCKF& ckf, vector<Sat> S);
//	bool Update2(vector<Sat> S, double dtr, CVect3 xyz);
//};


//class CSINSTDKF:public CKalman
//{
//public:
//	double meantdts, tdts, Pz0, innovation;//tdts���ܲ���ʱ�䣿
//	int iter, ifn, adptOKi, measRes, tdStep, maxStep, curOutStep, maxOutStep;
//	CMat Fk, Pk1;
//	CVect Pxz, Qk, Kk, Hi, tmeas;
//	CVect3 meanfn;
//	CSINS sins;
//
//	CSINSTDKF(void);
//	CSINSTDKF(int nq0, int nr0);
//	void Init(const CSINS& sins0);
//	void TDReset(void);
//	int TDUpdate(const CVect3* pwm, const CVect3* pvm, int nSamples, double ts, int nStep = 1);  // Time-Distributed Update
//	void MeasUpdate(const CVect& Hi, double Ri, double Zi);
//	void MarkovGyro(const CVect3& tauG, const CVect3& sRG, int stateeb = 9);
//	void MarkovAcc(const CVect3& tauA, const CVect3& sRA, int statedb = 12);
//	void SetYaw(double yaw, int statephi = 0, int statedvn = 3);
//	virtual void RTOutput(void) {};  // real-time output before KF update i.e. just after SINS update
//	virtual void Miscellanous(void) {};
//	virtual void SecretAttitude(void) {};
//};
//
//class CSINSGNSS :public CSINSTDKF	// sizeof(CSINSGNSS)~=21k bytes for 15-state KF
//{
//public:
//	double posGNSSdelay, vnGNSSdelay, yawGNSSdelay, dtGNSSdelay, kfts;
//	CVect3 lvGNSS;
//	double tmeas_gps;
//	CSINSGNSS(void);
//	CSINSGNSS(int nq0, int nr0, double ts);
//	void Init(const CSINS& sins0, int grade = -1);
//	virtual void SetFt(int nnq);
//	virtual void SetHk(int nnq);
//	virtual void Feedback(int nnq, double fbts);
//	virtual void SetMeas(void) {};
//	void SetMeasGNSS(const CVect3& pgnss = O31, const CVect3& vgnss = O31);
//
//	void SetMeasMC(void);
//	int Update(const CVect3* pwm, const CVect3* pvm, int nSamples, double ts, int nSteps = 5);
//};

class CBytemani {
public:
	friend unsigned short swap16(unsigned short ui16);
	friend unsigned int swap32(unsigned int ui32);
	friend unsigned long swap64(unsigned long ui64);
	friend unsigned char* int24(unsigned char *pchar3, int int32);
	friend unsigned char* swap24(unsigned char* puc3, unsigned char* pres);
	friend int diffint24(const unsigned char *pc1, const unsigned char *pc0);
};

class imuDataReader
{
public:
	imuDataReader();
	~imuDataReader();
	bool getDataFromFile(string filename, int n);
	vector<CVect3> pwm;
	vector<CVect3> pvm;
	vector<CVect3> vgps;
	vector<CVect3> pgps_xyz;
	vector<double> gpstime;
	vector<CVect3> pgps;
	vector<double> timu;
	vector<double> tgps;
	vector<int> sign;
};
typedef struct {
	unsigned short head, chksum; 
	float t, gx, gy, gz, ax, ay, az, magx, magy, magz, bar, pch, rll, yaw, ve, vn, vu,
		lon0, lon1, lat0, lat1, hgt, gpsve, gpsvn, gpsvu, gpslon0, gpslon1, gpslat0, gpslat1, gpshgt,
		gpsstat, gpsdly, tmp, rsv;
} PSINSBoard;

#pragma pack()

#endif // _PSINS_H
