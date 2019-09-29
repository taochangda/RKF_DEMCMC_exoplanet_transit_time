/*
This is the function for searching transit time using RKF78
Version 1.0
Tao Changda, June, 2019, All rights reserved.
*/

/*
two statements that need special attention:
1,
int k_v funciotn:using different formula for different number of planets,for exmaple:
-----begin----for two planets----------
temp_k[i][0] = -(mc + mp[0])*position_array[i][0] / dst[0] - mp[1] * ((position_array[i][0] - position_array[i][1]) / diff[0][1] + position_array[i][1] / dst[1]);
temp_k[i][1] = -(mc + mp[1])*position_array[i][1] / dst[1] - mp[0] * ((position_array[i][1] - position_array[i][0]) / diff[1][0] + position_array[i][0] / dst[0]);
-----end----for two planets----------

2,
IN transit_time function:
--begin--notice:it will be changed different sutiation such as: different fitting parameter,differnet number of planets,and so on
    p, e, w, cm, mp, ci, omega --
    for (i = 0; i < N_planet; i++)
    {
        p[i] = p_para[i]; 
        e[i] = p_para[i + N_planet];
        w[i] = p_para[i + 2 * N_planet];
        cm0[i] = p_para[i + 3 * N_planet];
        mp[i] = p_para[i + 4 * N_planet]; 
        
        ci[i] = pi*0.5;
        omega[i] =0.0;

--end--notice:it will be changed different sutiation such as: different fitting parameter,differnet number of planets,and so on
p, e, w, cm, mp, ci, omega --
*/




#include"mex.h" /*for MATLAB API*/
#include"math.h"
/*#include "stdio.h" */
#include"string.h"


#define N_transit 500 /*Number of transit time in tstop for all planets, set a litte big in facts ����ʱ��tstop�������������ǵ������ܵ����ֵ��ʵ�����ȡ��һ��*/
#define N_planet  7  /*Number of planet ��������*/
#define N_para 35  /*Number of fitting parameter ��ϵĲ�������*/
#define N_order 13  /*RKF78:13��RKF45:6*/

const double a[N_order+1][N_order-1] = { {0,0,0,0,0,0,0,0,0,0,0,0},{0.0740740740740741,0,0,0,0,0,0,0,0,0,0,0},{0.0277777777777778,0.0833333333333333,0,0,0,0,0,0,0,0,0,0},{0.0416666666666667,0,0.125000000000000,0,0,0,0,0,0,0,0,0},{0.416666666666667,0,-1.56250000000000,1.56250000000000,0,0,0,0,0,0,0,0},{0.0500000000000000,0,0,0.250000000000000,0.200000000000000,0,0,0,0,0,0,0},{-0.231481481481481,0,0,1.15740740740741,-2.40740740740741,2.31481481481482,0,0,0,0,0,0},{0.103333333333333,0,0,0,0.271111111111111,-0.222222222222222,0.0144444444444444,0,0,0,0,0},{2,0,0,-8.83333333333333,15.6444444444444,-11.8888888888889,0.744444444444445,3,0,0,0,0},{-0.842592592592593,0,0,0.212962962962963,-7.22962962962963,5.75925925925926,-0.316666666666667,2.83333333333333,-0.0833333333333333,0,0,0},{0.581219512195122,0,0,-2.07926829268293,4.38634146341463,-3.67073170731707,0.520243902439024,0.548780487804878,0.274390243902439,0.439024390243902,0,0},{0.0146341463414634,0,0,0,0,-0.146341463414634,-0.0146341463414634,-0.0731707317073171,0.0731707317073171,0.146341463414634,0,0},{-0.433414634146341,0,0,-2.07926829268293,4.38634146341463,-3.52439024390244,0.534878048780488,0.621951219512195,0.201219512195122,0.292682926829268,0,1},{0,0,0,0,0,0,0,0,0,0,0,0} };
/*a: RKF78 constant parameter, an all zeros line is added for processing expediently.  ������Ϊ�̶��������Ҷ����һ��ȫ0ֵ��Ϊ���ǳ������ķ���*/
const double pi = 3.14159265358979323846264338328; /*cycle rate Բ����*/

void k_v(double position_array[3][N_planet], double mc, double mp[N_planet], double temp_k[3][N_planet]);/*compute k-value of velocity function.���ٶȷ���(vx,vy,vz)Kֵ�ĺ���*/
void rkf78(double *h, double *t, double position_array[3][N_planet], double velocity_array[3][N_planet], double mc, double mp[N_planet],double err);/*RKF78 integration function.RKF78���ֺ���*/
void transit_time(double *tr, double *error_flag, double p_para[N_para], double mc, double Rs, double tstop, double time_precise,double err);/*search transit time funcion Ѱ������ʱ�̺���*/
/*--notice��for pointer variable, it will be changed if it is changed in the CALL fuction----*/
/*--ע�⣺�����ж��ڴ��ݵĲ�����ָ��(����)���͵ı������ں������и���ֵ������˺����壬��ЩֵҲ�Ǹ��º��----*/

/*
 * int main()
 * {...
 * return 0;
 * }
 */


/*-------------matlab API function bellow------------------------*/
/*---refer to matlab offical guide for more imformation of mexFunction----*/
/*---���ڽӿں��� mexFunction����ϸ�÷���ο�matlab�Ĺٷ��ĵ�----*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*input*/
    double *p_para; /*fitting parameter ��Ҫ��ϵ�N_para������*/
    double *mc; /*mass of star ��������*/
    double *Rs; /*radis of star ���ǰ뾶*/
    double *tstop; /*integration time ������ʱ��*/
    double *time_precise; /*minimal time interval for searching transit time Ѱ������ʱ�̷����ж�ʱ�䲽���ľ��ȿ��ƣ���������ƣ����ܻ����ʱ�䲽����С�������������࣬���Ĺ���ļ���ʱ��*/
    double *err; /*precision control parameter RKF���ȿ��Ʋ���*/
    
    /*outuput*/
    double *error_flag; /* 0: reasonable, 1: unreasonable ��־λ����time_precise���ʹ�ã�����ֵ0��ʾ������������ֵ1��ʾ������*/
    double *tr; /*transit time array, notices: elements is stored orderly base on row in C languge but on column in matlab ����ʱ�����飬C�����ǰ���˳��洢��MATLAB�ǰ���˳��洢����Ҫע��*/

    
    if (nrhs!=6)
        mexErrMsgTxt("Six inputs required.");
    if (nlhs != 2)
        mexErrMsgTxt("Two output required.");
    
    p_para = mxGetPr(prhs[0]);
    mc = mxGetPr(prhs[1]);
    Rs = mxGetPr(prhs[2]);
    tstop = mxGetPr(prhs[3]);
    time_precise = mxGetPr(prhs[4]);
    err = mxGetPr(prhs[5]);
        
    plhs[0] = mxCreateDoubleMatrix(N_transit, N_planet, mxREAL); /*notices: elements is stored orderly base on row in C languge but on column in matlab.C�����ǰ���˳��洢��MATLAB�ǰ���˳��洢����Ҫע��*/
    tr = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    error_flag = mxGetPr(plhs[1]);
    
    transit_time(tr, error_flag, p_para, *mc, *Rs, *tstop, *time_precise,*err); /*CALL function transit_time */
    
}
/*-------------matlab API function above------------------------*/



void transit_time(double *tr, double *error_flag, double p_para[N_para], double mc, double Rs, double tstop, double time_precise,double err)
{
    double t = 0.0; /*initail integration time ��ʼ����ʱ��*/
    double h = 0.02 * pi; /*initail integration time interval ��ʼ����ʱ����*/
    double t1 = 0.0;
    double h1 = 0.0;
    double p[N_planet], e[N_planet], ci[N_planet], w[N_planet], cm0[N_planet], mp[N_planet], omega[N_planet]; /*�������ǹ����������*/
    double ap[N_planet], cn[N_planet], ce[N_planet], d[N_planet], dd[N_planet]; /*�м����*/
    double position_array[3][N_planet], velocity_array[3][N_planet], temp_poistion_1[3][N_planet], temp_velocity_1[3][N_planet]; /*����λ�ú��ٶ�����*/
    double px, py, pz, qx, qy, qz, r; /*temp variable �м����*/
    double ce0, gap; /*temp variable �м����*/
    int num[N_planet]; /*array for record number of transit times ��¼����ʱ������������*/
    int flag[N_planet]; /**/
    int i, j,loop_index; /*loop index variable ѭ������*/
    
    /*
     * for (j = 0; j < N_planet * 3; j++)
     * {
     * temp_poistion_1[j % 3][j / 3] =0.0;
     * temp_velocity_1[j % 3][j / 3] = 0.0;
     * position_array[j % 3][j / 3] =0.0;
     * velocity_array[j % 3][j / 3] = 0.0;
     * }
     */
    
    /*-----initialization----------*/
    memset(position_array, 0, 3 * N_planet);
    memset(velocity_array, 0, 3 * N_planet);
    memset(temp_poistion_1, 0, 3 * N_planet);
    memset(temp_velocity_1, 0, 3 * N_planet);
    
    for (i = 0; i < N_planet; i++)
    {
        num[i] = 0;
        flag[i] = 0;
    }
     /*-----initialization----------*/
    
    
    
    /*--begin--notice:it will be changed different sutiation such as: different fitting parameter,differnet number of planets,and so on
    Ҫ�ģ�-��ͬ�����(������������һ������ϲ�����һ��)-- ����˳����ܲ�һ����p, e, w, cm, mp, ci, omega --*/
    for (i = 0; i < N_planet; i++)
    {
        p[i] = p_para[i]; /*period �������*/
        e[i] = p_para[i + N_planet];
        w[i] = p_para[i + 2 * N_planet];
        cm0[i] = p_para[i + 3 * N_planet];
        mp[i] = p_para[i + 4 * N_planet]; /*mass ����*/
        
        ci[i] = pi*0.5;
        omega[i] =0.0;

    /*--end--notice:it will be changed different sutiation such as: different fitting parameter,differnet number of planets,and so on
    Ҫ�ģ�-��ͬ�����(������������һ������ϲ�����һ��)-- ����˳����ܲ�һ����p, e, w, cm, mp, ci, omega --*/
        
        
          /*--begin---caculate initial position(x, y, z) and velocity(vx, vy, vz), base on parameter----------*/
        ap[i] = pow((mc + mp[i]), 1.0 / 3)*pow((p[i] / 365.2422), 2.0 / 3);
        cn[i] = sqrt((mc + mp[i]) / (ap[i] * ap[i] * ap[i])); /*cn[i] = sqrt((mc + mp[i]) / pow(ap[i], 3));*/
        ce[i] = cm0[i];
        do
        {
            ce0 = ce[i];
            ce[i] = ce0 - (ce0 - e[i] * sin(ce0) - cm0[i]) / (1.0 - e[i] * cos(ce0));
        } while (fabs(ce0 - ce[i]) > err);
        
        
        px = cos(omega[i])*cos(w[i]) - sin(omega[i])*sin(w[i])*cos(ci[i]);
        py = sin(omega[i])*cos(w[i]) + cos(omega[i])*sin(w[i])*cos(ci[i]);
        pz = sin(w[i])*sin(ci[i]);
        qx = -cos(omega[i])*sin(w[i]) - sin(omega[i])*cos(w[i])*cos(ci[i]);
        qy = -sin(omega[i])*sin(w[i]) + cos(omega[i])*cos(w[i])*cos(ci[i]);
        qz = cos(w[i])*sin(ci[i]);
        
        position_array[0][i] = ap[i] * (cos(ce[i]) - e[i])*px + ap[i] * sqrt(1 - e[i] * e[i])*sin(ce[i])*qx; /*x*/
        position_array[1][i] = ap[i] * (cos(ce[i]) - e[i])*py + ap[i] * sqrt(1 - e[i] * e[i])*sin(ce[i])*qy; /*y*/
        position_array[2][i] = ap[i] * (cos(ce[i]) - e[i])*pz + ap[i] * sqrt(1 - e[i] * e[i])*sin(ce[i])*qz; /*z*/
        r = sqrt(position_array[0][i] * position_array[0][i] + position_array[1][i] * position_array[1][i] + position_array[2][i] * position_array[2][i]);
        
        velocity_array[0][i] = -ap[i] * ap[i] * cn[i] / r * (sin(ce[i])*px - sqrt(1 - e[i] * e[i])*cos(ce[i])*qx); /*vx*/
        velocity_array[1][i] = -ap[i] * ap[i] * cn[i] / r * (sin(ce[i])*py - sqrt(1 - e[i] * e[i])*cos(ce[i])*qy); /*vy*/
        velocity_array[2][i] = -ap[i] * ap[i] * cn[i] / r * (sin(ce[i])*pz - sqrt(1 - e[i] * e[i])*cos(ce[i])*qz); /*vz*/
    }
       /*-end---caculate initial position(x, y, z) and velocity(vx, vy, vz), base on parameter----------*/
    
    *error_flag=0; 
    
    /*-begin---search transit time for all planets in "tstop" time ----------*/
    while (t < tstop)
    {
        if (h<time_precise) /*if time interval h <time_precise, stop searching, and this parameters(p_para) is abandoned,*/
        {                             /*h��С������ķѷǳ���ļ���ʱ�䡣�������������������Ҳ�ȷ���Ƿ����*/
         *error_flag=1;
         break;
        }
        
       /*Central Point of star is orignal Three-dimensional coordinates, if d=sqrt(x^2+y^2) is minimal and z>0,transit time is find
       �Ժ�������Ϊ����ԭ�㣬���߷���Ϊ z,������ά����ϵ(xyz)����d=sqrt(x^2+y^2)��С��z>0ʱ������Ϊ�ҵ�������ʱ��*/
        for (i = 0; i < N_planet; i++)
        {
            d[i] = sqrt(position_array[0][i] * position_array[0][i] + position_array[1][i] * position_array[1][i]);/*sqrt(x^2+y^2)*/
            if (d[i] < 5 * Rs&&position_array[2][i] > 0 && flag[i]==0)
            {
                gap = 1.0;
                h1 = h;
                t1 = t;
                
                /*
                 * memcpy(temp_poistion_1, position_array, sizeof(temp_poistion_1));
                 * memcpy(temp_velocity_1, velocity_array, sizeof(temp_velocity_1));
                 */
                                
                for (loop_index = 0; loop_index<3; loop_index++)
                {
                    for (j = 0; j < N_planet; j++)
                    {
                        temp_poistion_1[loop_index][j] = position_array[loop_index][j];
                        temp_velocity_1[loop_index][j] = velocity_array[loop_index][j] ;
                    }
                }
               
                
                dd[i] = d[i];
                while (gap > err&&d[i] > err)
                {
                    if (d[i] > dd[i])
                    {
                        h1 = -0.5*h1;
                    }
                    dd[i] = d[i];
                    rkf78(&h1, &t1, temp_poistion_1, temp_velocity_1, mc, mp,err); /*CALL rkf78 funciton*/
                    d[i] = sqrt(temp_poistion_1[0][i] * temp_poistion_1[0][i] + temp_poistion_1[1][i] * temp_poistion_1[1][i]);/*sqrt(x^2+y^2)*/
                    gap = fabs(dd[i] - d[i]);
                }
                
                
                
                /*-- record transit time to array "tr" �ҵ�����ʱ�̺��¼�������������tr��--*/
                *(tr + i * (N_transit)+num[i]) = (t1 - h1) / 2.0 / pi * 365.2422;      /*--  *(tr + i * (N_transit)+num[i]) = (t1 - h1) *58.130101555758664; --*/
                flag[i] = 1;
                num[i]+=1;
                
            }
            
            if (position_array[2][i] < 0 && flag[i]==1)
            {
                flag[i] = 0;
            }
        }
        /*search next transit time �ҵ�����ʱ�̺������� �������»��֣�����Ѱ��һ������ʱ��*/
        rkf78(&h, &t, position_array, velocity_array, mc, mp,err);  /*CALL rkf78 funciton*/
    }
    /*-end---search transit time for all planets in "tstop" time ----------*/
    
}


void rkf78(double *h, double *t, double position_array[3][N_planet], double velocity_array[3][N_planet], double mc, double mp[N_planet],double err)
{
    int j; /*loop  variable*/
    int loop_index; /*loop  variable*/
    int K_index; /*loop  variable*/
    double k_velocity[3][N_planet]; /*k-value for velocity �ٶ�kֵ����*/
    double K_po[N_order][N_planet * 3]; /*temp variable ����λ�õ���ʱ����������������ʹ�õ�*/
    double K_ve[N_order][N_planet * 3]; /*temp variable �����ٶȵ���ʱ����������������ʹ�õ�*/
    double temp_poistion[3][N_planet]; /*temp variable ����λ�õ���ʱ����*/
    double temp_velocity[3][N_planet]; /*temp variable �����ٶȵ���ʱ����*/
    double derr = 1.0; /*initial error of integration ����������ֵ*/
    double delta_po, delta_ve, max_po, max_ve;/*temp variable ��ʱ����*/
    double err_po[3 * N_planet]; /*error of position integration λ�õĻ������*/
    double err_ve[3 * N_planet]; /*error of position velocity �ٶȵĻ������*/
    
    /*
     * for (j = 0; j < N_planet * 3; j++)
     * {
     * temp_poistion[j % 3][j / 3] =0.0;
     * temp_velocity[j % 3][j / 3] = 0.0;
     * k_velocity[j % 3][j / 3] =0.0;
     * err_po[j]=0.0;
     * err_ve[j]=0.0;
     * for (K_index = 0; K_index < N_order; K_index++)
     * {
     * K_po[K_index][j] = 0.0;
     * K_ve[K_index][j] = 0.0;
     *
     * }
     *
     * }
     */
        
    /*--begin---initialization----*/    
    memset(k_velocity, 0, 3 * N_planet);
    memset(K_po, 0, N_order * 3 * N_planet);
    memset(K_ve, 0, N_order * 3 * N_planet);
    memset(temp_poistion, 0, 3 * N_planet);
    memset(temp_velocity, 0, 3 * N_planet);
    memset(err_po, 0, 3 * N_planet);
    memset(err_ve, 0, 3 * N_planet);
    /*--end---initialization----*/  
        
    
    while (derr > err)
    {
        
        for (loop_index = 0; loop_index<3; loop_index++)
        {
            for (j = 0; j < N_planet; j++)
            {
                temp_poistion[loop_index][j] = position_array[loop_index][j];
                temp_velocity[loop_index][j] = velocity_array[loop_index][j] ;
            }
        }

        for (K_index = 0; K_index < N_order; K_index++)
        {
            for (j = 0; j < N_planet * 3; j++)
            {
                K_po[K_index][j] = temp_velocity[j % 3][j / 3]; /*arrange newly for computing conveniently error �������У��������������*/
            }
            
            k_v(temp_poistion, mc, mp, k_velocity); /*CALL k_v function*/
            
            for (j = 0; j < N_planet * 3; j++)
            {
                K_ve[K_index][j] = k_velocity[j % 3][j / 3];  /*arrange newly for computing conveniently error �������У��������������*/
            }
            
            if (K_index<(N_order-1))   /* intergration using RKF78 ���ݱ�׼RKF78�������л���*/
            {
                for (j = 0; j < N_planet * 3; j++)
                {
                    delta_po = 0.0;
                    delta_ve = 0.0;
                    for (loop_index = 0; loop_index <= K_index; loop_index++)
                    {
                        delta_po += a[K_index + 1][loop_index] * K_po[loop_index][j];
                        delta_ve += a[K_index + 1][loop_index] * K_ve[loop_index][j];
                    }
                    
                    temp_poistion[j % 3][j / 3] = position_array[j % 3][j / 3] + (*h) * delta_po;
                    temp_velocity[j % 3][j / 3] = velocity_array[j % 3][j / 3] + (*h) * delta_ve;
                    
                }
                
            }
            
        }
        
        /*--begin--computing error of integration ���������----*/
        max_po=0;
        max_ve=0;
        for (j = 0; j < N_planet * 3; j++)
        {
            err_po[j] = (*h) *41.0 / 810.0 *fabs(K_po[0][j] + K_po[10][j] - K_po[11][j] - K_po[12][j]);
            err_ve[j] = (*h) *41.0 / 810.0*fabs(K_ve[0][j] + K_ve[10][j] - K_ve[11][j] - K_ve[12][j]);
            if (err_po[j]>max_po)
            {max_po=err_po[j];}
            if (err_ve[j]>max_ve)
            {max_ve=err_ve[j];}
        }
        /*--end--computing error of integration ���������----*/
        
        derr = max_ve;
        if (max_po > max_ve)
        {
            derr = max_po;
        }
        
        
        
        if (derr > err) /* decrease time interval if error is bigger than destination error(err) �������㣬��Сʱ������ٴλ���*/
        {
             (*h) = (*h) *0.5;
        }
        
    }
    
    /*--update position and velocity after integration �����������󣬸���ʱ�䡢λ�ú��ٶȵ�ֵ----*/
    (*t) +=(*h);
    for (j = 0; j < N_planet * 3; j++)
    {
        position_array[j % 3][j / 3] += (*h) * (272 * K_po[5][j] + 216 * (K_po[6][j] + K_po[7][j]) + 27 * (K_po[8][j] + K_po[9][j]) + 41 * (K_po[11][j] + K_po[12][j]))/840.0;
        velocity_array[j % 3][j / 3] += (*h) * (272 * K_ve[5][j] + 216 * (K_ve[6][j] + K_ve[7][j]) + 27 * (K_ve[8][j] + K_ve[9][j]) + 41 * (K_ve[11][j] + K_ve[12][j]))/840.0;
        
    }
    
    
}

void k_v(double position_array[3][N_planet], double mc, double mp[N_planet], double temp_k[3][N_planet])
{
    int i;
    int j;
    double dst[N_planet]; /*distance from original point ��ԭ�����*/
    double diff[N_planet][N_planet]; /*distance of two position ���������֮��ľ���*/
    double temp_var;
    
    /*
     * for (i = 0; i < N_planet; i++)
     * {
     * dst[i]=0.0;
     * for (j = 0; j < N_planet; j++)
     * {
     * diff[i][j]=0.0;
     * }
     * }
     */
    
    /*-initialization-*/
    memset(dst, 0, N_planet);
    memset(diff, 0, N_planet*N_planet);
    
    
    
    for (i = 0; i < N_planet; i++)
    {
        temp_var = (position_array[0][i] * position_array[0][i] + position_array[1][i] * position_array[1][i] + position_array[2][i] * position_array[2][i]);
        dst[i] = sqrt(temp_var*temp_var*temp_var); /*sqrt(x^2+y^2+z^2)*/
        
        for (j = 0; j < N_planet; j++)
        {
            
            temp_var = (position_array[0][i] - position_array[0][j])*(position_array[0][i] - position_array[0][j]) + (position_array[1][i] - position_array[1][j])*(position_array[1][i] - position_array[1][j]) + (position_array[2][i] - position_array[2][j])*(position_array[2][i] - position_array[2][j]);
            diff[i][j] = sqrt(temp_var*temp_var*temp_var); /*sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)*/
        }
    }
    
    for (i = 0; i < 3; i++) /*compute k-value of velocity base on Law of Universal Gravitation �����ٶȵ�kֵ������ţ��������������*/
    {
        temp_k[i][0] = -(mc + mp[0])*position_array[i][0] / dst[0] - mp[1] * ((position_array[i][0] - position_array[i][1]) / diff[0][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][0] - position_array[i][2]) / diff[0][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][0] - position_array[i][3]) / diff[0][3] + position_array[i][3] / dst[3])- mp[4] * ((position_array[i][0] - position_array[i][4]) / diff[0][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][0] - position_array[i][5]) / diff[0][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][0] - position_array[i][6]) / diff[0][6] + position_array[i][6] / dst[6]);
        temp_k[i][1] = -(mc + mp[1])*position_array[i][1] / dst[1] - mp[0] * ((position_array[i][1] - position_array[i][0]) / diff[1][0] + position_array[i][0] / dst[0])- mp[2] * ((position_array[i][1] - position_array[i][2]) / diff[1][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][1] - position_array[i][3]) / diff[1][3] + position_array[i][3] / dst[3])- mp[4] * ((position_array[i][1] - position_array[i][4]) / diff[1][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][1] - position_array[i][5]) / diff[1][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][1] - position_array[i][6]) / diff[1][6] + position_array[i][6] / dst[6]);
        temp_k[i][2] = -(mc + mp[2])*position_array[i][2] / dst[2] - mp[0] * ((position_array[i][2] - position_array[i][0]) / diff[2][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][2] - position_array[i][1]) / diff[2][1] + position_array[i][1] / dst[1])- mp[3] * ((position_array[i][2] - position_array[i][3]) / diff[2][3] + position_array[i][3] / dst[3])- mp[4] * ((position_array[i][2] - position_array[i][4]) / diff[2][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][2] - position_array[i][5]) / diff[2][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][2] - position_array[i][6]) / diff[2][6] + position_array[i][6] / dst[6]);
        temp_k[i][3] = -(mc + mp[3])*position_array[i][3] / dst[3] - mp[0] * ((position_array[i][3] - position_array[i][0]) / diff[3][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][3] - position_array[i][1]) / diff[3][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][3] - position_array[i][2]) / diff[3][2] + position_array[i][2] / dst[2])- mp[4] * ((position_array[i][3] - position_array[i][4]) / diff[3][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][3] - position_array[i][5]) / diff[3][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][3] - position_array[i][6]) / diff[3][6] + position_array[i][6] / dst[6]);
        temp_k[i][4] = -(mc + mp[4])*position_array[i][4] / dst[4] - mp[0] * ((position_array[i][4] - position_array[i][0]) / diff[4][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][4] - position_array[i][1]) / diff[4][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][4] - position_array[i][2]) / diff[4][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][4] - position_array[i][3]) / diff[4][3] + position_array[i][3] / dst[3]) - mp[5] * ((position_array[i][4] - position_array[i][5]) / diff[4][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][4] - position_array[i][6]) / diff[4][6] + position_array[i][6] / dst[6]);
        temp_k[i][5] = -(mc + mp[5])*position_array[i][5] / dst[5] - mp[0] * ((position_array[i][5] - position_array[i][0]) / diff[5][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][5] - position_array[i][1]) / diff[5][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][5] - position_array[i][2]) / diff[5][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][5] - position_array[i][3]) / diff[5][3] + position_array[i][3] / dst[3]) - mp[4] * ((position_array[i][5] - position_array[i][4]) / diff[5][4] + position_array[i][4] / dst[4]) - mp[6] * ((position_array[i][5] - position_array[i][6]) / diff[5][6] + position_array[i][6] / dst[6]);
        temp_k[i][6] = -(mc + mp[6])*position_array[i][6] / dst[6] - mp[0] * ((position_array[i][6] - position_array[i][0]) / diff[6][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][6] - position_array[i][1]) / diff[6][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][6] - position_array[i][2]) / diff[6][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][6] - position_array[i][3]) / diff[6][3] + position_array[i][3] / dst[3]) - mp[4] * ((position_array[i][6] - position_array[i][4]) / diff[6][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][6] - position_array[i][5]) / diff[6][5] + position_array[i][5] / dst[5]);
        
        /*-----begin----for two planets----------*/
        /*
         * temp_k[i][0] = -(mc + mp[0])*position_array[i][0] / dst[0] - mp[1] * ((position_array[i][0] - position_array[i][1]) / diff[0][1] + position_array[i][1] / dst[1]);
         * temp_k[i][1] = -(mc + mp[1])*position_array[i][1] / dst[1] - mp[0] * ((position_array[i][1] - position_array[i][0]) / diff[1][0] + position_array[i][0] / dst[0]);
         */
        /*-----end----for two planets----------*/
        
    }
    
}