/*
This is the function for searching transit time using RKF45.
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
/*#include "stdio.h"
 */
#include"string.h"

#define N_transit 500 /**/
#define N_planet  7  /*行星数量*/
#define N_para 35  /*参数个数*/
#define N_order 6  /*与RKF算法阶数有关,RKF78时取值13，RKF45时取值6*/


const double a[N_order+1][N_order-1] = {{0,0,0,0,0},{0.250000000000000,0,0,0,0},{0.0937500000000000,0.281250000000000,0,0,0},{0.879380974055530,-3.27719617660446,3.32089212562585,0,0},{2.03240740740741,-8,7.17348927875244,-0.205896686159844,0},{-0.296296296296296,2,-1.38167641325536,0.452972709551657,-0.275000000000000},{0,0,0,0,0}};
const double pi = 3.14159265358979323846264338328;


void k_v(double position_array[3][N_planet], double mc, double mp[N_planet], double temp_k[3][N_planet]);
void rkf45(double *h, double *t, double position_array[3][N_planet], double velocity_array[3][N_planet], double mc, double mp[N_planet],double err);
void transit_time(double *tr, double *error_flag, double p_para[N_para], double mc, double Rs, double tstop, double time_precise,double err);
/*
 * int main()
 * {...
 * return 0;
 * }
 */


/*-------------matlab API function bellow------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        /*input parameter*/
    double *p_para;
    double *mc;
    double *Rs;
    double *tstop;
    double *time_precise;
    double *err;
    
        /*outuput parameter*/
    double *error_flag;
    double *tr;
    
    
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
    
    
    plhs[0] = mxCreateDoubleMatrix(N_transit, N_planet, mxREAL);
    tr = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    error_flag = mxGetPr(plhs[1]);
    
    transit_time(tr, error_flag, p_para, *mc, *Rs, *tstop, *time_precise,*err);
    
}
/*-------------matlab API function above------------------------*/




void transit_time(double *tr, double *error_flag, double p_para[N_para], double mc, double Rs, double tstop, double time_precise,double err)
{
    double t = 0.0;
    double h = 0.02 * pi;
    double t1 = 0.0;
    double h1 = 0.0;
    double p[N_planet], e[N_planet], ci[N_planet], w[N_planet], cm0[N_planet], mp[N_planet], omega[N_planet];
    double ap[N_planet], cn[N_planet], ce[N_planet], d[N_planet], dd[N_planet];
    double position_array[3][N_planet], velocity_array[3][N_planet], temp_poistion_1[3][N_planet], temp_velocity_1[3][N_planet];
    double px, py, pz, qx, qy, qz, r;
    double ce0, gap;
    int num[N_planet];
    int flag[N_planet];
    int i, j,loop_index;
    
    /*
     * for (j = 0; j < N_planet * 3; j++)
     * {
     * temp_poistion_1[j % 3][j / 3] =0.0;
     * temp_velocity_1[j % 3][j / 3] = 0.0;
     * position_array[j % 3][j / 3] =0.0;
     * velocity_array[j % 3][j / 3] = 0.0;
     * }
     *
     *
     */
    
    
    memset(position_array, 0, 3 * N_planet);
    memset(velocity_array, 0, 3 * N_planet);
    memset(temp_poistion_1, 0, 3 * N_planet);
    memset(temp_velocity_1, 0, 3 * N_planet);
    
    for (i = 0; i < N_planet; i++)
    {
        num[i] = 0;
        flag[i] = 0;
    }
    
    
    
    /*---下面的要改：-不同的情况(比如行星数不一样，拟合参数不一样)-- 参数顺序可能不一样：p, e, w, cm, mp, ci, omega --*/
    for (i = 0; i < N_planet; i++)
    {
        p[i] = p_para[i];
        e[i] = p_para[i + N_planet];
        w[i] = p_para[i + 2 * N_planet];
        cm0[i] = p_para[i + 3 * N_planet];
        mp[i] = p_para[i + 4 * N_planet];
        
        ci[i] = pi*0.5;
        omega[i] =0.0;
    /*--上面的要改：不同的情况(比如行星数不一样，拟合参数不一样)------参数顺序可能不一样----*/
        
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
    
    *error_flag=0;
    
    while (t < tstop)
    {

        if (h<time_precise)
        {
         *error_flag=1;
         break;
        }

        
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
                
                /*
                 * for (j = 0; j < N_planet * 3; j++)
                 * {
                 * temp_poistion_1[j % 3][j / 3] = position_array[j % 3][j / 3];
                 * temp_velocity_1[j % 3][j / 3] = velocity_array[j % 3][j / 3];
                 * }
                 */
                
                
                dd[i] = d[i];
                while (gap > err&&d[i] > err)
                {
                    if (d[i] > dd[i])
                    {
                        h1 = -0.5*h1;
                    }
                    dd[i] = d[i];
                    rkf45(&h1, &t1, temp_poistion_1, temp_velocity_1, mc, mp,err);
                    d[i] = sqrt(temp_poistion_1[0][i] * temp_poistion_1[0][i] + temp_poistion_1[1][i] * temp_poistion_1[1][i]);/*sqrt(x^2+y^2)*/
                    gap = fabs(dd[i] - d[i]);
                }
                
                
                
                
                *(tr + i * (N_transit)+num[i]) = (t1 - h1) / 2.0 / pi * 365.2422;   /*--        *(tr + i * (N_transit)+num[i]) = (t1 - h1) *58.130101555758664; --*/
                flag[i] = 1;
                num[i]+=1;
                
            }
            
            if (position_array[2][i] < 0 && flag[i]==1)
            {
                flag[i] = 0;
            }
        }
        
        rkf45(&h, &t, position_array, velocity_array, mc, mp,err);
    }
    
}

void rkf45(double *h, double *t, double position_array[3][N_planet], double velocity_array[3][N_planet], double mc, double mp[N_planet],double err)
{
    int j;
    int loop_index;
    int K_index;
    double k_velocity[3][N_planet];
    double K_po[N_order][N_planet * 3];
    double K_ve[N_order][N_planet * 3];
    double temp_poistion[3][N_planet];
    double temp_velocity[3][N_planet];
    double derr = 1.0;
    double delta_po, delta_ve, max_po, max_ve;
    double err_po[3 * N_planet];
    double err_ve[3 * N_planet];
    
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
    
    
    memset(k_velocity, 0, 3 * N_planet);
    memset(K_po, 0, N_order * 3 * N_planet);
    memset(K_ve, 0, N_order * 3 * N_planet);
    memset(temp_poistion, 0, 3 * N_planet);
    memset(temp_velocity, 0, 3 * N_planet);
    memset(err_po, 0, 3 * N_planet);
    memset(err_ve, 0, 3 * N_planet);
    
    
    
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
                K_po[K_index][j] = temp_velocity[j % 3][j / 3];
            }
            
            k_v(temp_poistion, mc, mp, k_velocity);
            
            for (j = 0; j < N_planet * 3; j++)
            {
                K_ve[K_index][j] = k_velocity[j % 3][j / 3];
            }
            
            if (K_index<(N_order-1))
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
        
        max_po=0;
        max_ve=0;
        for (j = 0; j < N_planet * 3; j++)
        {
            err_po[j] = (*h) *fabs((16.0/135.0-25.0/216.0)*K_po[0][j] + (6656.0/12825.0-1408.0/2565.0)*K_po[2][j]+(28561.0/56430.0-2197.0/4104.0)*K_po[3][j] + (-9.0/50.0+0.2)*K_po[4][j]+2.0/55.0*K_po[5][j]);
            err_ve[j] = (*h) *fabs((16.0/135.0-25.0/216.0)*K_ve[0][j] + (6656.0/12825.0-1408.0/2565.0)*K_ve[2][j]+(28561.0/56430.0-2197.0/4104.0)*K_ve[3][j] + (-9.0/50.0+0.2)*K_ve[4][j]+2.0/55.0*K_ve[5][j]);
            if (err_po[j]>max_po)
            {max_po=err_po[j];}
            if (err_ve[j]>max_ve)
            {max_ve=err_ve[j];}
        }
        
        derr = max_ve;
        if (max_po > max_ve)
        {
            derr = max_po;
        }
        
        
        
        if (derr > err)
        {
            (*h) = (*h) *0.5;
        }
        
    }
    
    
    (*t) +=(*h);
    for (j = 0; j < N_planet * 3; j++)
    {
    	position_array[j % 3][j / 3] += (*h) * (16.0/135.0 * K_po[0][j] + 6656.0/12825.0 * K_po[2][j] +28561.0/56430.0* K_po[3][j]- 9.0/50.0 * K_po[4][j] + 2.0/55.0 * K_po[5][j]);
    	velocity_array[j % 3][j / 3] += (*h) * (16.0/135.0 * K_ve[0][j] + 6656.0/12825.0 * K_ve[2][j] +28561.0/56430.0* K_ve[3][j]- 9.0/50.0 * K_ve[4][j] + 2.0/55.0 * K_ve[5][j]);        
    }
    
    
}




void k_v(double position_array[3][N_planet], double mc, double mp[N_planet], double temp_k[3][N_planet])
{
    int i;
    int j;
    double dst[N_planet];
    double diff[N_planet][N_planet];
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
    
    memset(dst, 0, N_planet);
    memset(diff, 0, N_planet*N_planet);
    
    
    
    for (i = 0; i < N_planet; i++)
    {
        temp_var = (position_array[0][i] * position_array[0][i] + position_array[1][i] * position_array[1][i] + position_array[2][i] * position_array[2][i]);
        dst[i] = sqrt(temp_var*temp_var*temp_var);
        
        for (j = 0; j < N_planet; j++)
        {
            
            temp_var = (position_array[0][i] - position_array[0][j])*(position_array[0][i] - position_array[0][j]) + (position_array[1][i] - position_array[1][j])*(position_array[1][i] - position_array[1][j]) + (position_array[2][i] - position_array[2][j])*(position_array[2][i] - position_array[2][j]);
            diff[i][j] = sqrt(temp_var*temp_var*temp_var);
        }
    }
    
    for (i = 0; i < 3; i++)
    {
        temp_k[i][0] = -(mc + mp[0])*position_array[i][0] / dst[0] - mp[1] * ((position_array[i][0] - position_array[i][1]) / diff[0][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][0] - position_array[i][2]) / diff[0][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][0] - position_array[i][3]) / diff[0][3] + position_array[i][3] / dst[3])- mp[4] * ((position_array[i][0] - position_array[i][4]) / diff[0][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][0] - position_array[i][5]) / diff[0][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][0] - position_array[i][6]) / diff[0][6] + position_array[i][6] / dst[6]);
        temp_k[i][1] = -(mc + mp[1])*position_array[i][1] / dst[1] - mp[0] * ((position_array[i][1] - position_array[i][0]) / diff[1][0] + position_array[i][0] / dst[0])- mp[2] * ((position_array[i][1] - position_array[i][2]) / diff[1][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][1] - position_array[i][3]) / diff[1][3] + position_array[i][3] / dst[3])- mp[4] * ((position_array[i][1] - position_array[i][4]) / diff[1][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][1] - position_array[i][5]) / diff[1][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][1] - position_array[i][6]) / diff[1][6] + position_array[i][6] / dst[6]);
        temp_k[i][2] = -(mc + mp[2])*position_array[i][2] / dst[2] - mp[0] * ((position_array[i][2] - position_array[i][0]) / diff[2][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][2] - position_array[i][1]) / diff[2][1] + position_array[i][1] / dst[1])- mp[3] * ((position_array[i][2] - position_array[i][3]) / diff[2][3] + position_array[i][3] / dst[3])- mp[4] * ((position_array[i][2] - position_array[i][4]) / diff[2][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][2] - position_array[i][5]) / diff[2][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][2] - position_array[i][6]) / diff[2][6] + position_array[i][6] / dst[6]);
        temp_k[i][3] = -(mc + mp[3])*position_array[i][3] / dst[3] - mp[0] * ((position_array[i][3] - position_array[i][0]) / diff[3][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][3] - position_array[i][1]) / diff[3][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][3] - position_array[i][2]) / diff[3][2] + position_array[i][2] / dst[2])- mp[4] * ((position_array[i][3] - position_array[i][4]) / diff[3][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][3] - position_array[i][5]) / diff[3][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][3] - position_array[i][6]) / diff[3][6] + position_array[i][6] / dst[6]);
        temp_k[i][4] = -(mc + mp[4])*position_array[i][4] / dst[4] - mp[0] * ((position_array[i][4] - position_array[i][0]) / diff[4][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][4] - position_array[i][1]) / diff[4][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][4] - position_array[i][2]) / diff[4][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][4] - position_array[i][3]) / diff[4][3] + position_array[i][3] / dst[3]) - mp[5] * ((position_array[i][4] - position_array[i][5]) / diff[4][5] + position_array[i][5] / dst[5]) - mp[6] * ((position_array[i][4] - position_array[i][6]) / diff[4][6] + position_array[i][6] / dst[6]);
        temp_k[i][5] = -(mc + mp[5])*position_array[i][5] / dst[5] - mp[0] * ((position_array[i][5] - position_array[i][0]) / diff[5][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][5] - position_array[i][1]) / diff[5][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][5] - position_array[i][2]) / diff[5][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][5] - position_array[i][3]) / diff[5][3] + position_array[i][3] / dst[3]) - mp[4] * ((position_array[i][5] - position_array[i][4]) / diff[5][4] + position_array[i][4] / dst[4]) - mp[6] * ((position_array[i][5] - position_array[i][6]) / diff[5][6] + position_array[i][6] / dst[6]);
        temp_k[i][6] = -(mc + mp[6])*position_array[i][6] / dst[6] - mp[0] * ((position_array[i][6] - position_array[i][0]) / diff[6][0] + position_array[i][0] / dst[0])- mp[1] * ((position_array[i][6] - position_array[i][1]) / diff[6][1] + position_array[i][1] / dst[1])- mp[2] * ((position_array[i][6] - position_array[i][2]) / diff[6][2] + position_array[i][2] / dst[2])- mp[3] * ((position_array[i][6] - position_array[i][3]) / diff[6][3] + position_array[i][3] / dst[3]) - mp[4] * ((position_array[i][6] - position_array[i][4]) / diff[6][4] + position_array[i][4] / dst[4]) - mp[5] * ((position_array[i][6] - position_array[i][5]) / diff[6][5] + position_array[i][5] / dst[5]);
        
        /*---------下面是对于2颗行星的----------*/        
        /*
         * temp_k[i][0] = -(mc + mp[0])*position_array[i][0] / dst[0] - mp[1] * ((position_array[i][0] - position_array[i][1]) / diff[0][1] + position_array[i][1] / dst[1]);
         * temp_k[i][1] = -(mc + mp[1])*position_array[i][1] / dst[1] - mp[0] * ((position_array[i][1] - position_array[i][0]) / diff[1][0] + position_array[i][0] / dst[0]);
         */
        /*---------上面是对于2颗行星的----------*/
    }
    
}
