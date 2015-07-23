
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;

/**
 *
 * @author ngo
 */
public class WeatherPredict {

    /**
     * @param args the command line arguments
     */
    private static final int X=50;
    private static final int Y=40;
    private static final int K=5;
    private static final int R=287;
    private static final int CP=1005;
    private static final int M=1;
    private static final int N=50;
    private static final int[] PN= {1000, 900, 650, 400, 150};
    private static final int P0=1000;
    private static final float del=1;
//    float[][][] U_t= new float[X][Y][K];
//    float[][][] U_t1= new float[X][Y][K];
//    float[][][] U_t2= new float[X][Y][K];
    static float[][][][] U= new float[N][X][Y][K];
//    float[][][] V_t= new float[X][Y][K];
//    float[][][] V_t1= new float[X][Y][K];
//    float[][][] V_t2= new float[X][Y][K];
    static float[][][][] V= new float[N][X][Y][K];
//    float[][][] TE= new float[X][Y][K];
//    float[][][] TE1= new float[X][Y][K];
//    float[][][] TE2= new float[X][Y][K];
    static float[][][][] TE= new float[N][X][Y][K];
//    float[][][] Q_t= new float[X][Y][K];
//    float[][][] Q_t1= new float[X][Y][K];
//    float[][][] Q_t2= new float[X][Y][K];
    static float[][][][] Q= new float[N][X][Y][K];
//    float[][][] PS= new float[X][Y][K];
//    float[][][] PS1= new float[X][Y][K];
//    float[][][] PS2= new float[X][Y][K];
    static float[][][] PS= new float[N][X][Y];
    static float[][][] Wt= new float[X][Y][K];
    static float[][][] Ws= new float[X][Y][1];
    static float[][][] PHI= new float[X][Y][K];
    static float F;
    public void init(){
        
    }
    public static void main(String[] args) {
        // TODO code application logic here
        float M1= M/(4*del);
        float RC= R/CP;
        float P1, P2, P3, P4; 
        float W9;
        float R0, TV,V_TB;
        for(int n=0; n<N; n++){
            for(int k=0; k<K; k++){
                //giai doan du doan
                for (int i=1;i<X-1; i++){
                    for(int j=1; j<Y-1; j++){
                        U[n+2][i][j][k]=-M*(PHI[i][j][k]-PHI[i+1][j][k])/del
                            -M1*((U[n+1][i][j][k]+U[n+1][i+1][j][k])*(U[n+1][i+1][j][k]-U[n+1][i][j][k])+(U[n+1][i][j][k]+U[n+1][i-1][j][k])*(U[n+1][i][j][k]-U[n+1][i-1][j][k]))
                            +F*(V[n+1][i][j][k]+V[n+1][i][j-1][k]+V[n+1][i+1][j][k]+V[n+1][i+1][j-1][k])/4
                            -M1*((V[n+1][i][j][k]+V[n+1][i+1][j][k])*(U[n+1][i][j+1][k]-U[n+1][i][j][k])+(V[n+1][i][j-1][k]+V[n+1][i+1][j-1][k])*(U[n+1][i][j-i][k]-U[n+1][i][j][k]));
                        V[n+2][i][j][k]= -M*(PHI[i][j+1][k]-PHI[i][j][k])/del
                            - M1*((U[n+1][i][j][k]+U[n+1][i][j+1][k])*(V[n+1][i+1][j][k]-V[n+1][i][j][k])+(U[n+1][i-1][j][k]+U[n+1][i-1][j+1][k])*(V[n+1][i][j][k]-V[n+1][i][j][k]))
                            - M1*((V[n+1][i][j][k]+V[n+1][i][j+1][k])*(V[n+1][i][j+1][k]-V[n+1][i][j][k])+(V[n+1][i][j][k]+V[n+1][i][j-1][k])*(V[n+1][i][j][k]-V[n+1][i][j-1][k]))
                            -F*(U[n+1][i][j][k]+U[n+1][i-1][j][k]+U[n+1][i][j+1][k]+V[n+1][i-1][j+1][k])/4;
                        TE[n+2][i][j][k]=-2*M1*(U[n+1][i+1][j][k]*(TE[n+1][i+1][j][k]-TE[n+1][i][j][k])
                                -U[n+1][i][j][k]*(TE[n+1][i][j][k]-TE[n+1][i-1][j][k])
                                +V[n+1][i][j+1][k]*(TE[n+1][i][j+1][k]-TE[n+1][i][j][k])
                                -V[n+1][i][j][k]*(TE[n+1][i][j][k]-TE[n+1][i][j-1][k]));
                        Q[n+2][i][j][k]=-2*M1*(U[n+1][i][j][k]*(Q[n+1][i+1][j][k]-Q[n+1][i][j][k])
                                -U[n+1][i][j][k]*(Q[n+1][i][j][k]-Q[n+1][i-1][j][k])
                                +V[n-1][i][j+1][k]*(Q[n+1][i][j+1][k]-Q[n+1][i][j][k])
                                -V[n+1][i][j][k]*(Q[n+1][i][j][k]-Q[n+1][i][j-1][k]));
                        if(k==0){
                            P1= (PN[k+1]+PN[k])/2;                            
                            P2= (PN[k+2]+PN[k+1])/2;                                                        
                            //**1
                            U[n+2][i][j][k]=U[n+2][i][j][k]-(Wt[i+1][j][k]+Wt[i][j][k])*(U[n+1][i][j][k+1]-U[n+1][i][j][k])/(2*(PN[k+1]-PN[k]));
                            //U[n+1][i][j][k]=U[n+1][i][j][k]-(Wt[i+1][j][k]+Wt[i][j][k])*(U[n][i][j][k+1]-U[n][i][j][k])/(2*PN[k+1]-PN[k]);
                            //**2
                            V[n+2][i][j][k]=V[n+2][i][j][k]-(Wt[i][j+1][k]+Wt[i][j][k])*(V[n+1][i][j][k+1]-V[n+1][i][j][k])/(2*(PN[k+1]-PN[k]));
                            //V[n+1][i][j][k]=V[n+1][i][j][k]-(Wt[i][j+1][k]+Wt[i][j][k])*(V[n][i][j][k+1]-V[n][i][j][k])/(2*PN[k+1]-PN[k]);
                            //**3
                            TE[n+2][i][j][k]=TE[n+2][i][j][k]-Wt[i][j][k]*(TE[n+1][i][j][k+1]-TE[n+1][i][j][k])/(PN[k+1]-PN[k]);
//                            TE[n+1][i][j][k]=TE[n+1][i][j][k]-Wt[i][j][k]*(TE[n][i][j][k+1]-TE[n][i][j][k])/(PN[k+1]-PN[k]);
                            //**4
                            Q[n+2][i][j][k]= Q[n+2][i][j][k]-Wt[i][j][k]*(Q[n+1][i][j][k+1]-Q[n+1][i][j][k])/(PN[k+1]-PN[k]);
//                            Q[n+1][i][j][k]= Q[n+1][i][j][k]-Wt[i][j][k]*(Q[n][i][j][k+1]-Q[n][i][j][k])/(PN[k+1]-PN[k]);
                            //**5
                            PS[n+2][i][j]=-2*M1*(U[n+1][i+1][j][k]*(PS[n+1][i+1][j]-PS[n+1][i][j])
                                    -U[n+1][i][j][k]*(PS[n+1][i][j]-PS[n+1][i-1][j])
                                    +V[n+1][i][j+1][k]*(PS[n+1][i][j+1]-PS[n+1][i][j]));
//                            PS[n+1][i][j]=-2*M1*(U[n][i+1][j][k]*(PS[n][i+1][j]-PS[n][i][j])
//                                    -U[n][i][j][k]*(PS[n][i][j]-PS[n][i-1][j])
//                                    +V[n][i][j+1][k]*(PS[n][i][j+1]-PS[n][i][j]));
                                    //-);
                            W9= (Wt[i][j][k+2]*P1+Wt[i][j][k]*P2)/(P1+P2);
                            PS[n+2][i][j]=PS[n+2][i][j]+W9-(M*((U[n+1][i+1][j][k+1]-U[n+1][i][j][k])
                                    +(V[n+1][i][j+1][k+1]-V[n+1][i][j][k+1])+(U[n+1][i+1][j][k]-U[n+1][i][j][k])
                                    +(V[n+1][i][j+1][k]-V[n+1][i][j][k]))/(2*del))*(PS[n+1][i][j]-900);
                        }else if(k==K-1){
                            //***
                            P1= (PN[k]-PN[k-1])/2;                            
                            P2=25;                            
                            //***1
                            U[n+2][i][j][k]= U[n+2][i][j][k]
                                    -(Wt[i+1][j][k-1]+Wt[i][j][k-1])*(PN[k]-P2)/((2*(P1- P2)))*((U[n+1][i][j][k]-U[n+1][i][j][k-1])/(PN[k]-PN[k-1]));
//                            U[n][i][j][k]= U[n][i][j][k]
//                                    -(Wt[i+1][j][k-1]+Wt[i][j][k-1])*(PN[k]-P2)/((2*(P1- P2)))*((U[n][i][j][k]-U[n][i][j][k-1])/(PN[k]-PN[k-1]));
                            //***2
                            V[n+2][i][j][k]= V[n+2][i][j][k]
                                    -((Wt[i][j+1][k]+Wt[i][j][k-1])*(PN[k]-PN[k-1])/(2*(P1-P2)))*((V[n+1][i][j][k-1]-V[n+1][i][j][k-1])/(PN[k]-PN[k-1]));
//                            V[n+1][i][j][k]= V[n+1][i][j][k]
//                                    -((Wt[i][j+1][k]+Wt[i][j][k-1])*(PN[k]-PN[k-1])/(2*(P1-P2)))*((V[n][i][j][k-1]-V[n][i][j][k-1])/(PN[k]-PN[k-1]));
                            //***3
                            TE[n+2][i][j][k]=0;
//                            TE[n+1][i][j][k]=0;
                            Q[n+2][i][j][k]=0;
//                            Q[n+1][i][j][k]=0;
                        }else{
                        P1= (PN[k-1]+PN[k])/2;
                        P2= (PN[k]+PN[k+1])/2;
                        
                        U[n+2][i][j][k]= U[n+2][i][j][k]-(((Wt[i+1][j][k]+Wt[i][j][k])*(PN[k]-P2)
                                +(Wt[i+1][j][k+1]+Wt[i][j][k+1])*(P1-PN[k]))/(2*(P1-P2))
                                *((U[n+1][i][j][k-1]-U[n+1][i][j][k])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                +((U[n+1][i][j][k]-U[n+1][i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1])))/(PN[k-1]-PN[k+1]);
//                        U[n+1][i][j][k]=-M*(PHI[i][j][k]-PHI[i+1][j][k])/del
//                            -M1*((U[n][i][j][k]+U[n][i+1][j][k])*(U[n][i+1][j][k]-U[n][i][j][k])+(U[n][i][j][k]+U[n][i-1][j][k])*(U[n][i][j][k]-U[n][i-1][j][k]))
//                            +F*(V[n][i][j][k]+V[n][i][j-1][k]+V[n][i+1][j][k]+V[n][i+1][j-1][k])/4
//                            -M1*((V[n][i][j][k]+V[n][i+1][j][k])*(U[n][i][j+1][k]-U[n][i][j][k])+(V[n][i][j-1][k]+V[n][i+1][j-1][k])*(U[n][i][j-i][k]-U[n][i][j][k]));
                        
                        V[n+2][i][j][k]= V[n+2][i][j][k]-(((Wt[i][j+1][k]+Wt[i][j][k])*(PN[k]-P2)
                                +(Wt[i][j+1][k+1]+Wt[i][j][k+1])*(P1-PN[k]))/(2*(P1-P2))
                                *((V[n+1][i][j][k-1]-V[n+1][i][j][k+1])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                +((V[n+1][i][j][k]-V[n+1][i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1])))/(PN[k-1]-PN[k+1]);
                                
//                        V[n+1][i][j][k]= -M*(PHI[i][j+1][k]-PHI[i][j][k])/del
//                            - M1*((U[n][i][j][k]+U[n][i][j+1][k])*(V[n][i+1][j][k]-V[n][i][j][k])+(U[n][i-1][j][k]+U[n][i-1][j+1][k])*(V[n][i][j][k]-V[n][i][j][k]))
//                            - M1*((V[n][i][j][k]+V[n][i][j+1][k])*(V[n][i][j+1][k]-V[n][i][j][k])+(V[n][i][j][k]+V[n][i][j-1][k])*(V[n][i][j][k]-V[n][i][j-1][k]))
//                            -F*(U[n][i][j][k]+U[n][i-1][j][k]+U[n][i][j+1][k]+V[n][i-1][j+1][k])/4;
                        
                        TE[n+2][i][j][k]=TE[n+2][i][j][k]-Wt[i][j][k]*((TE[n+1][i][j][k-1]-TE[n+1][i][j][k])
                                *(PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(TE[n+1][i][j][k]-TE[n][i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
//                        TE[n+1][i][j][k]=-2*M1*(U[n][i+1][j][k]*(TE[n][i+1][j][k]-TE[n][i][j][k])
//                                -U[n][i][j][k]*(TE[n][i][j][k]-TE[n][i-1][j][k])
//                                +V[n][i][j+1][k]*(TE[n][i][j+1][k]-TE[n][i][j][k])
//                                -V[n][i][j][k]*(TE[n][i][j][k]-TE[n][i][j-1][k]));
                        
                        Q[n+2][i][j][k]= Q[n+2][i][j][k]-Wt[i][j][k]*((Q[n+1][i][j][k-1]-Q[n+1][i][j][k])*
                                (PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(Q[n+1][i][j][k]-Q[n+1][i][j][k+1])*
                                (PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
//                        Q[n+1][i][j][k]=-2*M1*(U[n][i][j][k]*(Q[n][i+1][j][k]-Q[n][i][j][k])
//                                -U[n][i][j][k]*(Q[n][i][j][k]-Q[n][i-1][j][k])
//                                +V[n][i][j+1][k]*(Q[n][i][j+1][k]-Q[n][i][j][k])
//                                -V[n][i][j][k]*(Q[n][i][j][k]-Q[n][i][j-1][k]));
                        }
                    }
                }
                //giai doan chuan doan
                for (int i=1; i<X-1; i++){
                    for(int j=1; j<Y-1; j++){
                        R0= (float) (TE[n+1][i][j][k]*Math.pow(PN[k]/P0, RC));
                        TV=(float) (TE[n+1][i][j][k]*(1+0.61*Q[n+1][i][j][k]));
                        if(k==0){
                            P1=(PN[k]-PN[k+1])/2;
                            P2=(PN[k+1]-PN[k+2])/2;
                            P3=(PN[k+1]+PN[k+2])/2;
                            P4=(PN[1]+PN[2])/2;
                            //**
                            V_TB=-M*((U[n+1][i+1][j][k]-U[n+1][i][j][k])
                                 +(V[n+1][i][j+1][k]-V[n+1][i][j][k])*P2
                                 +(U[n+1][i+1][j][k+1]-U[n+1][i][j][k+1])
                                 +(V[n+1][i][j+1][k]-V[n+1][i][j][k+1])*P1)/(del*(P1+P2));
                            float VTS= M*((U[n+1][i+1][j][k]-U[n+1][i][j][k])+(V[n+1][i][j+1][k]-V[n+1][i][j][k]))/del;
                            Ws[i][j][1]=Ws[i][j][1]+VTS*(P4- PN[k]);
                            //**1
                            Wt[i][j][k]=(-V_TB*(PN[k]-P3)+(Ws[i][j][1]*((PN[k+1]-P3)/(PN[k]-PN[k+1])))
                                    -(Wt[i][j][k+1]*((PN[k]-PN[k+1])/(PN[k+1]-P3))))/((PN[k+1]-P3)/(PN[k]-PN[k+1])-(PN[k]-PN[k+1])/(PN[k+1]-P3));
                            //**2
                            PHI[i][j][k]=(float) (PHI[i][j][k+1]+R/R0*TV*Math.pow((PN[k]/P0), RC)*(PN[k+1]-PN[k]));
                        } else if (k==K-1){
                            //***1
                            PHI[i][j][k]= (float) (PHI[i][j][k-1]+R/R0*TV*Math.pow(PN[k]/P0, RC)*(PN[k]-PN[k-1]));
                            Wt[i][j][k]=0;
                        }else{
                            PHI[i][j][k]=(float) ((R/R0*TV*Math.pow(PN[k]/P0, RC)*(PN[k-1]-PN[k+1])+PHI[i][j][k+1]*((PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                    -PHI[i][j][k+1]*((PN[k-1]-PN[k])/(PN[k]-PN[k-1])))/((PN[k]-PN[k+1])/(PN[k-1]-PN[k])-((PN[k-1]-PN[k])/(PN[k]-PN[k+1]))));
                            V_TB=-M*((U[n+1][i+1][j][k]-U[n+1][i][j][k])
                                    +(V[n+1][i][j+1][k]-V[n+1][i][j][k])
                                    +(U[n+1][i+1][j][k+1]-U[n+1][i][j][k+1])
                                    +(V[n+1][i][j+1][k+1]-V[n+1][i][j][k+1]))/(2*del);
                            Wt[i][j][k]=(-V_TB*(PN[k-1]-PN[k+1])
                                    +Wt[i][j][k-1]*((PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                    -Wt[i][j][k+1]*((PN[k-1]-PN[k])/(PN[k]-PN[k+1])))
                                    /((PN[k]-PN[k+1])/(PN[k-1]-PN[k])-(PN[k-1]-PN[k])/(PN[k]-PN[k+1]));
                        }
                        
                    }
                }
            }
        }
    }
}
